# Minigraph-Cactus pangenome construction and structural variant detection pipeline for _P. brassicae_


## 1. Pangenome graph construction with Minigraph-Cactus

```bash

## Running minigraph cactus (mount the directory with data to run minigraph cactus)

docker run -it -v /path/to/directory/minigraph_cactus:/data quay.io/comparative-genomics-toolkit/cactus:v2.9.9 bash

cactus-pangenome --help

cactus-pangenome /data/jobStore /data/seqfile.tsv --outDir /data/pb_pangenome --outName pb_pangenome --reference Pb3A --vcf --giraffe --gfa --gbz --odgi --chrom-vg --chrom-og --viz --draw --mgCores 12

# To rerun the command

toil clean /data/jobStore

## Make additional gfa by adding paths with VG

vg convert -fW pb_pangenome.gbz > pb_pangenome.paths.gfa

cut -f1 pb_pangenome.gfa | sort | uniq -c # checking paths is present

```

## 2. Structural variant VCF post-processing and SV-TYPE annotation

```bash

# 1) Flatten nested graph sites (vcfbub)

vcfbub --input pb_pangenome.vcf.gz --max-level 0 --max-ref-length 10000 > pb_45.flattened.vcf

# Make a proper BGZF-compressed, sorted copy
bcftools sort pb_pangenome.flattened.vcf.gz -Oz -o pb_pangenome.flattened.sorted.vcf.gz

# Index it (tbi index)
bcftools index -t pb_pangenome.flattened.sorted.vcf.gz

# Set reference paths
REF="path/to/Pb3A_renamed.fa"

# 2) Decompose/realign alleles with inversion detection (shorter k-mer & smaller min len makes INV detection more sensitive)

vcfwave pb_pangenome.flattened.vcf > pb_pangenome.flattened.waved.vcf  # without inversion SV sets

vcfwave --inv-kmer 15 --inv-min 40 pb_pangenome.flattened.vcf > pb_pangenome.flattened.waved.inv.vcf   # With inversion cut-off

# 3) Split multi-allelic to biallelic; bgzip & index
bcftools norm -m -any -Oz -o pb_pangenome.split.vcf.gz pb_pangenome.flattened.waved.vcf
bcftools index pb_45.split.vcf.gz

# 4) Left-normalize vs reference & recompute AC/AN/AF/NS
bcftools norm -f "$REF" -Oz -o pb_pangenome.split.norm.vcf.gz pb_45.split.vcf.gz
bcftools +fill-tags pb_pangenome.split.norm.vcf.gz -Oz -o pb_pangenome.split.norm.fill.vcf.gz -- -t AC,AN,AF,NS

## sorting and idexing

IN="pb_pangenome_200kb.flattened.sorted.waved.split.SORTED.norm.fillgaps.vcf.gz"
OUT="${IN%.vcf.gz}.SORTED.vcf.gz"

TMPDIR="/path/to/directory/tmp_bcftools_sort"
mkdir -p "$TMPDIR"

bcftools sort -T "$TMPDIR" -Oz -o "$OUT" "$IN"
bcftools index -f -c "$OUT"


# 5) Add SVTYPE/END/SVLEN + INV + (heuristic) DUP_TANDEM detection (script below)
python add_svtypes_advanced.py --vcf pb_45.split.norm.fill.vcf.gz --ref "$REF" --out pb_45.final.vcf
bgzip -f pb_45.final.vcf && bcftools index pb_45.final.vcf.gz

```

## 3. Custom SVTYPE, SVLEN and END annotation with DUP_TANDEM heuristic

```bash
#!/usr/bin/env python3

# Add SVTYPE SVLEN END and optional SVSUBTYPE to a VCF from minigraph cactus + vcfwave
# Advance SV file for all SV info tags SVTYPE/END/SVLEN (this script is advanced SV type addition with SV length including INV, INS, DEL should be ≥ 50.


Input
- VCF or VCF.GZ

Output
- Plain VCF only
  Then bgzip + sort + index with bcftools

Example
python add_sv_types_fix_inv.py --vcf in.vcf.gz --ref ref.fa --out out.vcf --sv-min 50 --inv-min 50
"""

import sys, gzip, argparse
from typing import Optional, Tuple

try:
    import pysam
except ImportError:
    sys.stderr.write("Install pysam with conda or mamba install -c bioconda pysam\n")
    sys.exit(1)

def open_in_auto(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def revcomp(s: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]

def safe_fetch(fa: pysam.FastaFile, chrom: str, start_0: int, end_0: int) -> str:
    try:
        clen = fa.get_reference_length(chrom)
    except Exception:
        return ""
    s = max(0, min(start_0, clen))
    e = max(0, min(end_0, clen))
    if e <= s:
        return ""
    try:
        return fa.fetch(chrom, s, e)
    except Exception:
        return ""

def parse_info(info: str) -> dict:
    if info == "." or not info:
        return {}
    d = {}
    for f in info.split(";"):
        if not f:
            continue
        if "=" in f:
            k, v = f.split("=", 1)
            d[k] = v
        else:
            d[f] = True
    return d

def info_to_str(d: dict) -> str:
    if not d:
        return "."
    parts = []
    for k, v in d.items():
        parts.append(k if v is True else f"{k}={v}")
    return ";".join(parts) if parts else "."

def infer_sv(
    chrom: str,
    pos1: int,
    ref: str,
    alt: str,
    ref_fa: pysam.FastaFile,
    dup_window: int,
    allele_type: Optional[str],
    inv_flag: bool,
    sv_min: int,
    inv_min: int
) -> Tuple[str, int, int, Optional[str]]:

    lr = len(ref)
    la = len(alt)
    end1 = pos1 + lr - 1 if lr > 0 else pos1
    svsub = None

    # 1) Trust vcfwave inversion flag FIRST
    # This prevents INV being swallowed by INS DEL OTHER logic
    if inv_flag and lr >= inv_min:
        return "INV", lr, end1, None

    # 2) Special alleles
    if alt in (".", "*"):
        return "OTHER", 0, end1, None
    if alt.startswith("<") or "[" in alt or "]" in alt:
        up = alt.upper().strip("<>")
        if "INV" in up and lr >= inv_min:
            return "INV", lr, end1, None
        if "DEL" in up:
            svlen = max(0, lr - la) if la > 0 else lr
            return ("DEL" if svlen >= sv_min else "OTHER"), int(svlen), end1, None
        if "INS" in up:
            return "INS", 0, pos1, None
        if "DUP" in up:
            if "TANDEM" in up:
                svsub = "DUP_TANDEM"
            return "DUP", 0, end1, svsub
        if "[" in alt or "]" in alt:
            return "BND", 0, end1, None
        return "OTHER", 0, end1, None

    # 3) Literal alleles INS DEL
    if la > lr:
        svlen = la - lr
        end1 = pos1
        svtype = "INS" if svlen >= sv_min else "OTHER"

        if svtype == "INS":
            ins_seq = None
            if lr > 0 and alt.startswith(ref):
                ins_seq = alt[lr:]
            elif lr == 1:
                ins_seq = alt[1:]

            if ins_seq and dup_window > 0:
                k = len(ins_seq)
                left = safe_fetch(ref_fa, chrom, (pos1 - 1) - k, (pos1 - 1))
                right = safe_fetch(ref_fa, chrom, pos1, pos1 + k)  # after anchor
                if ins_seq.upper() == left.upper() or ins_seq.upper() == right.upper():
                    svsub = "DUP_TANDEM"

        return svtype, svlen, end1, svsub

    if lr > la:
        svlen = lr - la
        svtype = "DEL" if svlen >= sv_min else "OTHER"
        return svtype, svlen, end1, None

    # 4) Same length replacements
    # Optional strict revcomp INV fallback if you ever want it
    if lr >= inv_min and la == lr and alt.upper() == revcomp(ref).upper():
        return "INV", lr, end1, None

    return "OTHER", 0, end1, None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dup-window", type=int, default=1000)
    ap.add_argument("--sv-min", type=int, default=50)
    ap.add_argument("--inv-min", type=int, default=50)
    args = ap.parse_args()

    try:
        ref_fa = pysam.FastaFile(args.ref)
    except Exception as e:
        sys.stderr.write(f"Cannot open reference FASTA {args.ref}: {e}\n")
        sys.exit(1)

    have_svtype = have_svlen = have_end = have_svsub = False
    hdr = []

    with open_in_auto(args.vcf) as fi, open(args.out, "wt") as fo:
        for line in fi:
            if line.startswith("##"):
                if line.startswith("##INFO=<ID=SVTYPE"):
                    have_svtype = True
                elif line.startswith("##INFO=<ID=SVLEN"):
                    have_svlen = True
                elif line.startswith("##INFO=<ID=END"):
                    have_end = True
                elif line.startswith("##INFO=<ID=SVSUBTYPE"):
                    have_svsub = True
                hdr.append(line)
                continue

            if line.startswith("#CHROM"):
                for h in hdr:
                    fo.write(h)
                if not have_svtype:
                    fo.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="INS DEL INV OTHER inferred. vcfwave INFO/INV is trusted first">\n')
                if not have_svlen:
                    fo.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="INS DEL use |len(ALT)-len(REF)|. INV uses span len(REF).">\n')
                if not have_end:
                    fo.write('##INFO=<ID=END,Number=1,Type=Integer,Description="DEL INV: POS+len(REF)-1. INS: POS">\n')
                if not have_svsub:
                    fo.write('##INFO=<ID=SVSUBTYPE,Number=1,Type=String,Description="DUP_TANDEM for large INS when inserted seq matches adjacent reference copy">\n')
                fo.write(line)
                break

        for line in fi:
            if not line or line[0] == "#":
                fo.write(line)
                continue

            F = line.rstrip("\n").split("\t")
            if len(F) < 8:
                fo.write(line)
                continue

            chrom, pos_s, vid, ref, alt, qual, flt, info = F[:8]
            try:
                pos1 = int(pos_s)
            except ValueError:
                fo.write(line)
                continue

            alt1 = alt.split(",")[0]
            info_d = parse_info(info)

            allele_type = (info_d.get("TYPE") or "").lower() if "TYPE" in info_d else None
            inv_flag = ("INV" in info_d)

            svtype, svlen, end1, svsub = infer_sv(
                chrom, pos1, ref, alt1, ref_fa,
                args.dup_window, allele_type, inv_flag,
                args.sv_min, args.inv_min
            )

            info_d["SVTYPE"] = svtype
            info_d["SVLEN"] = str(int(svlen))
            info_d["END"] = str(int(end1))
            if svsub:
                info_d["SVSUBTYPE"] = svsub

            F[7] = info_to_str(info_d)
            fo.write("\t".join(F) + "\n")

if __name__ == "__main__":
    main()

```

## 4. Variant annotation quality control and sanity checks

```bash
# To check the sanity

bcftools view -i 'INFO/SVTYPE="INV"' -H pb_pangenome.final.sensitive.vcf.gz | wc -l

# AF should be ≤1 per biallelic site; AN should reflect your non-missing calls
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/SVTYPE\n' pb_pangenome.final.vcf.gz | head

# counts by type
for t in INS DEL INV OTHER; do
  echo -n "$t: "; bcftools view -i "INFO/SVTYPE=\"$t\"" -H pb_45.final.vcf.gz | wc -l
done

```

## 5. VCF contig header correction and reindexing against reference genome

```bash

# rebuild contig header lines from reference
samtools faidx Pb3A_renamed.fa   # creates Pb3A_renamed.fa.fai if missing

# replace contig order in the VCF header with the FAI order
bcftools reheader -f Pb3A_renamed.fa.fai \
  -o pb_pangenome.reheaded.vcf.gz pb_pangenome.waved.biallelic.norm.fill.vcf.gz
bcftools index pb_45.reheaded.vcf.gz

# sort & index (after reheader)
bcftools sort -T ./tmp_sort \
  -Oz -o pb_pangenome.reheaded.sorted.vcf.gz pb_pangenome.reheaded.vcf.gz
bcftools index pb_pangenome.reheaded.sorted.vcf.gz

```

## 6. Variant type annotation and VCF conversion to conventional format

```bash

python map_type_to_svtype.py pb_pangenome.waved.biallelic.norm.fill.vcf.gz pb_pangenome.with_svtype.vcf
bgzip -f pb_pangenome.with_svtype.vcf
bcftools index pb_pangenome.with_svtype.vcf.gz

# Sanity check
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE\t%INFO/INV\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\n' \
  pb_45.with_svtype.vcf.gz | head

```

## 7. Variant class enumeration across the pangenome VCF

```bash
bcftools view -v snps pb_pangenome_final.svtype.fixed.vcf.gz | wc -l

bcftools view -v indels pb_pangenome_final.svtype.fixed.vcf.gz | wc -l


for t in SNP INDEL INS DEL INV OTHER; do
  echo -n "$t: "
  bcftools view -i "INFO/SVTYPE=\"$t\"" -H pb_pangenome_final.svtype.fixed.vcf.gz | wc -l
done

```


