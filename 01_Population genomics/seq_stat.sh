#!/bin/bash

#script to convert fasta to CDS, clean CDS and compute GC3 and pn/ps
# define file & directory names (full path)


outputdirCDS=$1
outprefix=$2
outputdir4fold=$3


## generate a new list with processed alignments
ls $outputdirCDS/ | grep ".fst.clean.fst.clean.fst" > $outprefix.list_CDS.txt


## compute summary statistics (-tstv = transition transervision ratio here fixed to 2 but can be set to another value)
echo "COMPUTING SUMMARY STATISTICS"
cd $outputdirCDS
../../../01_scripts/00_utils/seq_stat_coding -seq ../$outprefix.list_CDS.txt -f fasta -tstv 2 -code univ -o ../$outprefix.CDS.sumstats > ../$outprefix.CDS.sumstats.info
cd ..

# keep info of genes without premature stop codons
grep "stop" $outprefix.CDS.sumstats.info | awk '{print $1}' > $outprefix.CDS.withprematurestopcodons


### compute GC3s based on 4-fold degenerate sites (also includes prot alignments)
cd $outputdirCDS
echo $outputdirCDS
rm ../$outprefix.list_CDS.4fold

while read line; do
    python2 ../../../01_scripts/00_utils/script_python_sequencecodons4folddegenerateonly.py $line $line
    mv $line.sites4foldonly ../$outputdir4fold/
    echo "$outputdir4fold/$line.sites4foldonly" >> ../$outprefix.list_CDS.4fold # generate a list
done < ../$outprefix.list_CDS.txt
cd ..

# clean alignments
echo "A"
../../01_scripts/10_cleanAlignment -seq $outprefix.list_CDS.4fold -f fasta -n 4
ls $outputdir4fold/ | grep ".sites4foldonly.clean.fst" > $outprefix.list_CDS.4fold
echo "B"
# compute GC3s

echo "C"
cd $outputdir4fold
../../../01_scripts/00_utils/seq_stat_coding -seq ../$outprefix.list_CDS.4fold -f fasta -tstv 2 -code univ -o ../$outprefix.4fold.CDS.sumstats
cd ..
echo "D"

awk '{print $1" "$2"    "$3"    "$11}' $outprefix.4fold.CDS.sumstats | sed 's/GC3/GC3/g' > $outprefix.4fold.CDS.sumstats.GC3s

### merge pnps datasets & GC3s [require to exclude GC3s for seq with premature stop codons]
echo "name  Size    N   S   P   W   Ps  Pn  NbSS    D_Taj   GC3 name2    Size4fold  N4fold  GC3s    names3  PN  PS" > $outprefix.4fold.CDS.sumstats.final
less $outprefix.4fold.CDS.sumstats.GC3s | grep -v "Size" > tmp
less $outprefix.CDS.sumstats_PNPS | grep -v "Name"> tmp2 
while read line; do
    grep -v "$line" tmp > tmp3
    grep -v "$line" tmp2 > tmp4
    mv tmp3 tmp
    mv tmp4 tmp2

done < $outprefix.CDS.withprematurestopcodons
grep -v "Size" $outprefix.CDS.sumstats > tmp5
paste tmp5 tmp tmp2 >> $outprefix.4fold.CDS.sumstats.final # please check that this output is correct
rm tmp*


