# Parsre user input
results_file <- commandArgs(trailingOnly = T)[1]
prefix<- commandArgs(trailingOnly = T)[2]

# read results file
results <- read.table(results_file, h=T)

# Calculate metrics
clean_data <- results[,c("name","size","n", "s")]
clean_data$pns <- results$polNS / (results$size - results$nss)
clean_data$ps <- results$polS / results$nss

#output

write.table(file = paste0(prefix, ".final.cleaned"), 
	    x = clean_data,
	    sep = "\t",
	    col.names = T,
	    row.names = F,
	    quote = F)
