#!/usr/bin/R

library(pvclust)

# first column is codon names, rest are data
df <- read.delim("data.tsv", row.names=1)

# takes about 10 minutes for 100000
# ~12 hours for 1 million
result <- pvclust(df,method.dist="cor",method.hclust="average",nboot=10000)

plot(result)

pvrect(result,alpha=0.95)

seplot(result) 
