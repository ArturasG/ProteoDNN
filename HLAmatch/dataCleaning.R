
# Data cleaning for HLA project

seq = read.csv('~/Projects/ProteoDNN/HLAmatch/data/SeqData.csv', sep=',', header=T)
data = read.csv('~/Projects/ProteoDNN/HLAmatch/data/HLAtypes.csv', sep=',', header=T)

haplotypes = colnames(seq)[2:ncol(seq)]

