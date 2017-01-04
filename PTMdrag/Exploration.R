
# Exploration of the PTM dataset

# Starting with the oxidation
setwd('~/Projects/ProteoDNN/PTMdrag')

dataRaw = read.table('data/OxidationAll.txt', header=T, sep=',', stringsAsFactors = F)
dataRaw = dataRaw[,1:20]

# removing the early and late retention times (large machine error)
dataPep = dataRaw[dataRaw$X.L.Retention.Time..min. >= 20 & dataRaw$X.L.Retention.Time..min. <= 80,]

pepLen = function(sqn){
  length(strsplit(sqn,'')[[1]])
}

# Add peptide length
dataPep$pepLen = sapply(dataPep$Peptide, pepLen)
dataPep$RetDcuts = cut(dataPep$pepLen, breaks = c(5,11,15,20,60))

ggplot(data=dataPep, aes(x=Retention..Delta.))+
  geom_histogram(bins=100, fill='steelblue')+
  facet_grid(.~RetDcuts)+
  theme_bw()
  
ggplot(data=dataPep, aes(x=Retention..Delta., y=X.L.m.z))+
  geom_density_2d()+
  theme_bw()