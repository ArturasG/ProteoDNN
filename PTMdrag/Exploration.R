
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

ggplot(data=dataPep, aes(y=Retention..Delta., x=pepLen, col=as.factor(X.H.Charge)))+
  geom_density_2d()+
  facet_grid(.~X.H.Charge)+
  theme_bw()

ggplot(data=dataPep, aes(y=Retention..Delta., x=pepLen, z=pepLen))+
  stat_summary_hex(fun=function (x) log10(length(x)))+
  scale_fill_gradient(low = 'blue', high = 'red')+
  facet_grid(.~X.H.Charge)+
  theme_bw()


ggplot(data=dataPep, aes(y=Retention..Delta., x=X.H.Mass, z=pepLen))+
  stat_summary_hex(fun=length)+
  scale_fill_gradient(low = 'blue', high = 'red')+
  facet_grid(.~X.H.Charge)+
  theme_bw()

ggplot(data=dataPep, aes(y=Retention..Delta., x=X.H.Mass, z=Mass.Difference..Delta.))+
  stat_summary_hex(fun=mean)+
  scale_fill_gradient(low = 'blue', high = 'red')+
  facet_grid(.~X.H.Charge)+
  theme_bw()

ggplot(data=dataPep, aes(y=X.H.Mass, x=pepLen, z=Retention..Delta.))+
  stat_summary_hex(fun=mean)+
  scale_fill_gradient(low = 'blue', high = 'red')+
  facet_grid(.~X.H.Charge)+
  theme_bw()

# Calculating physicochemical properties of peptides for learning
require(Peptides)

calcAAComp = function(seq){
  aacomp(seq)[,2]*0.01
}

calcPepProps = function(seq){
  c(lengthpep(seq=seq), # peptide length
    calcAAComp(seq=seq), # amino acid composition
    mw(seq=seq), # molecular weight
    charge(seq),
    pI(seq=seq), # isoelectric point
    aindex(seq=seq), # aliphatic index
    instaindex(seq=seq), # instability index
    boman(seq=seq), # boman index
    hydrophobicity(seq=seq), # hydrophobicity index
    hmoment(seq=seq) # hydrophobicity moment
  )
}

AACompNames = c('Tiny', 'Small', 'Aliphatic', 'Aromatic', 'NonPolar', 'Polar', 'Charged', 'Basic', 'Acidic')
propNames = c('Length',AACompNames, 'Mol_weight', 'Charge', 'pI', 'AliphaticI', 'InstabilityI', 'BomanI', 'Hydrophobicity', 'Hidrophob_moment')

pepProps = do.call('rbind', lapply(dataPep$Peptide, calcPepProps))
colnames(pepProps) = propNames

#trying a random forest
library(randomForest)

set.seed(2002)
N = nrow(pepProps)
intrain = sample(1:N, floor(0.8*N))

trX = pepProps[intrain,]
trY = dataPep$Retention..Delta.[intrain]
teX = pepProps[-intrain,]
teY = dataPep$Retention..Delta.[-intrain]

RFmod = randomForest(trX, trY)
preds = predict(RFmod, newdata = unname(teX))

save(list=c('pepProps','RFmod','preds'), file = '~/Projects/ProteoDNN/ProteoDNN/PTMdrag/RFmodel.RData')
MSE = mean((preds-teY)^2)

# finding closest peptides in properties (looking at outliers)

# distance
euclDist = function(x,y) sqrt(sum((x - y)^2))

findClosestN = function(x, Y, n, retIdx = T){
  dists = apply(Y, 1, function(rw) euclDist(rw,x))
  closest = order(dists)[1:n]
  if(retIdx) return(closest)
  else return(Y[closest,])
}


pRFmod = parallelRandomForest::randomForest(trX, trY)
preds = predict(RFmod, newdata = unname(teX))


