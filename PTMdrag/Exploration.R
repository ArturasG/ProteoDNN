
# Exploration of the PTM dataset

# Starting with the oxidation
setwd('~/Projects/ProteoDNN/ProteoDNN/PTMdrag')

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

calcPPss = function(seqs){
  AACompNames = c('Tiny', 'Small', 'Aliphatic', 'Aromatic', 'NonPolar', 'Polar', 'Charged', 'Basic', 'Acidic')
  propNames = c('Length',AACompNames, 'Mol_weight', 'Charge', 'pI', 'AliphaticI', 'InstabilityI', 'BomanI', 'Hydrophobicity', 'Hidrophob_moment')
  
  out = matrix(0, length(seqs), 18)
  Np = floor(length(seqs)/100)
  for( i in 1:length(seqs)) {
    out[i,] = calcPepProps(seqs[i])
    if(i %% Np == 0) cat(sprintf('%d%%..', i/Np))
    }
  colnames(out) = propNames
  out
}

pt_s1 = proc.time()
pepProps = do.call('rbind', lapply(dataPep$Peptide[1:5000], calcPepProps))
colnames(pepProps) = propNames
pt_e1 = -pt_s1 + proc.time()
pt_e1

pt_s2 = proc.time()
pepProps = calcPPss(dataPep$Peptide)
pt_e2 = -pt_s2 + proc.time()
pt_e2


save(list=c('pepProps','dataPep'), file='~/Projects/ProteoDNN/pepProps.RData')

#trying a random forest
library(randomForest)

set.seed(2004)
N = nrow(pepProps)
intrain = sample(1:N, floor(0.8*N))

trX = pepProps[intrain,]
trY = dataPep$Retention..Delta.[intrain]
teX = pepProps[-intrain,]
teY = dataPep$Retention..Delta.[-intrain]

RFmod = randomForest(trX, trY, ntree = 50)
preds = predict(RFmod, newdata = unname(teX))

save(list=c('pepProps','RFmod','preds'), file = '~/Projects/ProteoDNN/ProteoDNN/PTMdrag/RFmodel.RData')
MRSE = mean((preds-teY)^2)
MSE = mean((teY-mean(teY))^2)

MRSE_ = mean((preds[teY >0 & teY <20]-teY_)^2)
MSE_ = mean((teY_-mean(teY_))^2)

Errors = cbind(teY - preds, abs(teY-mean(teY)),abs(preds-teY))
Errors = cbind(teY,preds,Errors)
colnames(Errors) = c('Errors','teY','preds','AbsE','AbsRE')
ErrorsDF = as.data.frame(Errors)
ErrorsDF$charge = as.factor(round(pepProps[-intrain,12]))

ggplot(data=ErrorsDF, aes(x=AbsE, y=AbsRE, col=charge))+
  geom_point()+
  facet_grid(.~charge)+
  theme_bw()

tempErr = data.frame(Errors = c(Errors[,3],Errors[,4]),
                     ErrInt = c(as.character(cut(Errors[,3],breaks=c(0,0.5,1,2,5,10,100))), as.character(cut(Errors[,4],breaks=c(0,0.5,1,2,5,10,100)))),
                     ErrType = c(rep('AbsE',nrow(Errors)),rep('AbsRE',nrow(Errors))))

ggplot(ErrorsDF, aes(x=teY-mean(teY)))+
  geom_histogram(fill='steelblue', col='black', bins = 100)+
  theme_bw()+
  xlab('Ret_delta errors')

ggplot(ErrorsDF, aes(x=teY-preds))+
  geom_histogram(fill='steelblue', col='black', bins = 100)+
  theme_bw()+
  xlab('Ret_delta errors')


ggplot(data = tempErr, aes(x=factor(ErrInt, ordered=T, levels = c('(0,0.5]','(0.5,1]','(1,2]','(2,5]','(5,10]','(10,100]')) ))+
  geom_bar(stat='count')+
  facet_grid(.~ErrType)+
  theme_bw()+
  xlab('Retention delta errors')+
  ylab('Frequency')

plot(pepProps[-intrain,1], Errors[,4])

# --------- pred the retention time as well

trYRT = dataPep$X.L.Retention.Time..min.[intrain]
teYRT = dataPep$X.L.Retention.Time..min.[-intrain]

RFmodRT = randomForest(trX, trYRT, ntree = 50)
predsRT = predict(RFmodRT, newdata = unname(teX))

MRSE_RT = mean((predsRT-teYRT)^2)
MSE_RT = mean((teYRT-mean(teYRT))^2)

leftHump = which((teYRT-predsRT) < (-10))
errHumps = rep(0,length(predsRT))
errHumps[leftHump] = 1

errD = teYRT - predsRT
errDDF = data.frame(errDelta = errD, pepLen = teX[,'Length'], pepLenGrp = cut(teX[,'Length'],breaks = c(0,7,15,25,35,45,100))) 
ggplot(data=errDDF, aes(x=errDelta)) + geom_histogram(bins = 100, fill='steelblue', col='black') + theme_bw() + xlab('Ret_time prediction error')
ggplot(errDDF, aes(x=errDelta)) + geom_histogram(bins = 100) + facet_grid(.~pepLenGrp)+ theme_bw() + xlab('Ret_time prediction error')

corsMain = apply(teX[-leftHump,],2,function(x) cor(x,teYRT[-leftHump]))
corsHump = apply(teX[leftHump,],2,function(x) cor(x,teYRT[leftHump]))
cors = data.frame(vars = rep(names(corsMain),each=1), cors = c(corsMain, corsHump), split = rep(c('main','hump'),each=18))

ggplot(cors, aes(x=vars, y=cors, fill=as.factor(split)))+
  geom_bar(stat='identity', position='dodge', col='black')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  xlab('Variables')+
  ylab('Correlations')+
  labs(fill='Group')+
  #scale_fill_discrete(label='Split')

# looking at the sequences in hump
seqH = dataPep[-intrain,1][leftHump]
seqH_ = dataPep[-intrain,1][-leftHump]

calcAAfreq = function(seq){
  AAs = c('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')
  store = rep(0,20)
  seq = strsplit(toupper(seq), '')[[1]]
  for(i in seq) {
    pos = which(AAs == i)  
    store[pos] = store[pos] + 1
  }
  #store/sum(store)
  store[store!=0] = 1
  store
}

freqs_ = do.call('rbind',lapply(dataPep[-intrain,1], calcAAfreq))
plot2GroupsBars(colMeans(freqs_[errHumps==0,]),colMeans(freqs_[errHumps==1,]), AAs, c('Main','Hump'))

teXScaled = scale(teX)
plot2GroupsBars(colMeans(teXScaled[errHumps==0,]),colMeans(teXScaled[errHumps==1,]), colnames(teX), c('Main','Hump'), flipXlabs = T) + ggtitle('Scaled Values')

# ------------------------------------------

plot2GroupsBars = function(x, y, vars, labels=c('grp1','grp2'), flipXlabs = F){
  dat = data.frame(obs = c(x,y), labs = c(vars, vars), grp = c(rep(labels[1], length(x)),rep(labels[2], length(y))))
  p = ggplot(dat, aes(y=obs, x=as.factor(labs), fill=as.factor(grp)))+
        geom_bar(stat='identity', position='dodge', col='black')+
        theme_bw()+
        xlab('')+
        ylab('Observations')+
        labs(fill='Group')
  if(flipXlabs) p = p + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  p
}

#-------------------- try diff dataset --------------------  
# function to do everything:

doTheRF_PTM = function(datapath){
  dataRaw = read.table(datapath, header=T, sep=',', stringsAsFactors = F)
  dataRaw = dataRaw[,1:20]
  
  # removing the early and late retention times (large machine error)
  dataPep = dataRaw[dataRaw$X.L.Retention.Time..min. >= 20 & dataRaw$X.L.Retention.Time..min. <= 80,]
  
  pepProps = calcPPss(dataPep$Peptide)
  
  set.seed(100)
  N = nrow(pepProps)
  intrain = sample(1:N, floor(0.8*N))
  
  trX = pepProps[intrain,]
  trY = dataPep$Retention..Delta.[intrain]
  trY_RT = dataPep$X.L.Retention.Time..min.[intrain]
  teX = pepProps[-intrain,]
  teY = dataPep$Retention..Delta.[-intrain]
  teY_RT = dataPep$X.L.Retention.Time..min.[-intrain]
  
  RFmod = randomForest(trX, trY, ntree = 50)
  RFmodRT = randomForest(trX, trY_RT, ntree = 50)
  preds = predict(RFmod, newdata = unname(teX))
  predsRT = predict(RFmodRT, newdata = unname(teX))
  
  return(list(mod = RFmod, modRT = RFmodRT, preds = preds, predsRT = predsRT, teX = teX, teY=teY, teYRT = teY_RT))
  
}

Formylation = doTheRF_PTM('~/DataAnalysisProjects/PTM_Andrew/Formylation.txt')

Carbamylation = doTheRF_PTM('~/DataAnalysisProjects/PTM_Andrew/Carbamylation.txt')

Carbamylation_A = doTheRF_PTM('~/DataAnalysisProjects/PTM_Andrew/CarbaMylation_A.txt')
Carbamylation_D = doTheRF_PTM('~/DataAnalysisProjects/PTM_Andrew/CarbaMylation_D.txt')





#-------------------- "cross-validate the bad ones out" -------------------
















# --------------------Lets look at the largest errors -------------------- 

largeErrors = which(Errors[,'AbsRE']>=5)
medErrors = which(Errors[,'AbsRE']>1 & Errors[,'AbsRE']<5)










# finding closest peptides in properties (looking at outliers)

# distance
euclDist = function(x,y) sqrt(sum((x - y)^2))

findClosestN = function(x, Y, n, retIdx = T){
  dists = apply(Y, 1, function(rw) euclDist(rw,x))
  closest = order(dists)[1:n]
  if(retIdx) return(closest)
  else return(Y[closest,])
}


