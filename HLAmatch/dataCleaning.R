
# Data cleaning for HLA project

seq = read.csv('~/Projects/ProteoDNN/HLAmatch/data/SeqData.csv', sep=',', header=T, stringsAsFactors = F)
data = read.csv('~/Projects/ProteoDNN/HLAmatch/data/HLAtypes.csv', sep=',', header=T, stringsAsFactors = F)

data_gen = data[,c(1,9:21)]
data_gen_ = data_gen
data_gen_[,2:13] = do.call('rbind', 
                          lapply(1:nrow(data_gen_), 
                                 function(i) sapply(data_gen[i,2:13], 
                                                    function(hap) paste(substr(hap,1,1),
                                                                        substr(hap,3,4),
                                                                        substr(hap,5,6), 
                                                                        sep='_'))))


# modify variant codes to match the codes in the data
seq = seq[,2:ncol(seq)]
# remove missing seqs 
remsSeq = which(apply(seq,2,function(x) sum(is.na(x)))>0)
seq = seq[,-remsSeq]
variants = colnames(seq)
variants_mod = sapply(variants, function(x) {y = strsplit(x,'[.]')[[1]]; paste(y, collapse='_')})
colnames(seq) = variants_mod

# Patiends with id 343 and 188 had to be removed because it contained unclear/missing variants B_40_xx/C_15_01
remsSampl = which(apply(data_gen_[,2:13], 1, function(x) !all(x %in% variants_mod)))
data_gen_ = data_gen_[-remsSampl,]

# Expand the data to include all permutations of variants between two haplotypes
expanded_data = do.call(rbind, 
                        lapply(1:nrow(data_gen_), 
                               function(i) cbind(rep(data_gen_[i,1],8),
                                                 gen_pos_permutations(data_gen_[i,2:7]),
                                                 matrix(rep(data_gen_[i,8:13],each=8),nrow=8),
                                                 rep(data_gen_[i,14]))))

expanded_data = do.call(rbind, 
                        lapply(1:nrow(expanded_data), 
                               function(i) cbind(rep(expanded_data[i,1],8),
                                                 matrix(rep(expanded_data[i,2:7],each=8),nrow=8),
                                                 gen_pos_permutations(expanded_data[i,8:13]),
                                                 rep(expanded_data[i,14]))))

patient_ids_outcomes = expanded_data[,c(1,14)]
colnames(patient_ids_outcomes) = c('ID', 'outcome')
expanded_data_seq = matrix(unlist(expanded_data[,2:13]),ncol=12)

patient_seqs = do.call(rbind, lapply(1:nrow(expanded_data_seq), 
                                     function(i) convert_variants_to_seq(expanded_data_seq[i,],seq)))

write.table(patient_ids_outcomes, file='Projects/ProteoDNN/HLAmatch/data/meta.csv', sep='\t', row.names=F, col.names = T)
write.table(patient_seqs, file = 'Projects/ProteoDNN/HLAmatch/data/patientSeqData.csv', sep='\t', row.names = F, col.names=F)
#==================== functions ==================== 

gen_pos_permutations = function(datRow){
  tabl = list(c(1,2,3,4,5,6),
             c(1,2,3,4,6,5),
             c(1,2,4,3,5,6),
             c(1,2,4,3,6,5),
             c(2,1,3,4,5,6),
             c(2,1,3,4,6,5),
             c(2,1,4,3,5,6),
             c(2,1,4,3,6,5))
  do.call(rbind, lapply(1:8, function(i) datRow[tabl[[i]]]))
}

convert_variants_to_seq = function(vars, seq){
  do.call(c, lapply(vars, function(x) seq[,x]))
}