## script for generating synthetic data based on rates
## it also includes some simple benchmarking with rule-based machine learning and decision trees
library(readxl)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)
library(IRanges)

path <- "//Users//matga374//Desktop//simulating_adar_edits//GSE169710_Data_S3_A-to-I_sites_editing_rate_identified_by_RNA-seq.xlsx"

# human brain
hs_edits <- read_excel(path, sheet = 6)
hs_edits <- as.data.frame(hs_edits)
# mice brain
mm_edits <- read_excel(path, sheet = 7)

seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, hs_edits[,1], start=hs_edits[,2]-100+1, end=hs_edits[,2]+100) # get the sequence at a given position
names(seq)
seq <- as.character(seq)
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("ensembl_gene_id","start_position","end_position","strand","hgnc_symbol","chromosome_name")
filters <- c("chromosome_name","start","end")

# filtration
hs_edits_er <- as.matrix(sapply(hs_edits[,3:4], as.numeric))  
hs_edits$seqs <- seq
hs_edits2 <- hs_edits[-which(is.na(rowMeans(hs_edits_er, na.rm = T))),]
rts <- as.data.frame(sapply(hs_edits2[,3:4], as.numeric))
nind <- c(which(rowMeans(rts) == 1),which(rowMeans(rts) == 0)) # we want only those that are either 1 rate or 0 (edit or no edit)
#hs_edits2 <- hs_edits2[which(hs_edits2[,1] %in% c("chr1","chr17")),]
hs_edits2 <- hs_edits2[nind,]

## loop over sites
values <- list(chromosome=gsub("chr","",hs_edits2[,1]),start=as.character(hs_edits2[,2]),end=as.character(hs_edits2[,2]))
all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes2 <- all.genes[-which(all.genes[,5] == ""),]

pb <- txtProgressBar(min = 1, max = length(nind), style = 3)
lst <- list()
j = 0
n = dim(hs_edits2)[1]
for(i in 1:n){
  
  coords1 <- which(all.genes2[,2] <= hs_edits2[i,2] & all.genes2[,3] >= hs_edits2[i,2])
  tmp_ag <- all.genes2[coords1, ]
  coords2 <- which(tmp_ag[,6] == gsub("chr","",hs_edits2[i,1]))
  all.genes3 <- tmp_ag[coords2,]
  
  if (length(coords2) >= 1){
    j = j + 1
    hgnc_symbol <- all.genes3[1,"hgnc_symbol"] # always first row as the are duplicates 
    er <- mean(suppressWarnings(as.numeric(hs_edits2[i,3:4])), na.rm = TRUE)
    sqnc <- unname(as.character(hs_edits2$seqs[i]))
    #sqnc <- reverse(chartr("ATGC","TACG", sqnc)) # reverse complementary, update: we don't use it now
    pos <- hs_edits2[i,2]
    chr <- hs_edits2[i,1]
    lst[[j]] <- c(chr, pos, hgnc_symbol, sqnc, er)
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

dfout <- as.data.frame(do.call(rbind, lst))
colnames(dfout) <- c("chr","pos","gene","seqs","labels")

# remove singles
gnsn <- which(dfout$gene %in% names(table(unlist(dfout$gene))[table(unlist(dfout$gene)) == 1]))
dfout <- dfout[-gnsn,]
dfout2 <- dfout[order(dfout$gene),]
seqs2 <- tolower(dfout2$seqs)

substr(seqs2, 100, 100) <- toupper(substr(seqs2, 100, 100))
dfout2$seqs <- seqs2
#dfout3 <- dfout2[sample(1:nrow(dfout2)), ]

dfout2 <- dfout2[which(substr(seqs2, 100, 100)=="A"),] # only these with A in the edit site

library(stringi)
library(tidyverse)

# generating some set of unique k-mers, in the future should be done based on drawing it from sequences itself
comb1 <- c("A","T","C","G")
comb2 <- combn(c("A","A","T","T","C","C","G","G"), 2)
comb23 <- combn(c("A","A","T","T","C","C","G","G"), 3)
comb2 <- c(apply(comb2, 2, paste, collapse = ""), apply(comb23, 2, paste, collapse = ""))
comb3 <- combn(c("A","A","A","T","T","T","C","C","C","G","G","G"), 3)
comb32 <- combn(c("A","A","A","T","T","T","C","C","C","G","G","G"), 2)
comb3 <- c(apply(comb3, 2, paste, collapse = ""), apply(comb32, 2, paste, collapse = ""))
comb4 <- combn(c(rep("A",4),rep("T",4),rep("C",4),rep("G",4)), 4)
comb4 <- apply(comb4, 2, paste, collapse = "")
comb5 <- combn(c("TA","TA","GA","GA","CA","CA","AA","TT","CC","GG","GC","TC","AC","GT","AT","CT"), 2) # repetitive ones
comb5 <- apply(comb5, 2, paste, collapse = "")
content <- c(comb1, comb2, comb3, comb4, comb5, stri_reverse(comb2), stri_reverse(comb3), stri_reverse(comb4), stri_reverse(comb5)) %>% unique

full_content <- list()
for(i in 1:length(dfout2$seqs)){
  sequence <- dfout2$seqs[i] %>% toupper
  # GC content added separately
  num_g <- str_count(sequence, "G")
  num_c <- str_count(sequence, "C")
  gc_content <- (num_g + num_c) / str_length(sequence)
  # other content
  x_content <- c()
  for(x in 1:length(content)){
    num_x <- str_count(sequence, content[x])
    x_content[x] <- num_x / str_length(sequence)
  }
  full_content[[i]] <- c(gc_content, x_content)
}

df_content <- do.call(rbind,full_content)
colnames(df_content) <- c("GCcont",content)
decision <- dfout2$labels
decision[decision=="1"] <- "edit"
decision[decision=="0"] <- "noedit"
df_content <- as.data.frame(df_content)
df_content$decision <- decision

## separate part, run after trees
rowsn <- c(which(df_content$decision=="edit")[1:500], which(df_content$decision=="noedit")[1:500])
mat <- df_content[rowsn, names(kmers_imp[1:20])]

library(pheatmap)
edits_rows <- data.frame(case=df_content$decision[rowsn])
rownames(edits_rows) <- rownames(mat)
pheatmap(as.matrix(mat), scale = "column", annotation_row = edits_rows,show_rownames = F, cutree_rows = 5)

### checking the position of some kmers

full_posits_edit <- list()
full_posits_noedit <- list()

kmer <- "AG"
midpo <- c()
for(i in 1:length(dfout2$seqs)){
  sequence <- dfout2$seqs[i] %>% toupper
  midpo[i] <- substr(dfout2$seqs[i],101,101)
  num_x <- gregexpr(kmer, sequence)[[1]] %>% as.numeric
  if(dfout2$labels[i] == "1"){
    full_posits_edit[[i]] <- num_x
    full_posits_noedit[[i]] <- NULL
  }else{
    full_posits_edit[[i]] <- NULL
    full_posits_noedit[[i]] <- num_x
  }

}
#hist(unlist(full_posits_edit), col="#B06EE0", xlab="position [bp]", breaks=100, main=paste0("k-mer: ", kmer, "cond: edit"))
#abline(v = 100, col='red', lwd = 1) # marking the place of edit

#hist(unlist(full_posits_noedit), col="#B06EE0", xlab="position [bp]", breaks=100, main=paste0("k-mer: ", kmer, "cond: noedit"))
#abline(v = 100, col='red', lwd = 1) # marking the place of edit

#plot(density(unlist(full_posits_edit)),col='#B06EE0')
#lines(density(unlist(full_posits_noedit)),col='#E0B16E')

tabd <- data.frame(values=c(unlist(full_posits_edit), unlist(full_posits_noedit)), condition=c(rep("edit",length(unlist(full_posits_edit))), rep("noedit",length(unlist(full_posits_noedit)))))
ggplot(data=tabd, aes(x=values, fill=condition)) +
  geom_density(alpha=0.5) +
  theme_minimal()+
  xlab("position [bp]")+
  ggtitle(paste0("k-mer: ", kmer))+
  geom_vline(xintercept = 100)

# RULE BASED MODELS
#library(R.ROSETTA)
#n1 <- which(df_content$decision=="edit")
#n0 <- sample(which(df_content$decision=="noedit"), length(n1), replace=F)
#df_content2 <- df_content[c(n0,n1),]

#ruleModelUS <- rosetta(df_content, underSample=TRUE, underSampleNum=2, underSampleSize=1000,  discreteMethod = "EqualFrequency",
#                       discreteParam = 2)
#ruleModelUS$quality

#ruleModelUS <- rosetta(df_content2)
#ruleModelUS$quality

#ruleModelUS

# DECISION TREES

library(rpart)
library(rpart.plot)

# prepare balanced data
# repeat n times
nt <- 1000
accuracy_train <- c()
bacc <- c()
vars_imp <- list()
pb <- txtProgressBar(min = 1, max = nt, style = 3)
for(i in 1:nt){
  n1 <- which(df_content$decision=="edit")
  n0 <- sample(which(df_content$decision=="noedit"), length(n1), replace=F)
  df_content2 <- df_content[c(n0,n1),]
  
  tree <- rpart(decision~., data = df_content2, method = "class", cp=0.001)
  
  confMat <- table(df_content2$decision, predict(tree, type = "class"))
  accuracy_train[i] <- sum(diag(confMat))/sum(confMat)
  
  conf.matrix <- round(prop.table(table(df_content2$decision, predict(tree, type="class")), 2), 2)
  #printcp(tree)
  #tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  
  bestcp <- tree$cptable[which.min(tree$cptable[,"xerror"]),"CP"]
  tree.pruned <- prune(tree, cp = bestcp)
  
  confMat <- table(df_content2$decision, predict(tree.pruned, type = "class"))
  
  cfs <- confusionMatrix(as.factor(df_content2$decision), predict(tree.pruned, type = "class"))
  bacc[i] <- unname(cfs$byClass["Balanced Accuracy"])
  
  vars_imp[[i]] <- tree.pruned$variable.importance
  setTxtProgressBar(pb, i)
}
close(pb)

kmers <- lapply(vars_imp, names) %>% unlist %>% unique
kmers_imp <- c()
for(i in 1:length(kmers)){
  kmers_imp[i] <- lapply(vars_imp, function(x) x[kmers[i]]) %>% unlist %>% mean(na.rm=T)
  names(kmers_imp[i]) <- kmers[i]
}
names(kmers_imp) <- kmers
kmers_imp <- kmers_imp %>% sort(decreasing = T)
barplot(kmers_imp[1:20], xlab="k-mer",ylab="avg. importance", col="#4ABDDE")

data <- data.frame(accuracy_train, bacc)
boxplot(data, names=c("training trees","balanced for pruned trees"), ylab="accuracy", col=c("#EEEE1F","#EE831F"))

# plots
plot(performance(fit.pred,"acc"))

rpart.plot(tree.pruned)

barplot(tree.pruned$variable.importance[1:5])


rpart.plot(tree.pruned, extra=104, box.palette="GnBu",
           branch.lty=3, shadow.col="gray", nn=TRUE)

library(ROCR)
fit.pr = predict(tree.pruned, type="prob")[,2]
fit.pred = prediction(fit.pr, df_content2$decision)
fit.perf = performance(fit.pred,"tpr","fpr")
plot(fit.perf,lwd=2,col="blue",
     main="ROC:  Classification Trees on Titanic Dataset")
abline(a=0,b=1)

###### GBM ##################################
# this code below was used only for testing
# it got similar results as trees above, but it was tricky to run classification
# as most of algorithms work for regression (continouous decision)
library(gbm)
library(caret) 

df_content$decision <- abs(as.numeric(as.factor(df_content$decision))-2)

n1 <- which(df_content$decision==1)
n0 <- sample(which(df_content$decision==0), length(n1), replace=F)
df_content2 <- df_content[c(n0,n1),]


parts = createDataPartition(df_content2$decision, p = 2/3, list = F)
train = df_content2[parts, ]
test = df_content2[-parts, ]

model_gbm = gbm(decision~.,
                data = train,
                distribution = "poisson",
                cv.folds = 10,
                n.trees = 10000)       # 500 tress to be built

summary(model_gbm)

#use model to make predictions on test data
pred_test = predict.gbm(object = model_gbm,
                        newdata = test,
                        n.trees = 500)

pred_test

# Give class names to the highest prediction value.
class_names = round(pred_test) #colnames(pred_test)[apply(pred_test, 1, which.max)]
result = data.frame(test$decision, class_names)

print(result)

conf_mat = confusionMatrix(as.factor(as.numeric(test$decision)), as.factor(as.numeric(class_names)))
print(conf_mat)

