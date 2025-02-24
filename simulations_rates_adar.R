## script for generating synthetic data based on rates
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
  #sqnc <- reverse(chartr("ATGC","TACG", sqnc)) # reverse complementary, update: we don't use it now*
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

substr(seqs2, 101, 101) <- toupper(substr(seqs2, 101, 101))
dfout2$seqs <- seqs2
#dfout3 <- dfout2[sample(1:nrow(dfout2)), ]

# saving files in different settings
write.csv(dfout2, "//Users//matga374//Desktop//simulated_data//simulated_data_1or0_nt_all.csv")

dfout2A <- dfout2[which(substr(seqs2, 101, 101)=="A"),] # only ones with A
dfout2T <- dfout2[which(substr(seqs2, 101, 101)=="T"),] # only ones with T

write.csv(dfout2A, "//Users//matga374//Desktop//simulated_data//simulated_data_1or0_nt_A.csv")
write.csv(dfout2T, "//Users//matga374//Desktop//simulated_data//simulated_data_1or0_nt_T.csv")

# now balanced, cause in DL training model was biased when the data was not balanced
n1 <- which(dfout2$labels=="1")
n0 <- sample(which(dfout2$labels=="0"), length(n1), replace=F)

dfout3 <- dfout2[c(n0,n1),]
write.csv(dfout3, "//Users//matga374//Desktop//simulated_data//simulated_data_1or0_nt_balanced.csv") # this is the data set I mostly focus on

# *this code needs to be redone cause all the simulated data are generate include reverse complementary sequences
