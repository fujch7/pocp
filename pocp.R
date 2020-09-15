#!/usr/bin/env Rscript
#A R script to calculate pairwise POCP value using Blast
#Created by Jiacheng Fu, email:fujch@foxmail.com 
#If you find this useful, please cite: https://github.com/fujch7/pocp

library(data.table)
library(seqinr)
# set cutoff
idCutoff <- 40 # the value could be reset
queryCoverage <- 50 # 50 or 70


# =========================================================================
# =======================    Compute POCP    ==============================
# =========================================================================

if (!dir.exists('database')) {
  dir.create('database')
} else {
  unlink("database", recursive = TRUE)
}

# makeblastdb
genome.files <- list.files('input')
for (gn in genome.files) {
  header.file <- strsplit(gn,'.faa',fixed = T)[[1]][1]
  commond.makedb <- paste0('makeblastdb -dbtype prot -parse_seqids -in input/',
                           gn, ' -out database/', header.file)
  system(commond.makedb)
}

# pairwise blast
genome.comn <- combn(genome.files,2)
blast.comm1 <- 'blastp -num_threads 2 -query input/'
blast.comm2 <- paste0(' -db database/')
blast.comm3 <- paste0(' -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore slen" -out result/')
pocp.vector <- c()
for (i in (1:dim(genome.comn)[2]) ) {
  a.genome <- genome.comn[,i][1]
  b.genome <- genome.comn[,i][2]
  a.header <- strsplit(a.genome,'.faa',fixed = T)[[1]][1]
  b.header <- strsplit(b.genome,'.faa',fixed = T)[[1]][1]
  a.genome.seq <- read.fasta(paste0('input/', a.genome),'AA')
  b.genome.seq <- read.fasta(paste0('input/', b.genome),'AA')
  a.total <- length(a.genome.seq)
  b.total <- length(b.genome.seq)
  
  print(paste0('-- Blasting: ',a.header,' - VS - ',b.header))
  # blast forward
  result.forward <- paste0(a.header,'_VS_',b.header,'.tab')
  if (!file.exists(paste0('result/',result.forward))) {
    system(paste0(blast.comm1, a.genome,
                  blast.comm2, b.header,
                  blast.comm3, result.forward))

  } else {
    print(paste0(a.header,'_VS_',b.header,'.tab',' already exists!'))
  }
  
  # blast backward
  result.backward <- paste0(b.header,'_VS_',a.header,'.tab')
  if (!file.exists(paste0('result/',result.backward))) {
    system(paste0(blast.comm1, b.genome,
                  blast.comm2, a.header,
                  blast.comm3, result.backward))
  } else {
    print(paste0(b.header,'_VS_',a.header,'.tab',' already exists!'))
  }
  # computePOCP
  df.forward <- read.table(paste0('result/',result.forward),
                           header = F,sep = '\t',
                           stringsAsFactors = F)
  df.forward <- df.forward[which(df.forward$V3 > idCutoff & df.forward$V12 < 1e-5),]
  df.forward <- df.forward[which(df.forward$V4 > round((df.forward$V9*queryCoverage)/100,2)),]
  dt <- as.data.table(df.forward)
  setkey(dt, V1)
  df.forward <- dt[J(unique(V1)), mult = 'first']
  C1 <- dim(df.forward)[1]
  df.backward <- read.table(paste0('result/',result.backward),
                            header = F,sep = '\t',
                            stringsAsFactors = F)
  df.backward <- df.backward[which(df.backward$V3 > idCutoff & df.backward$V12 < 1e-5),]
  df.backward <- df.backward[which(df.backward$V4 > round((df.backward$V9*queryCoverage)/100,2)),]
  dt <- as.data.table(df.backward)
  setkey(dt, V1)
  df.backward <- dt[J(unique(V1)), mult = 'first']
  C2 <- dim(df.backward)[1]
  pocp <- round(((C1 + C2)/(a.total + b.total))*100, 2)
  pocp.vector <- append(pocp.vector, paste0(a.header,'  vs  ',b.header,'\t',pocp,'%'))
  print(paste0('-- Pair blast done: ',a.header,' - VS - ',b.header))
  print(paste0('-- The POCP : ', pocp, '%'))
  print(paste0('-- C1: ',C1))
  print(paste0('-- C2: ',C2))
  print(paste0('-- T1: ',a.total))
  print(paste0('-- T2: ',b.total))
}
write(pocp.vector, 'resultPOCP_blastp.txt')
