

library("mmgenome")
options(scipen = 8)



if (file.exists("data.RData")){
  load(file="data.RData")
  stop("Data already generated, loading it directly!")  
}



S1 <- read.csv("data/S1.csv", header = T)               
S2 <- read.csv("data/S2.csv", header = T)

gc <- read.delim("data/assembly.gc.tab", header = T)
kmer <- read.delim("data/assembly.kmer.tab", header = T)
ess <- read.table("data/assembly.orfs.hmm.id.txt", header = F)
cons.tax <- read.delim("data/tax.txt", header = T)

a16S<-read.delim("data/rRNA/16S.csv", header = T, sep = ";")
a23S<-read.delim("data/rRNA/23S.csv", header = T, sep = ";")

assembly <- readDNAStringSet("data/assembly.fa", format = "fasta")

network <- read.delim("data/network.txt", header = T)

colnames(kmer)[1] = "scaffold"
colnames(ess) = c("scaffold","orf","hmm.id")



d <- cbind.data.frame(gc$contig,                       
                      gc$gc, 
                      S1$Reference.length,
                      S1$Average.coverage, 
                      S2$Average.coverage)
colnames(d) = c("scaffold", "gc", "length", "S1", "S2")



rda <- rda(kmer[,2:ncol(kmer)])
d<-cbind(d,scores(rda,choices=1:3)$sites)



tax <- clean_tax(cons.tax, occurence=20, expand="Proteobacteria", clean = T)
d <- merge(d,tax[,c(1,2)], by = "scaffold", all = T)



t16S <- import_rrna(data=a16S, type="16S")
t23S <- import_rrna(data=a23S, type="23S")
d <- merge(x = d, y = t16S , by = "scaffold", all = T)
d <- merge(x = d, y = t23S , by = "scaffold", all = T)



g <- import_network(network=network, data=d, nconnections = 10)
d <- merge_network(graph = g, data = d)



e <- merge(ess, d, by = "scaffold", all.x = T)



rgb.p<- colorRampPalette(c('red','green','blue'))
gc.trans<-adjustcolor(rgb.p(max(d$gc)-min(d$gc)),alpha.f=0.3)



rm(list = c("S1","S2","a16S","a23S","cons.tax","ess","gc","kmer","network","t16S","t23S","rda","rgb.p","tax"))
save.image(file="data.RData")


