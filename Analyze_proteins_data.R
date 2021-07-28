#1. Package import ====

#Packages to APP protein
install.packages("seqinr",dependencies = TRUE)
install.packages("ape",dependencies = TRUE)
library(seqinr)
library(ape)
path1<-"C:/Users/usuario/Downloads/APP_RNA.fasta"
app.Hs<-seqinr::read.fasta(path)
app.Hs
path2<-"C:/Users/usuario/Downloads/APP_PROT.fasta"

app.Hs<-seqinr::read.fasta(path2,seqtype = "AA")
app.Hs
typeof(app.Hs)
class(app.Hs)
class(app.Hs$NP_000475.1)
barplot(sort(table(app.Hs)),col=1:15)
barplot(sort(table(app.Hs)/length(app.Hs$NP_000475.1)),col=1:15)

#Many proteins
path3<-"C:/Users/usuario/Downloads/sequence.fasta"
app.Hs3<-seqinr::read.fasta(path3,seqtype = "AA") #Lista
names(app.Hs3)<-c("mus","rattus","danio","gallina")
str(app.Hs3)
p_mus<-sort(table(app.Hs3$mus))/length(app.Hs3$mus)
p_danio<-sort(table(app.Hs3$danio))/length(app.Hs3$danio)
p_rattus<-sort(table(app.Hs3$rattus))/length(app.Hs3$rattus)
p_gallina<-sort(table(app.Hs3$gallina))/length(app.Hs3$gallina)


p_mus<-data.frame(a=names(p_mus),valor=as.vector(p_mus))
p_danio<-data.frame(a=names(p_danio),valor=as.vector(p_danio))
p_rattus<-data.frame(a=names(p_rattus),valor=as.vector(p_rattus))
p_gallina<-data.frame(a=names(p_gallina),valor=as.vector(p_gallina))
dim(p_mus)
dim(p_danio)
dim(p_rattus)
dim(p_gallina)
db<-merge(x=p_mus,y=p_danio, by="a")
db<-merge(x=db,y=p_rattus, by="a")
db<-merge(x=db,y=p_gallina, by="a")
db
rownames(db)<-db[,1]
db<-db[,2:ncol(db)]
colnames(db)<-c("mus","rattus","danio","gallina")
db<-t(db)
db
summary(db)
boxplot(y~x)

#Protein insulin
path4<-"C:/Users/usuario/Downloads/insulina.fasta"
insulina<-seqinr::read.fasta(path4,seqtype = "AA") #Lista
names(insulina)<-c("mus","rattus","danio","gallina")
insulina

p_insulina_mus<-sort(table(insulina$mus))/length(insulina$mus)
p_insulina_danio<-sort(table(insulina$danio))/length(insulina$danio)
p_insulina_rattus<-sort(table(insulina$rattus))/length(insulina$rattus)
p_insulina_gallina<-sort(table(insulina$gallina))/length(insulina$gallina)

p_insulina_mus<-data.frame(a=names(p_insulina_mus),valor=as.vector(p_insulina_mus))
p_insulina_danio<-data.frame(a=names(p_insulina_danio),valor=as.vector(p_insulina_danio))
p_insulina_rattus<-data.frame(a=names(p_insulina_rattus),valor=as.vector(p_insulina_rattus))
p_insulina_gallina<-data.frame(a=names(p_insulina_gallina),valor=as.vector(p_insulina_gallina))

db1<-merge(x=p_insulina_mus,y=p_insulina_rattus, by="a")
db1<-merge(x=db1,y=p_insulina_danio, by="a")
db1<-merge(x=db1,y=p_insulina_gallina, by="a")

db2<- data.frame(org=rep(c('mus','rattus','danio','gallina'),each=20),
                 prot=rep('amyloide',80),
                 AA=rep(db[,1],4),
                 value=c(db[,2],db[,3],db[,4],db[,5])
                 
)
db2

db3<- data.frame(org=rep(c('mus','danio','rattus','gallina'),each=20),
                 prot=rep('insulina',80),
                 AA=rep(db1[,1],4),
                 value=c(db1[,2],db1[,3],db1[,4],db1[,5])
                 
)
db3
value_I<-c(raton$value,rata$value,danio$value,gallina$value)

df<-data.frame(organismo=c(organimo_APP,organismo_I)),
proteina=c(proteina_APP,proteina_I)
AA=c(AA_APP,AA_I)
value=c(value_APP,value_I)
head(df)
boxplot(db3$value~db3$prot+db3$AA,las=2)

df
#Analysis doble variable 
fit<-lm(df$value~df$proteina*df$AA)
fit<-aov(df$value~df$proteina*df$AA)

tukeyobj<-TukeyHSD(fit)
tukeyobj[3]
View(tukeyobj[[3]])
TukeyHSD(fit)
summary(fit)

seqinr::choosebank(infobank=T)
seqinr::choosebank("swissprot")
APP<-seqinr::query("APP","sp=Homo sapiens AND K=Amyloid")
print(APP)
length(APP)
class(APP)
summary(APP)
attributes(APP)
seqinr::getName(APP)
seqinr::getAnnot(APP$req[3])
seqinr::getSequence(APP$req[3],as.string=T)
table(seqinr::getSequence(APP$req[3]))


names(APP)
APP$call
APP$name
APP$nelem
APP$typelist
APP$socket
APP$req
class(APP$req)
length(APP$req)
names(APP$req)
APP$req[1]
APP$req[2]


##Using APE====
#APE
library(ape)
ape
help("read.GenBank")
ids<-c("NC_045512.2","NC_045512.2","NC_004718.3","NC_045512.2","NC_004718.3") #lista
mseq<-ape::read.GenBank(ids,as.character = T)
mseq
structure(mseq)
str(mseq)
names(mseq)
mseq$NC_045512.2
seqinr::rho(mseq$NC_045512.2) #ProporciÃ³n de dinucleÃ³tidos  que hay en la secuencia 
plot(seqinr::rho(mseq$NC_045512.2)) #Si la variaiones o mutaciones suceden al azar 
barplot(seqinr::rho(mseq$NC_045512.2)-1, main=names(mseq),space=0.1,width=0.5,las=2, col=1.7,ylab="frecuency",xlab="dinucleotido")
## UTR sirve para marcar un producto, para ver si hay sobrerepresiÃ³n de una expresiÃ³n en especÃ­fico
sort(seqinr::rho(mseq$NC_045512.2),decreasing=T)

seqinr::zscore(sequence=mseq$NC_045512.2,modele="codon",simulations=1000,exact=T)
seqinr::rho(mseq$NC_045512.2,wordsize=8)
sort(seqinr::rho(mseq$NC_045512.2,wordsize=8),decreasing=T)
head(sort(seqinr::rho(mseq$NC_045512.2,wordsize=8),decreasing=T))
barplot(seqinr::rho(mseq$NC_045512.2,wordsize=8))
seqinr::GC(mseq$NC_045512.2) #muy burdo para anlizar un genoma sencillo
for (win.size in c(20,50,100,500)){ #tamaÃ±o de la ventana
  gc.perc<-vector()
  k<-1
  for (i in seq(from=1,to=length(mseq[[1]])-win.size,by=win.size/10)){
    j<-i+win.size-1
    gc.perc[k]<- seqinr::GC(mseq[[1]][i:j])
    k<-k+1
  }
  gc.perc <- gc.perc[!is.na(gc.perc)]
  plot(gc.perc*100,type="l", ylim = c(0,100),
       main = win.size, ylab="GC%", xlab="", col=1:15)
}

## Analyze proteins

seqinr::AAstat(seqinr::getSequence(insulina$rattus))


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
library(GenomeGraphs)
mart<-useMart("ensembl",dataset = "hsapiens_gene_ensembl")
listDatasets(mart)
gene<-makeGene(id="ENSG00000142192",type = "ensembl_gene_id",biomart = mart)
gene <- makeGene(id = "ENSG00000142192",type = "ensembl_gene_id", biomart = mart)

view(biomaRt::listFilters(mart))
searchFilters(mart,"ensembl")

gdPlot(list(gene,transcript))
transcript<-makeTranscript(id = "ENSG00000142192", type="ensembl_gene_id", biomart=mart)
anotacion<-makeAnnotationTrack(start=c(25940000,26020000,26140000),end=c(25980000,26100000,26170000),feature=c("Regulation","Deletion","Regulation"),dp=DisplayPars(Regulation="gray",Deletion="black"))
gdPlot(list(makeTitle("APP"),makeIdeogram(chromosome=21),gene,transcript, makeGenomeAxis(),anotacion))

#3. Database edition ====
sumary(ca) #General Properties
dim(ca) #dimension
ca<-ca[,-c(1,2)] # col
head(ca) # 6 rows
str(ca)
ca$Chromosome<as.factor(ca$chromosome) 
table(ca$chromosome) #Absolute count of each chromosome
barplot(table(ca$Chromosome))
table(ca$Alt)
table(ca$Ref)
ca[is.na(ca$Gene),]
ca<-ca[1:44,]
dim(ca)
table(ca$Zygosity,ca$Ref)

