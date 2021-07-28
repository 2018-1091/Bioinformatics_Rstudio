## Initially, it was decided to analyze the ratio by nucleotides of the SARS-CoV-2 sequence. For this, the proteins S, ORF1, N, ORF7a, ORF3a, M and ORF6 were chosen.

nombres<-c("S_protein","ORF1ab","N_protein","ORF7a","ORF3a","M_protein","ORF6")
genome_sequence<-c("NC_045512.2")
genes<-ape::read.GenBank(genome_sequence, as.character = T)

##Indicating the region of each protein
S_protein<-toupper(genes$NC_045512.2[21563:25384])
ORF1ab<-toupper(genes$NC_045512.2[266:21555])
N_protein<-toupper(genes$NC_045512.2[28274:29533])
ORF7a<-toupper(genes$NC_045512.2[27394:27759])
ORF3a<-toupper(genes$NC_045512.2[25393:26220])
M_protein<-toupper(genes$NC_045512.2[26523:27191])
ORF6<-toupper(genes$NC_045512.2[27202:27387])

##Function to find nucleotide ratio
proteins<-list(S_protein,ORF1ab,N_protein,ORF7a,ORF3a,M_protein,ORF6)
names(proteins)<-nombres

##Finding the proportion
proporcion<-data.table()
longitud<-c(length(proteins$S_protein),length(proteins$ORF1ab),length(proteins$N_protein),length(proteins$ORF7a),length(proteins$ORF3a),length(proteins$M_protein),length(proteins$ORF6))

for (i in 1:7){
  a<-table(proteins[i])/longitud[i]
  proporcion<-rbind(proporcion,a)
}
proporcion$proteina<-rep(nombres,each = 4)

##Sorting table
tabla<-matrix(proporcion$N, ncol = 7,nrow = 4)
rownames(tabla)<-c("A","C","G","T")
colnames(tabla)<-nombres

## Descriptive analysis - verifying that the proportions of GC and AT are the same in the same protein
barplot(tabla, beside = TRUE,col = 2:5)

China<-list(China_Mprotein,China_Nprotein,China_Sprotein)##Streams of China (Wuhan, Honkong, Beijing, Guangzhou [Canton], zhejiang)
USA<-list(USA_Mprotein,USA_Nprotein,USA_Sprotein) ##US sequences (Florida, California, Texas, Connecticut, Massachusett)

China_S_promedio<-data.table()
USA_S_promedio<-data.table()
for (i in 1:5){
  a<-table(China$Sprotein[i])/1273
  b<-table(USA$Sprotein[i])/1273
  China_S_promedio<-rbind(China_S_promedio,a)
  USA_S_promedio<-rbind(USA_S_promedio,b)
}
datos<-rbind(China_S_promedio, USA_S_promedio)
datos$pais<-rep(c("China","USA"),each=100)
datos$proteina<-rep(c("Sprotein"),nrow(datos))
datos

boxplot(datos$N~datos$pais+datos$V1,las=3, cex=0.8)

datos1<-subset(datos, subset = pais=="China" & V1=="Y")
datos2<-subset(datos, subset = pais=="USA" & V1=="Y")
t.test(datos1$N,datos2$N)

#Anova
fit<-aov(datos$N~datos$V1+datos$pais)
tukeyobj<-TukeyHSD(fit)
summary(fit)

tukeyobj

Finally, it was decided to analyze the N and S protein in different types of viruses belonging to the Coronaviridae family.

codigos<-c("NC_045512.2", "NC_034440.1","NC_006577.2") #Severe acute respiratory syndrome coronavirus 2 / Bat coronavirus / Human coronavirus HKU1

sequence<-ape::read.GenBank(codigos, as.character = T)
S_protein_COVID19<-toupper(sequence$NC_045512.2[21563:25384])
S_protein_BAT<-toupper(sequence$NC_034440.1[21232:25269])
S_protein_HKU1<-toupper(sequence$NC_006577.2[22942:27012])
N_protein_COVID19<-toupper(sequence$NC_045512.2[28274:29533])
N_protein_BAT<-toupper(sequence$NC_034440.1[28311:29558])
N_protein_HKU1<-toupper(sequence$NC_006577.2[28320:29645])
lista<-list(S_protein_COVID19,S_protein_BAT,S_protein_HKU1, N_protein_COVID19, N_protein_BAT, N_protein_HKU1)

titulos<-c("S_protein_SARS_COV2","S_protein_BAT","S_protein_HKU1", "N_protein_SARS_COV2", "N_protein_BAT", "N_protein_HKU1")
value=data.table()
for (gen in 1:length(lista)){
  a=round(table(lista[[gen]])/length(lista[[gen]]),digits = 3)
  value<-c(value,as.vector(a))
  print(a)
  barplot(a, col = c("royalblue", "seagreen", "purple", "grey"),las=1, main=titulos[gen])
  }
  
  datos3<-data.table()
for (i in 1:length(value)){
  datos3<-rbind(datos3,value[[i]] )
}

datos3$proteina<-rep(c("S_protein","N_protein"), each=12)
datos3$virus<-rep(c("SARS_CoV2","BAT coronavirus","HKU1 coronavirus","SARS_CoV2","BAT coronavirus","HKU1 coronavirus"),each=4)
datos3$base<-rep(c("A","C","G","T"),6)
datos3

boxplot(datos3$x~datos3$virus+datos3$base, cex=0.8)

fit1<-aov(datos3$x~datos3$virus+datos3$base)
tukeyobj1<-TukeyHSD(fit1)
summary(fit1)

tukeyobj1


#Packages:
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
#install.packages("forcats")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readr")
#install.packages("tidyr")
#install.packages("survminer")
#BiocManager::install("GEOquery")
#BiocManager::install("limma")
#BiocManager::install("pheatmap")
#BiocManager::install("hgu133plus2.db")
#BiocManager::install("pd.hg.u133.plus.2")

library(limma)
library(GEOquery)
library(gplots)
library(annotate)
library(topGO)
library(hgu133plus2.db)
library(ggplot2)
library(pd.hg.u133.plus.2)


data <- getGEO("data164805", dataMatrix =TRUE, AnnotGPL=FALSE)

```{r}
data <- data[[1]]
exprs(data) <- log2(exprs(data))
boxplot(exprs(data),outline=FALSE)

```
```{r}
sampleInfo <- pData(data)
sampleInfo
```

```{r}
design <- model.matrix(~0+sampleInfo$source_name_ch1)
design
```
```{r}
colnames(design) <- c("Normal","COVID_no_severo","COVID_severo")
```

```{r}
summary(exprs(data))

## calculate median expression level
cutoff <- median(exprs(data))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(data) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
data <- data[keep,]
```

```{r}
fit <- lmFit(exprs(data), design)
head(fit$coefficients)
```

```{r}
contrasts <- makeContrasts(COVID_severo - COVID_no_severo, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
```

```{r}
fit2 <- eBayes(fit2)
topTable(fit2)
```

```{r}
anno <- fData(data)
anno
```

```{r}
anno_genes<-cbind(anno$ID, anno$ORF)
anno_genes<-data.frame(anno_genes)

#tab<-topTable(fit2,number = 3e6)
tab$gene<-""
for (i in 1:nrow(tab)){
  for (j in 1:nrow(anno_genes)){
    if (row.names(tab)[i]==anno_genes$X1[j]){
      tab$gene[i]<-anno_genes$X2[j]
    }
  }
}

}```


