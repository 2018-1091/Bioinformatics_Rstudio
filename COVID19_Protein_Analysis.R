#Colors
colores1 <- c(rgb( 121/255, 159/255,167/255,alpha = 1),
              rgb( 206/255, 232/255,221/255,alpha = 1),
              rgb( 250/255, 240/255,230/255,alpha = 1),
              rgb( 245/255, 213/255,202/255,alpha = 1),
              rgb( 219/255, 154/255,150/255,alpha = 1)
)

colores2 <- c(rgb( 121/255, 159/255,167/255,alpha = .5),
              rgb( 206/255, 232/255,221/255,alpha = .5),
              rgb( 250/255, 240/255,230/255,alpha = .5),
              rgb( 245/255, 213/255,202/255,alpha = .5),
              rgb( 219/255, 154/255,150/255,alpha = .5)
)

## Installing packages with dependencies

install.packages("seqinr", dependencies = T)
install.packages("ape", dependencies = T)
library(seqinr)
library(ape)

## 1. Loading a sequence manually ====

path <- "D:/Usuarios/Administrador/Desktop/Bioinfo_avanzada/Clase2/sequence.fasta"
actin.necator <- seqinr::read.fasta(path, seqtype = "AA")


attr(actin.necator$XP_013298257.1,"Annot")
tabla <- table(actin.necator$XP_013298257.1)
sort(table(actin.necator$XP_013298257.1))
sort(table(actin.necator$XP_013298257.1),decreasing = T)
plot(sort(table(actin.necator$XP_013298257.1),decreasing = T), col=colores1)
barplot(sort(table(actin.necator$XP_013298257.1),decreasing = T), col=colores1)
longitud <- length(actin.necator$XP_013298257.1)
freqre <- tabla/longitud
barplot(freqre)

#2. Loading multiple sequences ====

path <- "D:/Usuarios/Administrador/Desktop/Bioinfo_avanzada/Clase2/app_Hs_multiple_seq.txt"
app.Hs <- seqinr::read.fasta(path, as.string = T, seqtype = "AA")

#3. Importing sequences from R ====

seqinr::choosebank(infobank = T)
seqinr::choosebank("swissprot")
app2 <- query("app2", "K=Amyloid AND sp=Homo sapiens")
summary(app2)
app2$call
app2$nelem
app2$req
app2$req[[2]]
seqinr::getSequence(app2$req[[2]])

#4. Using ape ====

acs.num <- c("MT419843.1","MT450923.1","MT419850.1",
             "MT419838.1","MT419849.1","MT419840.1")
library(ape)
cov2 <- ape::read.GenBank(acs.num, as.character = T)
cov2

For this part a package called CovidMutaciones () was used. Which was downloaded from the GitHub repository (devtools :: install_github ("MSQ-123 / CovidMutations")).

```{r}
library(CovidMutations)
#The example data:
data("nucmer")
head(nucmer)
# Fix IUPAC codes
nucmer<-nucmer[!nucmer$qvar%in%c("B","D","H","K","M","N","R","S","V","W","Y"),]
nucmer<- mergeEvents(nucmer = nucmer)## This will update the nucmer object
```

+ Complete annotation of all mutations identified by this study. The columns are described here. Sample: GISAID sample ID; refpos: position in the reference genome NC_045512.2; refvar: nucleotide composition of the reference at the refpos coordinate (a "." indicates an insertion); qvar: variant in the query sample (a "." indicates a deletion); qlength: length of the query genome (the reference genome is always 29,903 nucleotides long)

To provide effects of each SNP, insertion, and deletion in the virus genome.
```{r}
data("refseq")
data("gff3")

head(gff3)
annot <- setNames(gff3[,10],gff3[,9])  #annot: subset the gene and its product.
outdir <- tempdir()
head(indelSNP(nucmer = nucmer,
         saveRda = FALSE,
         refseq = refseq,
         gff3 = gff3,
         annot = annot,
         outdir = outdir))
```

The package also graphs the mutation statistics
```{r}
data("covid_annot")
head(covid_annot)
covid_annot <- as.data.frame(covid_annot)
```
