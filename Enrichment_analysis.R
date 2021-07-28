title: "Enrichment analysis"
author: "Ximena Fernandez"
date: "2/13/2021"
---

#ANALYSIS OF ENHANCEMENT
## To find or define the functionalities that have the greatest biological relevance in a pool of genes of interest. To do this, we use the fold-change, p-value (significance) and the fold discovery rate value.
# 1. Connect to the BioMart database. Install and load the liberia.

```{r}
install.packages('RCurl',dependencies = T, force=T)
BiocManager::install("RCurl", force = TRUE);
BiocManager::install("RSQLite", force = TRUE);

library(RSQLite)
library(RCurl)
library(biomaRt)
biomaRt::listMarts()
```

#We have 4 databases, from which we are going to use ENSEMBL_MART_ENSEMBL, which contains the genome of the main species. To load it, the useMart function will be used. This function has 2 parameters: mart and dataset. The ensembl database will be loaded first.
```{r}
mart<- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL')

biomaRt::listDatasets(mart=mart)
```

# We will work with the genome of the chimpace scientific name pan troglodytes and in ensembl (ptroglodytes_gene_ensembl) which is defined in the dataset parameter of the useMart function.
```{r}
mart<- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'ptroglodytes_gene_ensembl')
mart
str(mart)
```

#The attributes are all the variables that describe each gene of the selected organism, in this case they are 2958. But there are only 391 filters to select this information. If we want to find a filter related to a specific filter.

```{r}
searchFilters(mart=mart, pattern='GO')
```

# To create a data group with information from these databases we use the getBM function, which has four parameters:
+ _Attributes_: here you define the output columns that the generated object will contain.
+ _Filters_: this will be the column in which my selected terms will be searched, they can be genes.
+ _Values_: this parameter contains the values ​​to be searched in the column that is defined as a filter.
+ _Mart_: It is the mart object from which data will be taken.

#As an example we are going to create a database of GO terms of two randomly selected genes 100 and 230 corresponding to the ensemble ID. First we look for how the ID column is defined in the biomart database, for which the search attribute function is used.


```{r}
biomaRt::searchAttributes(mart=mart, pattern = 'ensembl')
```

The column that defines the ensembl ID is (ensembl_gene_id), this way you can search for columns for other data.
```{r}
geoannot<-getBM(attributes = c('ensembl_gene_id', 'go','external_gene_name'),
                filters = 'ensembl_gene_id',
                values = c(100,230),
                mart = mart)
```


```{r}
BiocManager::install('org.Hs.eg.db')
BiocManager::install('AnnotationDbi')
BiocManager::install('org.Pt.eg.db',dependencies = T, force = T)
library('AnnotationDbi')
library('org.Hs.eg.db')
library('org.Pt.eg.db')
keys(org.Pt.eg.db)
keytypes(org.Pt.eg.db) #es como ver el listMarts de biomart
```


```{r}
inmunome<-read.csv("C:/Users/usuario/Downloads/human_gene_set.txt")
View(inmunome)

goinmunome<-select(org.Hs.eg.db,keys = as.character(inmunome$EntrezGeneID),columns = "GO", keytype = "ENTREZID")
View(goinmunome) 

path <-select(org.Hs.eg.db,keys = as.character(inmunome$EntrezGeneID),columns = "PATH", keytype = "ENTREZID")
View(path)

protein <-select(org.Hs.eg.db,keys = as.character(inmunome$EntrezGeneID),columns = "UNIPROT", keytype = "ENTREZID")
View(protein)
```

```{r}
keytypes(org.Hs.eg.db)
```

```{r}
goinmunome1<-select(org.Hs.eg.db,keys = as.character(inmunome$EntrezGeneID),columns = c("UNIPROT","ALIAS"), keytype = "ENTREZID")
View(goinmunome1)
```

                    
```{r}
target<-read.csv("C:/Users/usuario/Downloads/human_gene_selected_set.txt")
install.packages(GOstats)
BiocManager::install('GOstats')
BiocManager::install('matrixStats')
library(Category)
library(matrixStats)
library(topGO)
library(GOstats)

a<-1
gene2go<-c()
genevec<-c()

for (geneid in inmunome$EntrezGeneID ) {
  genevec[a]<-geneid
  gene2go[a]<-list(goinmunome[goinmunome$ENTREZID==geneid,"GO"])
  a<-a+1
}
names(gene2go)<-genevec
#head(names(gene2go))
#head(inmunome)
genelist<-rep(0,nrow(inmunome))
genelist
genelist[inmunome$EntrezGeneID%in%target$EntrezGeneID]<-1 #si esta el gen en los genes objetivo
names(genelist)<-genevec
genelist<-factor(genelist)
```
For the enrichment analysis I compare my expected frequency with the one obtained.

```{r}
library(topGO)

GOdata.MF<-new("topGOdata",ontology="MF",description="Esta es una comparacion de genes de inmunidad innata",allGenes=genelist,annot=annFUN.gene2GO,gene2GO=gene2go)
```



```{r}
resultFisher.MF.classic <- runTest(GOdata.MF, algorithm =
                                     "classic", statistic = "fisher")
allRes.MF.classic <- GenTable(GOdata.MF, classicFisher =
                                resultFisher.MF.classic,topNodes=20)
allRes.MF.classic 
```
```{r}
resultFisher.MF.classic

```

```{r}
#GO_hyper_param<-new("GOHyperGParams", geneIds=target$EntrezGeneID, universeGeneIds=inmunome$EntrezGeneID, annotation="org.Hs.eg.db", pvalueCutoff=0.01, testDirection="over")

GO_hyper_param<-new("GOHyperGParams",
                    geneIds=target$EntrezGeneID,
                    universeGeneIds=inmunome$EntrezGeneID,
                    annotation="org.Hs.eg.db",
                    pvalueCutoff=0.01,
                    testDirection="over",
                    ontology="BP")

hypertest<-hyperGTest(GO_hyper_param)
summary(hypertest)
```

```{r}

GEOkeg<- new("KEGGHyperGParams",
             geneIds=target$EntrezGeneID,
                     universeGeneIds=inmunome$EntrezGeneID,
                     annotation="org.Hs.eg.db",
                     pvalueCutoff=0.01,
                     testDirection="over")

test_over_1<-hyperGTest(GEOkeg)
test_over_1

BiocManager::install("hgu95av2.db")
library(hgu95av2.db)
affylib<-'hgu95av2.db'
data(geneList)
```

```{r}
GOdata<-new("topGOdata", allGenes=geneList, geneSel=function(p){return(p<0.001)}, nodeSize=10, annot=annFUN.db, affyLib="hgu95av2.db", ontology="BP")

resultFisher<-runTest(GOdata, algorithm="classic", statistic = "Fisher")
resultFisher_ks<-runTest(GOdata, algorithm = "elim", statistic = "KS")

final_result<-GenTable(GOdata, Fisher=resultFisher, KS=resultFisher_ks, topNodes=20)

#final_result<-GenTable(GOdata,Fisher=resultFisher,KS=resultFisher_ks,orderBy=KS,topNodes=20)


```
