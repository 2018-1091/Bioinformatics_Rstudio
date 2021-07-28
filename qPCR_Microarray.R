BiocManager::install("ggplot2",dependencies=T)
BiocManager::install("rmdl",dependencies=T)

library(ggplot2)
ggplot(data=qPCR_intensities, aes(x=Cycles, y=Fluorescence))+geom_line()+ggtitle("Fluorescencia")+labs(x="Ciclo", y="Intensidad")
a<-1:50
b<-2^a
c<-1.5^a
d<-2.2^a
par(mfrow=c(1,1))
plot(x=a,y=b,type="l",xlim=c(35,50))
plot(x=a,y=c,type="l",xlim=c(35,50))
plot(x=a,y=d,type="l",xlim=c(35,50))

plot(x=a,y=b,type="l",xlim=c(35,50))
lines(x=a,y=c,col="red")

library(rdml)
install.packages("RDML",dependencies=T)
install.packages("qpcR",dependencies=T)
install.packages("lubridate",dependencies=T)
install.packages("MASS",dependencies=T)
install.packages("minpack.lm",dependencies=T)
install.packages("rgl",dependencies=T)
BiocManager::install("rgl",dependencies=T)
install.packages("robustbase",dependencies=T)

library(rgl)
library(minpack.lm)
library(MASS)
library(robustbase)
library(RDML)
library(qpcR)

data<-RDML::RDML$new("C:/Users/usuario/Downloads/qPCR_data_naive_18s_IL2.rdml")
class(data)
str(data)
summary(data)
structure(data)

RDML::ASDendrogram(data)
#Dendrogram, sample type: tissue, data types: adp, mdp

View(RDML::AsTable(data))
d_fluor<-RDML::GetFData(data)
View(d_fluor)


matplot(log((d_fluor[,2:79])),type="l")
dilutions
matplot(log(dilutions), type="l",col=c(1,1,1,2,2,2,3,3,3,4,4),lty=1)
legend("topleft",legend=c("No diluido","Diluido x5","Diluido x25","Diluido x155"),col=c(1,2,3,4),lty=1)


#cyber joins a chain of two complementary strands, all reaching the same point because the reaction saturates
abline(v=c(16,18,20,22,24),lty=2,col="gray60")
abline(h=c(2),lty=2,col="gray60")

#14,17,20,23 (18 y 22 ideal)(pendiente=eficiencia)(eficiencia teórica=2)
#14+log5(dil)


# 4assumptions:
#1. all samples must have the same efficiency
#2. The treshold must cross all lines in its exponential phase
#The ct value must be calculated from a function that adjusts to the intensity
# Must have replicas

library(qpcR)
no_dil<-qpcR::pcrfit(data=cbind(c(1:40),dilutions),
                     cyc=1,
                     fluo=2:4)
no_dil

dil_x5<-qpcR::pcrfit(data=cbind(c(1:40),dilutions),
                     cyc=1,
                     fluo=5:7)

dil_x25<-qpcR::pcrfit(data=cbind(c(1:40),dilutions),
                     cyc=1,
                     fluo=8:9)

dil_x125<-qpcR::pcrfit(data=cbind(c(1:40),dilutions),
                     cyc=1,
                     fluo=10:11)
par(mfrow=c(1,1))
ef_nod<-qpcR::efficiency(no_dil)
ef_x5<-qpcR::efficiency(dil_x5)
ef_x25<-qpcR::efficiency(dil_x25)
ef_x125<-qpcR::efficiency(dil_x125)
df<-cbind(no_dil=unlist(ef_nod),dil_x5=unlist(ef_x5),dil_x25=unlist(ef_x25),dil_x125=unlist(ef_x125))
view(df)
mean(df[1,1:4])
x<-c(0,5,25,125)
plot(x,df[1,1:4],ylim=c(0,2),type="l")


trh_no<-qpcR::efficiency(no_dil,thresh=ef_nod$cpD2)
trh_x5<-qpcR::efficiency(no_dil,thresh=ef_x5$cpD2)
trh_x25<-qpcR::efficiency(no_dil,thresh=ef_x25$cpD2)
trh_x125<-qpcR::efficiency(no_dil,thresh=ef_x125$cpD2)

df_2<-cbind(no_dil=unlist(trh_no),dilx5=unlist(trh_x5),dilx25=unlist(trh_x25),dil_x125=unlist(trh_x125))
view(df_2)
