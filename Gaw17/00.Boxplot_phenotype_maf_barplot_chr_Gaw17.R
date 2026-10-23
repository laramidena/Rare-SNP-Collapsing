if(!require(dplyr))
  install.packages("dplyr")
library(dplyr)
if(!require(psych))
  install.packages("psych")
library(psych)
library(readr)
library(ggplot2)
library(grDevices)

#### LEITURA DA BASE DE DADOS DO CROMOSSOMO 1 E DO FENÓTIPO ####
snp_info <- read_csv("TG 1/Dados_independentes/Dados_independentes/snp_info")
unr_phen <- read_csv("TG 1/Dados_independentes/Dados_independentes/unr_phen.1")
c13_snps <- read_csv("TG 1/Dados_independentes/Dados_independentes/c13_snps.unr")

#View(c13_snps)

#boxplot(unr_phen$Q1, col = "lightblue", xlab = "Q1")
#boxplot.stats(unr_phen$Q1)$out
############# RETIRANDO SNPS SEM VARIAÇÃO NO CROMOSSOMO 13
mono_13<-numeric(ncol(c13_snps)-1)
for (i in 1:(ncol(c13_snps)-1)) mono_13[i]<-length(table(c13_snps[,i+1]))
checa<-table(mono_13)
j<-1
excl<-numeric(checa[1])
for (i in 1:(ncol(c13_snps)-1)){
  if (mono_13[i]==1) {
    excl[j]<-i
    j<-j+1}}
excl<-excl+1
c13_snps<-c13_snps[,-excl] 
View(c13_snps)
ncol(c13_snps) #NÃO REMOVEU SNPs SEM VARIAÇÃO, MAS AQUI JÁ TÁ SEM A COLUNA ID



#### CLASSIFICACAO DAS BASE DE DADOS EM 0, 1 E 2 ####
dados_13<-matrix(0,nrow=nrow(c13_snps),ncol = ncol(c13_snps))
for (i in 1:nrow(dados_13)){
  for (j in 1: ncol(dados_13)){
    if (c13_snps[i,j]=="A/A"|c13_snps[i,j]=="T/T") dados_13[i,j]=2
    if (c13_snps[i,j]=="G/G"|c13_snps[i,j]=="C/C") dados_13[i,j]=0
    if (c13_snps[i,j]!="A/A" & c13_snps[i,j]!="T/T" &
        c13_snps[i,j]!="G/G" & c13_snps[i,j]!="C/C") dados_13[i,j]=1
  }
}
View(dados_13)
ncol(dados_13)
table(dados_13[,7])

#### EXPORTA O BANCO DE DADOS PADRONIZADO SEM O IDENTIFICADOR
#write.csv(dados_13, "/dados_cr_13_rev.csv", row.names = FALSE) 
#write.table(dados_2,"dados_c1.txt",sep=",")
nomes_13 = colnames(c13_snps)
View(nomes_13)
colnames(dados_13)<- nomes_13


#### VENDO A FREQ ALÉLICA DE TDS OS SNPS DO CR13
freq_alelica_SNPs_cr13 <-matrix(0,nrow=2,ncol = ncol(c13_snps))
row.names(freq_alelica_SNPs_cr13)<- c('freq alelo A','freq alelo a')
colnames(freq_alelica_SNPs_cr13)<- colnames(c13_snps)
View(freq_alelica_SNPs_cr13)
for (i in 1:ncol(dados_13)){
  freq_alelica_SNPs_cr13[1,i]<- (0*length(dados_13[which(dados_13[,i]=='2'),i])/(2*nrow(dados_13))) + 
    (1*length(dados_13[which(dados_13[,i]=='1'),i])/(2*nrow(dados_13))) +
    (2*length(dados_13[which(dados_13[,i]=='0'),i])/(2*nrow(dados_13)))
  
  freq_alelica_SNPs_cr13[2,i]<- (2*length(dados_13[which(dados_13[,i]=='2'),i])/(2*nrow(dados_13))) + 
    (1*length(dados_13[which(dados_13[,i]=='1'),i])/(2*nrow(dados_13))) +
    (0*length(dados_13[which(dados_13[,i]=='0'),i])/(2*nrow(dados_13)))
  
}
View(freq_alelica_SNPs_cr13)

#separando os SNPs em raros (MAF<5%) e comuns e colocando os nomes deles

indice_SNPs_raros_cr_13<-which(freq_alelica_SNPs_cr13[,1:ncol(dados_13)]<0.05)
for (i in 1:length(indice_SNPs_raros_cr_13)){
  if (indice_SNPs_raros_cr_13[i]%%2==0) indice_SNPs_raros_cr_13[i] = indice_SNPs_raros_cr_13[i]/2
  else indice_SNPs_raros_cr_13[i] = (indice_SNPs_raros_cr_13[i]+1)/2
}
indice_SNPs_raros_cr_13

SNPs_raros_cr13<- dados_13[,indice_SNPs_raros_cr_13]
View(SNPs_raros_cr13)
colnames(SNPs_raros_cr13)<- nomes_13[indice_SNPs_raros_cr_13]
SNPs_comuns_cr13<- dados_13[,-c(indice_SNPs_raros_cr_13)]
View(SNPs_comuns_cr13)

### FAZENDO O MAF DE CADA SNP DO CR 13 PRA FAZER BOXPLOT
maf_13<- NULL
for (i in 1:ncol(freq_alelica_SNPs_cr13)){
  maf_13[i]<- min(freq_alelica_SNPs_cr13[1,i],freq_alelica_SNPs_cr13[2,i])
}
maf_13


######## FENÓTIPO
ytr_1 = NULL #cria uma lista vazia

for (i in 1:nrow(unr_phen)){
  ytr_1[i] = unr_phen$Q1[i]
}
ytr_1


######## JUNTANDO AS INFORMAÇÕES DE TODOS OS 22 CROMOSSOMOS 


###### ANÁLISE DESCRITIVA

#### AGREGANDO A INFORMAÇÃO DE QUAL O CROSSOMO PARA CADA FREQUÊNCIA DO ALELO MENOR

### PARTE 1
MAF_1<- cbind(maf_1, rep(1,length(maf_1)))
MAF_2<- cbind(maf_2, rep(2,length(maf_2)))
MAF_3<- cbind(maf_3, rep(3,length(maf_3)))
MAF_4<- cbind(maf_4, rep(4,length(maf_4)))
MAF_5<- cbind(maf_5, rep(5,length(maf_5)))
MAF_6<- cbind(maf_6, rep(6,length(maf_6)))
MAF_7<- cbind(maf_7, rep(7,length(maf_7)))
MAF_8<- cbind(maf_8, rep(8,length(maf_8)))
MAF_9<- cbind(maf_9, rep(9,length(maf_9)))
MAF_10<- cbind(maf_10, rep(10,length(maf_10)))
MAF_11<- cbind(maf_11, rep(11,length(maf_11)))
MAF_tds<- rbind(MAF_1,MAF_2,MAF_3,MAF_4,MAF_5,MAF_6,
                MAF_7,MAF_8,MAF_9,MAF_10,MAF_11)
MAF_tds<- as.data.frame(MAF_tds)
colnames(MAF_tds)<- c('MAF','Cr')
red_to_lilac_palette_1 <- colorRampPalette(c("lightcoral","orchid4"))(12)

### CONSTRUINDO BOXPLOT PARA O MAF DOS CROMOSSOMOS DE 1 A 11
boxplot(MAF_tds$MAF~MAF_tds$Cr, xlab = "Cromossomo", ylab= "MAF",
        col=red_to_lilac_palette_1)

### CONSTRUINDO GRÁFICO DE BARRAS PARA OS CROMOSSOMOS DE 1 A 11
SNPs_1<- rbind(ncol(SNPs_comuns_cr1),ncol(SNPs_raros_cr1))
SNPs_2<- rbind(ncol(SNPs_comuns_cr2),ncol(SNPs_raros_cr2))
SNPs_3<- rbind(ncol(SNPs_comuns_cr3),ncol(SNPs_raros_cr3))
SNPs_4<- rbind(ncol(SNPs_comuns_cr4),ncol(SNPs_raros_cr4))
SNPs_5<- rbind(ncol(SNPs_comuns_cr5),ncol(SNPs_raros_cr5))
SNPs_6<- rbind(ncol(SNPs_comuns_cr6),ncol(SNPs_raros_cr6))
SNPs_7<- rbind(ncol(SNPs_comuns_cr7),ncol(SNPs_raros_cr7))
SNPs_8<- rbind(ncol(SNPs_comuns_cr8),ncol(SNPs_raros_cr8))
SNPs_9<- rbind(ncol(SNPs_comuns_cr9),ncol(SNPs_raros_cr9))
SNPs_10<- rbind(ncol(SNPs_comuns_cr10),ncol(SNPs_raros_cr10))
SNPs_11<- rbind(ncol(SNPs_comuns_cr11),ncol(SNPs_raros_cr11))
SNPs_10<- rbind(ncol(SNPs_comuns_cr10),ncol(SNPs_raros_cr10))
SNPs_11<- rbind(ncol(SNPs_comuns_cr11),ncol(SNPs_raros_cr11))
SNPs<- cbind(SNPs_1,SNPs_2,SNPs_3,SNPs_4,SNPs_5,SNPs_6,
             SNPs_7,SNPs_8,SNPs_9,SNPs_10,SNPs_11)
rownames(SNPs)<- c("SNPs comuns", "SNPs raros")
colnames(SNPs)<- c("Cr 1", "Cr 2","Cr 3", "Cr 4", "Cr 5","Cr 6",
                   "Cr 7","Cr 8","Cr 9","Cr 10","Cr 11")
bp<- barplot(SNPs,horiz=F,beside=T, col=c("#a1e9f0", "#d9b1f0"),ylim=c(0,2500),
             xlab="Cromossomo",ylab="Frequência absoluta")
legend("topright",inset = c(-0.13, -0.115),"#a1e9f0", "#d9b1f0",
       legend=c("SNPs comuns","SNPs raros"), fill = c("#a1e9f0", "#d9b1f0"),bty="n")
text(x=c(bp), y=SNPs,labels=paste0(round(prop.table(SNPs, margin = 2) * 100,2),"%")
     ,pos=3)

### PARTE 2
MAF_12<- cbind(maf_12, rep(12,length(maf_12)))
MAF_13<- cbind(maf_13, rep(13,length(maf_13)))
MAF_14<- cbind(maf_14, rep(14,length(maf_14)))
MAF_15<- cbind(maf_15, rep(15,length(maf_15)))
MAF_16<- cbind(maf_16, rep(16,length(maf_16)))
MAF_17<- cbind(maf_17, rep(17,length(maf_17)))
MAF_18<- cbind(maf_18, rep(18,length(maf_18)))
MAF_19<- cbind(maf_19, rep(19,length(maf_19)))
MAF_20<- cbind(maf_20, rep(20,length(maf_20)))
MAF_21<- cbind(maf_21, rep(21,length(maf_21)))
MAF_22<- cbind(maf_22, rep(22,length(maf_22)))
MAF_tds_2<- rbind(MAF_12,MAF_13,MAF_14,MAF_15,MAF_16,
                  MAF_17,MAF_18,MAF_19,MAF_20,MAF_21,MAF_22)
MAF_tds_2<- as.data.frame(MAF_tds_2)
colnames(MAF_tds_2)<- c('MAF','Cr')


### CONSTRUINDO BOXPLOT PARA O MAF PARA OS CROMOSSOMOS DE 12 A 22
red_to_lilac_palette <- colorRampPalette(c("paleturquoise4","palegreen" ))(12)
boxplot(MAF_tds_2$MAF~MAF_tds_2$Cr, xlab = "Cromossomo", ylab= "MAF",
        col=red_to_lilac_palette)


### CONSTRUINDO GRÁFICO DE BARRAS PARA OS CROMOSSOMOS DE 12 A 22
SNPs_12<- rbind(ncol(SNPs_comuns_cr12),ncol(SNPs_raros_cr12))
SNPs_13<- rbind(ncol(SNPs_comuns_cr13),ncol(SNPs_raros_cr13))
SNPs_14<- rbind(ncol(SNPs_comuns_cr14),ncol(SNPs_raros_cr14))
SNPs_15<- rbind(ncol(SNPs_comuns_cr15),ncol(SNPs_raros_cr15))
SNPs_16<- rbind(ncol(SNPs_comuns_cr16),ncol(SNPs_raros_cr16))
SNPs_17<- rbind(ncol(SNPs_comuns_cr17),ncol(SNPs_raros_cr17))
SNPs_18<- rbind(ncol(SNPs_comuns_cr18),ncol(SNPs_raros_cr18))
SNPs_19<- rbind(ncol(SNPs_comuns_cr19),ncol(SNPs_raros_cr19))
SNPs_20<- rbind(ncol(SNPs_comuns_cr20),ncol(SNPs_raros_cr20))
SNPs_21<- rbind(ncol(SNPs_comuns_cr21),ncol(SNPs_raros_cr21))
SNPs_22<- rbind(ncol(SNPs_comuns_cr22),ncol(SNPs_raros_cr22))
SNPs_v2<- cbind(SNPs_12,SNPs_13,SNPs_14,SNPs_15,SNPs_16,
                SNPs_17,SNPs_18,SNPs_19,SNPs_20,SNPs_21,SNPs_22)
rownames(SNPs_v2)<- c("SNPs comuns", "SNPs raros")
colnames(SNPs_v2)<- c("Cr 12","Cr 13", "Cr 14", "Cr 15","Cr 16",
                      "Cr 17","Cr 18","Cr 19","Cr 20","Cr 21","Cr 22")
bp_2<- barplot(SNPs_v2,horiz=F,beside=T, col=c("#a1e9f0", "#d9b1f0"),ylim=c(0,2500),
               xlab="Cromossomo",ylab="Frequência absoluta")
legend("topright",inset = c(-0.13, -0.115),"#a1e9f0", "#d9b1f0", 
       legend=c("SNPs comuns","SNPs raros"), fill = c("#a1e9f0", "#d9b1f0"),bty="n")
text(x=c(bp_2), y=SNPs_v2,labels=paste0(round(prop.table(SNPs_v2, margin = 2) *
                                                100,2),"%"),pos=3)


