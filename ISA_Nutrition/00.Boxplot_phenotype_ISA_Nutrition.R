#if(!require(dplyr))
#  install.packages("dplyr")
library(compiler)
library(dplyr)
#if(!require(psych))
#  install.packages("psych")
library(psych)
library(readr)
#library(ggplot2)
library(grDevices)
library(SKAT)
library(glmnet)
library(rpart)
library(rpart.plot)
library(haven)
#### LEITURA DA BASE DE DADOS DO CROMOSSOMO 1 E DO FENÓTIPO ####


fenotipos<-read.csv(file="~/Mestrado/Defesa/Dados_ISA/Bases/banco_minicurso_841.csv",header=TRUE) ## essa base será utilizada apenas para fazermos match de indivíduos entre bases, mais abaixo no match
#View(fenotipos)
#### MATRIZ DE CORRELAÇÃO ENTRE OS INDIVÍDUOS

grm_v2<-read.table("~/Mestrado/Defesa/Dados_ISA/Bases/isa.grm841x841.domicilio.grm",sep="\t")
#View(grm_v2)

## MATRIZ DE CORRELAÇÃO (PARA EXCLUIR IND MUITO CORRELACIONADOS)
relatedness_threshold <- 0.354
related_ids <- grm_v2[which((grm_v2$V4 > relatedness_threshold) & (grm_v2$V1!=grm_v2$V2)),] #ind relacionados (mesma familia)
diagonal_ids <- grm_v2[which((grm_v2$V4 > relatedness_threshold) & (grm_v2$V1==grm_v2$V2)),] #todos ind (diag princ)
unrelated_ids <- diagonal_ids[!(diagonal_ids$V2 %in% related_ids$V2), ] #base sem ind relacionados

## COMP PRINC DE ANCESTRALIDADE (ordena as duas bases pelo ID e faz merge delas)
comp<-data.frame(read_dta("~/Mestrado/Defesa/Dados_ISA/Bases/bancoISA_GWAS_n841_2024.dta")[,c(1,3,5,19,24,26,28,30,32,33)]) # dessa base vem os fenótipos
comp2<-data.frame(read_dta("~/Mestrado/Defesa/Dados_ISA/Bases/CPancestralidadeISA_2024.dta")[,c(1,4,5)]) # dessa base vem os componentes principais de ancestralidade
comp<-comp[order(comp[,1]),]
comp2<-comp2[order(comp2[,1]),]
comp_tot<-cbind(comp[,-1],comp2[,-1])
#View(comp_tot) ## merge das duas bases sem o ID agora (tem os componentes principais de ancestralidade - CP1 e CP2, sexo, idade e IMC)
#

## FAZ O MATCH DOS INDIVÍDUOS
match<-NULL
for (i in 1:nrow(fenotipos)){
  if (!is.na(fenotipos[i,13]) & !is.na(fenotipos[i,19])) {match[i]<-which(round(fenotipos[i,13],2)==round(comp_tot[,4],2) & round(fenotipos[i,19],4)==round(comp_tot[,7],4) & fenotipos[i,3]==comp_tot[,2])} else {if (!is.na(fenotipos[i,13])) match[i]<-which(round(fenotipos[i,13],2)==round(comp_tot[,4],2) & (is.na(comp_tot[,7])) & fenotipos[i,3]==comp_tot[,2]) else match[i]<-which(round(fenotipos[i,19],4)==round(comp_tot[,7],4) & (is.na(comp_tot[,4])) & fenotipos[i,3]==comp_tot[,2])}}
#	

comp_tot<-as.matrix(comp_tot[match,])
#View(comp_tot) #base final com fenótipos e com comp princ de ancestralidade e idade, sexo, etc


## PARA PUXAR O FENÓTIPO NA ORDEM CERTA
feno<-3

#### SNPS COMUNS CR 8
snp_comum<-paste("~/Mestrado/Defesa/Dados_ISA/Bases/SNPs_chr",8,"_total_imput_fam.txt",sep="")
SNPs_comuns<-read.csv(file=snp_comum, header=TRUE, sep=" ")
#View(SNPs_comuns)
names(SNPs_comuns) <- sub("_.*$", "", names(SNPs_comuns))



relevantes<- c("rs12679789","rs17601051","rs72719740","rs940079","rs79159586")
SNPs_comum_sign<- SNPs_comuns[,relevantes]

#SNPs_comum_sign<- SNPs_comuns[,grep(paste(relevantes, collapse = "|"), names(SNPs_comuns))]
#View(SNPs_comum_sign)

SNPs_comum_sign_v2<- as.matrix(SNPs_comum_sign)

##### BASE DOS SNPS RAROS DO CR 22 E SNPS COMUNS DO CR 8 JÁ ESTÃO NA MESMA ORDEM QUE FENÓTIPOS, COMP PRINCIP E MATRIZ


#
# excluindo os missings em desfecho
#
#View(comp_tot)
rownames(comp_tot) <- NULL

enableJIT(3)
for (cr in c(22)){
  # SUBSTITUIR AQUI
  
  base_snp<-paste("~/Mestrado/Defesa/Dados_ISA/Bases/SNPs_chr",cr,"_total_imput_fam.txt",sep="")
  c13_snps<-read.csv(file=base_snp, header=TRUE, sep=" ")
  #View(c13_snps)
  names(c13_snps) <- sub("_.*$", "", names(c13_snps))
  
  
  #### CLASSIFICACAO DAS BASE DE DADOS EM 0, 1 E 2 ####
  library(compiler)
  enableJIT(3)
  maff<-NULL
  for(i in 1:ncol(c13_snps)){
    n0 <- sum(c13_snps[,i]==0)
    n1 <- sum(c13_snps[,i]==1)
    n2 <- sum(c13_snps[,i]==2)
    m0 <- 2*n0+n1
    m1 <- 2*n2+n1
    maff[i]<-min(m0,m1)/(m0+m1)}
  #
  raros<-which(maff<=0.05)
  length(raros)
  SNPs_raros<-c13_snps[,raros]
  #View(SNPs_raros)
  dados_13<-as.matrix(SNPs_raros)
  #View(dados_13)
  
  #### PEGAMOS SÓ OS IND NÃO CORRELACIONADOS 
  fica<- unrelated_ids$V1
  #View(fica)
  
  Q2<-comp_tot[fica,feno] ##### é o fenótipo na escala de log, ficar atenta a isso
  dados_13<-dados_13[fica,] 
  cov_amb<-comp_tot[fica,c(1,2,10,11)] ###### checar se ficou certinho a variável de sexo, idade e componentes principais de ancestralidade 
  SNPs_comum_sign_v2<- SNPs_comum_sign_v2[fica,]
  #View(SNPs_comum_sign_v2)
  #View(cov_amb)
  
  var_comum_sign <- cbind(cov_amb,SNPs_comum_sign_v2)
  #View(var_comum_sign)
  ### aqui já estão ind indep
  #View(dados_13)
  
  #
  ##### aqui já removeu ind muito correlacionados (vamos remover só os que tem NA na resposta agora)
  ##### NÃO PRECISA MAIS DA MATRIZ DE CORRELACAO À PARTIR DAQUI
  
  excl_NA<-which(is.na(Q2)) # deve excluir menos de 10 indivíduos
  
  #View(comp_tot)
  
  ## PREPARA A BASE REMOVENDO OS INDIVÍDUOS COM MISSING NA RESPOSTA (IMC)
  Q2<-log(Q2[-excl_NA]) ##### é o fenótipo na escala de log, ficar atenta a isso
  dados_13<-dados_13[-excl_NA,] 
  var_comum_sign<-var_comum_sign[-excl_NA,] ###### checar se ficou certinho a variável de sexo, idade e componentes principais de ancestralidade 
  
  #View(dados_13)
  SNPs_raros_cr13<- dados_13 #aqui é pra ter o nome das colunas
  
  #View(SNPs_raros_cr13)
  
  SNPs_raros_cr13<-as.data.frame(SNPs_raros_cr13)
  
  mono_13<-numeric(ncol(SNPs_raros_cr13)-1)
  for (i in 1:(ncol(SNPs_raros_cr13)-1)) mono_13[i]<-length(table(SNPs_raros_cr13[,i+1]))
  checa<-table(mono_13)
  j<-1
  excl<-numeric(checa[1])
  for (i in 1:(ncol(SNPs_raros_cr13)-1)){
    if (mono_13[i]==1) {
      excl[j]<-i
      j<-j+1}}
  excl<-excl+1
  #length(excl)
  SNPs_raros_cr13<-SNPs_raros_cr13[,-excl] 
  #View(SNPs_raros_cr13)
  
  
  ### SEPARANDO OS SNPS RAROS DE ACORDO COM O GENE
  gene<- matrix(0,2,ncol(SNPs_raros_cr13))
  tab_nomes<-as.matrix(colnames(SNPs_raros_cr13)) ## PRA PEGAR SÓ OS NOMES DOS SNPS RAROS E PUXAR INFO DOS GENES DELES
  colnames(tab_nomes)<-c("SNP")
  #View(tab_nomes)
  tab_nomes<-as.data.frame(tab_nomes)
  
  #tab_nomes$SNP_original <- tab_nomes$SNP
  #tab_nomes$SNP <- sub("(_[A-Za-z]|_\\.)$", "", tab_nomes$SNP)
  
  #View(tab_nomes)
  
  
  #### ARRUMAR AQUI  - COLOCAR O NOME DA BASE QUE TEM DEPARA DE SNP E GENE
  lista_genes <- readRDS("~/Mestrado/Defesa/Dados_ISA/lista_chr_SNPs_genes_0726.rds/lista_chr_SNPs_genes_0726.rds")
  snp_info<- lista_genes[[cr]]$annot[,c(1,7,10)]
  #View(snp_info)
  snp_info$crom<- cr
  
  snp_info$gene <- ifelse(is.na(snp_info$gene_id), snp_info$nearest_gene, snp_info$gene_id)
  
  #View(snp_info)
  
  
  
  tabela_full_cr13<- tab_nomes %>% left_join(snp_info,by='SNP')
  #View(tabela_full_cr13)
  
  sum(is.na(tabela_full_cr13$gene))# ve se tem alguém sem gene marcado na base como NA
  
  remove_snps<- tabela_full_cr13$SNP[which(is.na(tabela_full_cr13$gene))]
  remove_snps
  
  idx_remove <- which(is.na(tabela_full_cr13$gene))
  
  #teste<- tabela_full_cr13[-idx_remove,c(1,5)]
  #View(teste)
  
  
  tabela_gene_cr13<- tabela_full_cr13[-idx_remove,c(1,5)] ## PEGA SÓ NOME DO SNP E GENE
  
  SNPs_raros_cr13 <- SNPs_raros_cr13[, !colnames(SNPs_raros_cr13) %in% remove_snps]
  
  #View(SNPs_raros_cr13)
  #View(tabela_gene_cr13)
  
  
  rownames(tabela_gene_cr13) <- NULL
  
  
  
  #### FENÓTIPO
  ytr_1 = Q2
}

boxplot(
  ytr_1,
  xlab="log-IMC",
  xaxt = "n",   # remove nomes padrão
  col = "lightblue"
)
