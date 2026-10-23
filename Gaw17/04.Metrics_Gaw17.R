library(dplyr)
library(readxl)
library(readr)


############# REAL
simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
#View(simula_gaw17)



snp_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/snp_info")
#View(snp_info)
snp_info<-snp_info[,4]
genes <- unique(snp_info)
colnames(genes)<- "Genes"
#View(genes)


genes_n_sign<- anti_join(genes, simula_gaw17, by = c("Genes"))

base_final_sens_grm_com_cov_cont<- NULL
base_final_esp_grm_com_cov_cont<- NULL
base_final_f1_grm_com_cov_cont<- NULL
f1_grm_com_cov_cont<- NULL

base_final_sens_gpcm_com_cov_cont<- NULL
base_final_esp_gpcm_com_cov_cont<- NULL
base_final_f1_gpcm_com_cov_cont<- NULL
f1_gpcm_com_cov_cont<- NULL

base_final_sens_grm_com_cov_categ<- NULL
base_final_esp_grm_com_cov_categ<- NULL
base_final_f1_grm_com_cov_categ<- NULL
f1_grm_com_cov_categ<- NULL

base_final_sens_gpcm_com_cov_categ<- NULL
base_final_esp_gpcm_com_cov_categ<- NULL
base_final_f1_gpcm_com_cov_categ<- NULL
f1_gpcm_com_cov_categ<- NULL


num_genes_sign_grm_com_cov_cont<- NULL

num_genes_sign_gpcm_com_cov_cont<- NULL

num_genes_sign_grm_com_cov_categ<- NULL

num_genes_sign_gpcm_com_cov_categ<- NULL


for (repli in (c(1:10))){
  ############ BANCO DE DADOS REPLICA
  gene_raro_com_cov_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/resultados_replicas_full/gene_raro_com_cov_repli",repli,".txt"))
  #View(gene_raro_com_cov_repli)
  
  
  #View(genes_n_sign)
  ############# GRM ESCORE CONT
  grm_sign_com_cov_repli_escore_cont<- gene_raro_com_cov_repli[,c(1,2,6,8)]
  #View(grm_sign_com_cov_repli_escore_cont)
  #str(grm_sign_com_cov_repli_escore_cont)
  grm_sign_com_cov_repli_escore_cont$`Valor-p GRM C COV` <-
    as.numeric(grm_sign_com_cov_repli_escore_cont$`Valor-p GRM C COV`)
  
  
  grm_sign_com_cov_repli_escore_cont <- grm_sign_com_cov_repli_escore_cont %>%
    filter(`Valor-p GRM C COV` < `Nivel de sign grm`) %>%
    arrange(`Valor-p GRM C COV`)
  grm_sign_com_cov_repli_escore_cont$replica<-c(rep(repli,nrow(grm_sign_com_cov_repli_escore_cont)))
  
  base<-cbind(repli[1],nrow(grm_sign_com_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_grm_com_cov_cont<- rbind(num_genes_sign_grm_com_cov_cont,base)
  
  especif_grm_com_cov_repli_escore_cont <- inner_join(grm_sign_com_cov_repli_escore_cont, genes_n_sign, by = c("Genes"))
  especif_grm_com_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_grm_com_cov_repli_escore_cont))/(3196)),
                                                   nrow(especif_grm_com_cov_repli_escore_cont)))
  
  
  #View(especif_grm_com_cov_repli_escore_cont)
  base_final_esp_grm_com_cov_cont<- rbind(base_final_esp_grm_com_cov_cont,especif_grm_com_cov_repli_escore_cont)
  #View(base_final_esp_grm_com_cov_cont)
  
  ## F1-SCORE
  
  
  fp_grm_com_cov_cont<- nrow(inner_join(grm_sign_com_cov_repli_escore_cont, genes_n_sign, by = c("Genes")))
  tp_grm_com_cov_cont<- nrow(inner_join(grm_sign_com_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_grm_com_cov_cont<- nrow(anti_join(simula_gaw17, grm_sign_com_cov_repli_escore_cont, by = c("Genes","Cromossomo")))
  
  f1_grm_com_cov_cont$f1<- (2*tp_grm_com_cov_cont)/((2*tp_grm_com_cov_cont) + fn_grm_com_cov_cont + fp_grm_com_cov_cont)
  
  
  base_final_f1_grm_com_cov_cont<- rbind(base_final_f1_grm_com_cov_cont,f1_grm_com_cov_cont)
  
  
  sensib_grm_com_cov_repli_escore_cont <- inner_join(grm_sign_com_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_grm_com_cov_repli_escore_cont$S<- c(rep(nrow(sensib_grm_com_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_grm_com_cov_repli_escore_cont)))
  #View(sensib_grm_com_cov_repli_escore_cont)
  
  
  base_final_sens_grm_com_cov_cont<- rbind(base_final_sens_grm_com_cov_cont,sensib_grm_com_cov_repli_escore_cont)
  #View(base_final_sens_grm_com_cov_cont)
  
  ############# GPCM ESCORE CONT
  
  gpcm_sign_com_cov_repli_escore_cont<- gene_raro_com_cov_repli[,c(1,3,7,8)]
  #View(gpcm_sign_com_cov_repli_escore_cont)
  str(gpcm_sign_com_cov_repli_escore_cont)
  
  gpcm_sign_com_cov_repli_escore_cont$`Valor-p GPCM C COV`<- as.numeric(gpcm_sign_com_cov_repli_escore_cont$`Valor-p GPCM C COV`)
  
  gpcm_sign_com_cov_repli_escore_cont <- gpcm_sign_com_cov_repli_escore_cont %>%
    filter(`Valor-p GPCM C COV` < `Nível de sign gpcm`) %>%
    arrange(`Valor-p GPCM C COV`)
  gpcm_sign_com_cov_repli_escore_cont$replica<-c(rep(repli,nrow(gpcm_sign_com_cov_repli_escore_cont)))
  #View(gpcm_sign_com_cov_repli_escore_cont)
  
  base<-cbind(repli[1],nrow(gpcm_sign_com_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_gpcm_com_cov_cont<- rbind(num_genes_sign_gpcm_com_cov_cont,base)
  
  
  especif_gpcm_com_cov_repli_escore_cont <- anti_join(gpcm_sign_com_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  especif_gpcm_com_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_gpcm_com_cov_repli_escore_cont))/(3196)),
                                                    nrow(especif_gpcm_com_cov_repli_escore_cont)))
  
  base_final_esp_gpcm_com_cov_cont<- rbind(base_final_esp_gpcm_com_cov_cont,
                                           especif_gpcm_com_cov_repli_escore_cont)
  #View(base_final_esp_gpcm_com_cov_cont)
  
  ## F1-SCORE
  
  
  fp_gpcm_com_cov_cont<- nrow(inner_join(gpcm_sign_com_cov_repli_escore_cont, genes_n_sign, by = c("Genes")))
  tp_gpcm_com_cov_cont<- nrow(inner_join(gpcm_sign_com_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_gpcm_com_cov_cont<- nrow(anti_join(simula_gaw17, gpcm_sign_com_cov_repli_escore_cont, by = c("Genes","Cromossomo")))
  
  f1_gpcm_com_cov_cont$f1<- (2*tp_gpcm_com_cov_cont)/((2*tp_gpcm_com_cov_cont) + fn_gpcm_com_cov_cont + fp_gpcm_com_cov_cont)
  
  
  base_final_f1_gpcm_com_cov_cont<- rbind(base_final_f1_gpcm_com_cov_cont,f1_gpcm_com_cov_cont)
  
  
  
  sensib_gpcm_com_cov_repli_escore_cont <- inner_join(gpcm_sign_com_cov_repli_escore_cont,
                                                       simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_gpcm_com_cov_repli_escore_cont$S<- c(rep(nrow(sensib_gpcm_com_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_gpcm_com_cov_repli_escore_cont)))
  
  base_final_sens_gpcm_com_cov_cont<- rbind(base_final_sens_gpcm_com_cov_cont,
                                            sensib_gpcm_com_cov_repli_escore_cont)
  
  #View(base_final_sens_gpcm_com_cov_cont)
  
  
  ############# GRM ESCORE CATEG
  grm_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_repli[,c(1,4,6,8)]
  grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM C COV`<- as.numeric(grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM C COV`)
  
  
  grm_sign_com_cov_repli_escore_categ <- grm_sign_com_cov_repli_escore_categ %>%
    filter(`Valor-p Arvore GRM C COV` < `Nivel de sign grm`) %>%
    arrange(`Valor-p Arvore GRM C COV`)
  #View(grm_sign_com_cov_repli_escore_categ)
  
  
  grm_sign_com_cov_repli_escore_categ$replica<-c(rep(repli,nrow(grm_sign_com_cov_repli_escore_categ)))
  #View(grm_sign_com_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(grm_sign_com_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_grm_com_cov_categ<- rbind(num_genes_sign_grm_com_cov_categ,base)
  
  
  especif_grm_com_cov_repli_escore_categ <- anti_join(grm_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  #View(especif_grm_com_cov_repli_escore_categ)
  especif_grm_com_cov_repli_escore_categ$E<- c(rep(((3196-nrow(especif_grm_com_cov_repli_escore_categ))/(3196)),
                                                    nrow(especif_grm_com_cov_repli_escore_categ)))
  #View(especif_grm_com_cov_repli_escore_categ)
  base_final_esp_grm_com_cov_categ<- rbind(base_final_esp_grm_com_cov_categ,
                                           especif_grm_com_cov_repli_escore_categ)
  #View(base_final_esp_grm_com_cov_categ)
  
  
  ## F1-SCORE
  
  
  fp_grm_com_cov_categ<- nrow(inner_join(grm_sign_com_cov_repli_escore_categ, genes_n_sign, by = c("Genes")))
  tp_grm_com_cov_categ<- nrow(inner_join(grm_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_grm_com_cov_categ<- nrow(anti_join(simula_gaw17, grm_sign_com_cov_repli_escore_categ, by = c("Genes","Cromossomo")))
  
  f1_grm_com_cov_categ$f1<- (2*tp_grm_com_cov_categ)/((2*tp_grm_com_cov_categ) + fn_grm_com_cov_categ + fp_grm_com_cov_categ)
  
  
  base_final_f1_grm_com_cov_categ<- rbind(base_final_f1_grm_com_cov_categ,f1_grm_com_cov_categ)
  
  
  
  
  
  sensib_grm_com_cov_repli_escore_categ <- inner_join(grm_sign_com_cov_repli_escore_categ,
                                                       simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_grm_com_cov_repli_escore_categ$S<- c(rep(nrow(sensib_grm_com_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_grm_com_cov_repli_escore_categ)))
  #View(sensib_grm_com_cov_repli_escore_categ)
  
  
  base_final_sens_grm_com_cov_categ<- rbind(base_final_sens_grm_com_cov_categ,sensib_grm_com_cov_repli_escore_categ)
  #View(base_final_sens_grm_com_cov_categ)
  
  
  
  ############# GPCM ESCORE CATEG
  
  gpcm_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_repli[,c(1,5,7,8)]
  #View(gpcm_sign_com_cov_repli_escore_categ)
  str(gpcm_sign_com_cov_repli_escore_categ)
  
  gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM C COV` <- as.numeric(gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM C COV`)
  
  gpcm_sign_com_cov_repli_escore_categ <- gpcm_sign_com_cov_repli_escore_categ %>%
    filter(`Valor-p Arvore GPCM C COV` < `Nível de sign gpcm`) %>%
    arrange(`Valor-p Arvore GPCM C COV`)
  gpcm_sign_com_cov_repli_escore_categ$replica<-c(rep(repli,nrow(gpcm_sign_com_cov_repli_escore_categ)))
  #View(gpcm_sign_com_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(gpcm_sign_com_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_gpcm_com_cov_categ<- rbind(num_genes_sign_gpcm_com_cov_categ,base)
  
  
  especif_gpcm_com_cov_repli_escore_categ <- anti_join(gpcm_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  #View(especif_gpcm_com_cov_repli_escore_categ)
  especif_gpcm_com_cov_repli_escore_categ$E<- (3196-nrow(especif_gpcm_com_cov_repli_escore_categ))/(3196)
  
  base_final_esp_gpcm_com_cov_categ<- rbind(base_final_esp_gpcm_com_cov_categ,
                                            especif_gpcm_com_cov_repli_escore_categ)
  
  ## F1-SCORE
  
  
  fp_gpcm_com_cov_categ<- nrow(inner_join(gpcm_sign_com_cov_repli_escore_categ, genes_n_sign, by = c("Genes")))
  tp_gpcm_com_cov_categ<- nrow(inner_join(gpcm_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_gpcm_com_cov_categ<- nrow(anti_join(simula_gaw17, gpcm_sign_com_cov_repli_escore_categ, by = c("Genes","Cromossomo")))
  
  f1_gpcm_com_cov_categ$f1<- (2*tp_gpcm_com_cov_categ)/((2*tp_gpcm_com_cov_categ) + fn_gpcm_com_cov_categ + fp_gpcm_com_cov_categ)
  
  
  base_final_f1_gpcm_com_cov_categ<- rbind(base_final_f1_gpcm_com_cov_categ,f1_gpcm_com_cov_categ)
  
  
  
  
  #View(especif_gpcm_com_cov_repli_escore_categ)
  sensib_gpcm_com_cov_repli_escore_categ <- inner_join(gpcm_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_gpcm_com_cov_repli_escore_categ$S<- c(rep(nrow(sensib_gpcm_com_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_gpcm_com_cov_repli_escore_categ)))
  
  
  base_final_sens_gpcm_com_cov_categ<- rbind(base_final_sens_gpcm_com_cov_categ,
                                            sensib_gpcm_com_cov_repli_escore_categ)
  #View(base_final_sens_gpcm_com_cov_categ)
}

#View(base_final_sens_grm_com_cov_cont)
#View(base_final_esp_grm_com_cov_cont)
#View(base_final_f1_grm_com_cov_cont)

#View(base_final_sens_gpcm_com_cov_cont)
#View(base_final_esp_gpcm_com_cov_cont)
#View(base_final_f1_gpcm_com_cov_cont)

#View(base_final_sens_grm_com_cov_categ)
#View(base_final_esp_grm_com_cov_categ)
#View(base_final_f1_grm_com_cov_categ)

#View(base_final_sens_gpcm_com_cov_categ)
#View(base_final_esp_gpcm_com_cov_categ)
#View(base_final_f1_gpcm_com_cov_categ)

#View(num_genes_sign_grm_com_cov_cont)
#View(num_genes_sign_gpcm_com_cov_cont)
#View(num_genes_sign_grm_com_cov_categ)
#View(num_genes_sign_gpcm_com_cov_categ)


##### EXPORTANDO
options(scipen = 999)

bases <- c(
  "base_final_sens_grm_com_cov_cont",
  "base_final_sens_gpcm_com_cov_cont",
  "base_final_sens_grm_com_cov_categ",
  "base_final_sens_gpcm_com_cov_categ",
  "base_final_esp_grm_com_cov_cont",
  "base_final_esp_gpcm_com_cov_cont",
  "base_final_esp_grm_com_cov_categ",
  "base_final_esp_gpcm_com_cov_categ",
  "num_genes_sign_grm_com_cov_cont",
  "num_genes_sign_gpcm_com_cov_cont",
  "num_genes_sign_grm_com_cov_categ",
  "num_genes_sign_gpcm_com_cov_categ"
)

bases_2<- c(
  "base_final_f1_grm_com_cov_cont",
  "base_final_f1_gpcm_com_cov_cont",
  "base_final_f1_grm_com_cov_categ",
  "base_final_f1_gpcm_com_cov_categ"
)

for (nome in bases_2) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0(nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}

############################ SEM COV


base_final_sens_grm_sem_cov_cont<- NULL
base_final_esp_grm_sem_cov_cont<- NULL
base_final_f1_grm_sem_cov_cont<- NULL
f1_grm_sem_cov_cont<- NULL

base_final_sens_gpcm_sem_cov_cont<- NULL
base_final_esp_gpcm_sem_cov_cont<- NULL
base_final_f1_gpcm_sem_cov_cont<- NULL
f1_gpcm_sem_cov_cont<- NULL

base_final_sens_grm_sem_cov_categ<- NULL
base_final_esp_grm_sem_cov_categ<- NULL
base_final_f1_grm_sem_cov_categ<- NULL
f1_grm_sem_cov_categ<- NULL

base_final_sens_gpcm_sem_cov_categ<- NULL
base_final_esp_gpcm_sem_cov_categ<- NULL
base_final_f1_gpcm_sem_cov_categ<- NULL
f1_gpcm_sem_cov_categ<- NULL

num_genes_sign_grm_sem_cov_cont<- NULL

num_genes_sign_gpcm_sem_cov_cont<- NULL

num_genes_sign_grm_sem_cov_categ<- NULL

num_genes_sign_gpcm_sem_cov_categ<- NULL

for (repli in (c(1:10))){
  ############ BANCO DE DADOS REPLICA
  gene_raro_sem_cov_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/resultados_replicas_full/gene_raro_sem_cov_repli",repli,".txt"))
  #View(gene_raro_sem_cov_repli)
  
  
  
  ############# GRM ESCORE CONT
  grm_sign_sem_cov_repli_escore_cont<- gene_raro_sem_cov_repli[,c(1,2,6,8)]
  #View(grm_sign_sem_cov_repli_escore_cont)
  #str(grm_sign_sem_cov_repli_escore_cont)
  grm_sign_sem_cov_repli_escore_cont$`Valor-p GRM ML` <-
    as.numeric(grm_sign_sem_cov_repli_escore_cont$`Valor-p GRM ML`)
  
  
  
  grm_sign_sem_cov_repli_escore_cont <- grm_sign_sem_cov_repli_escore_cont %>%
    filter(`Valor-p GRM ML` < `Nível sign grm`) %>%
    arrange(`Valor-p GRM ML`)
  grm_sign_sem_cov_repli_escore_cont$replica<-c(rep(repli,nrow(grm_sign_sem_cov_repli_escore_cont)))
  
  
  base<-cbind(repli[1],nrow(grm_sign_sem_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_grm_sem_cov_cont<- rbind(num_genes_sign_grm_sem_cov_cont,base)
  
  
  especif_grm_sem_cov_repli_escore_cont <- inner_join(grm_sign_sem_cov_repli_escore_cont, genes_n_sign, by = c("Genes"))
  especif_grm_sem_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_grm_sem_cov_repli_escore_cont))/(3196)),
                                                  nrow(especif_grm_sem_cov_repli_escore_cont)))
  
  #View(especif_grm_sem_cov_repli_escore_cont)
  base_final_esp_grm_sem_cov_cont<- rbind(base_final_esp_grm_sem_cov_cont,especif_grm_sem_cov_repli_escore_cont)
  #View(base_final_esp_grm_sem_cov_cont)
  
  ### F1-SCORE
  
  fp_grm_sem_cov_cont<- nrow(inner_join(grm_sign_sem_cov_repli_escore_cont, genes_n_sign, by = c("Genes")))
  tp_grm_sem_cov_cont<- nrow(inner_join(grm_sign_sem_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_grm_sem_cov_cont<- nrow(anti_join(simula_gaw17, grm_sign_sem_cov_repli_escore_cont, by = c("Genes","Cromossomo")))
  
  f1_grm_sem_cov_cont$f1<- (2*tp_grm_sem_cov_cont)/((2*tp_grm_sem_cov_cont) + fn_grm_sem_cov_cont + fp_grm_sem_cov_cont)
  
  
  base_final_f1_grm_sem_cov_cont<- rbind(base_final_f1_grm_sem_cov_cont,f1_grm_sem_cov_cont)
  
  
  
  sensib_grm_sem_cov_repli_escore_cont <- inner_join(grm_sign_sem_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_grm_sem_cov_repli_escore_cont$S<- c(rep(nrow(sensib_grm_sem_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_grm_sem_cov_repli_escore_cont)))
  #View(sensib_grm_sem_cov_repli_escore_cont)
  
  
  base_final_sens_grm_sem_cov_cont<- rbind(base_final_sens_grm_sem_cov_cont,sensib_grm_sem_cov_repli_escore_cont)
  #View(base_final_sens_grm_sem_cov_cont)
  
  ############# GPCM ESCORE CONT
  
  gpcm_sign_sem_cov_repli_escore_cont<- gene_raro_sem_cov_repli[,c(1,3,7,8)]
  #View(gpcm_sign_sem_cov_repli_escore_cont)
  str(gpcm_sign_sem_cov_repli_escore_cont)
  
  gpcm_sign_sem_cov_repli_escore_cont$`Valor-p GPCM ML`<- as.numeric(gpcm_sign_sem_cov_repli_escore_cont$`Valor-p GPCM ML`)
  
  gpcm_sign_sem_cov_repli_escore_cont <- gpcm_sign_sem_cov_repli_escore_cont %>%
    filter(`Valor-p GPCM ML` < `Nível sign gpcm`) %>%
    arrange(`Valor-p GPCM ML`)
  gpcm_sign_sem_cov_repli_escore_cont$replica<-c(rep(repli,nrow(gpcm_sign_sem_cov_repli_escore_cont)))
  #View(gpcm_sign_sem_cov_repli_escore_cont)
  
  base<-cbind(repli[1],nrow(gpcm_sign_sem_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_gpcm_sem_cov_cont<- rbind(num_genes_sign_gpcm_sem_cov_cont,base)
  
  especif_gpcm_sem_cov_repli_escore_cont <- anti_join(gpcm_sign_sem_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  especif_gpcm_sem_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_gpcm_sem_cov_repli_escore_cont))/(3196)),
                                                   nrow(especif_gpcm_sem_cov_repli_escore_cont)))
  
  base_final_esp_gpcm_sem_cov_cont<- rbind(base_final_esp_gpcm_sem_cov_cont,
                                           especif_gpcm_sem_cov_repli_escore_cont)
  #View(base_final_esp_gpcm_sem_cov_cont)
  
  
  ### F1-SCORE
  
  fp_gpcm_sem_cov_cont<- nrow(inner_join(gpcm_sign_sem_cov_repli_escore_cont, genes_n_sign, by = c("Genes")))
  tp_gpcm_sem_cov_cont<- nrow(inner_join(gpcm_sign_sem_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_gpcm_sem_cov_cont<- nrow(anti_join(simula_gaw17, gpcm_sign_sem_cov_repli_escore_cont, by = c("Genes","Cromossomo")))
  
  f1_gpcm_sem_cov_cont$f1<- (2*tp_gpcm_sem_cov_cont)/((2*tp_gpcm_sem_cov_cont) + fn_gpcm_sem_cov_cont + fp_gpcm_sem_cov_cont)
  
  
  base_final_f1_gpcm_sem_cov_cont<- rbind(base_final_f1_gpcm_sem_cov_cont,f1_gpcm_sem_cov_cont)
  
  
  
  
  sensib_gpcm_sem_cov_repli_escore_cont <- inner_join(gpcm_sign_sem_cov_repli_escore_cont,
                                                      simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_gpcm_sem_cov_repli_escore_cont$S<- c(rep(nrow(sensib_gpcm_sem_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_gpcm_sem_cov_repli_escore_cont)))
  
  base_final_sens_gpcm_sem_cov_cont<- rbind(base_final_sens_gpcm_sem_cov_cont,
                                            sensib_gpcm_sem_cov_repli_escore_cont)
  
  #View(base_final_sens_gpcm_sem_cov_cont)
  
  
  ############# GRM ESCORE CATEG
  grm_sign_sem_cov_repli_escore_categ<- gene_raro_sem_cov_repli[,c(1,4,6,8)]
  grm_sign_sem_cov_repli_escore_categ$`Valor-p Arvore GRM`<- as.numeric(grm_sign_sem_cov_repli_escore_categ$`Valor-p Arvore GRM`)
  
  
  grm_sign_sem_cov_repli_escore_categ <- grm_sign_sem_cov_repli_escore_categ %>%
    filter(`Valor-p Arvore GRM` < `Nível sign grm`) %>%
    arrange(`Valor-p Arvore GRM`)
  #View(grm_sign_sem_cov_repli_escore_categ)
  
  
  grm_sign_sem_cov_repli_escore_categ$replica<-c(rep(repli,nrow(grm_sign_sem_cov_repli_escore_categ)))
  #View(grm_sign_sem_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(grm_sign_sem_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_grm_sem_cov_categ<- rbind(num_genes_sign_grm_sem_cov_categ,base)
  
  
  
  especif_grm_sem_cov_repli_escore_categ <- anti_join(grm_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  #View(especif_grm_sem_cov_repli_escore_categ)
  especif_grm_sem_cov_repli_escore_categ$E<- c(rep(((3196-nrow(especif_grm_sem_cov_repli_escore_categ))/(3196)),
                                                   nrow(especif_grm_sem_cov_repli_escore_categ)))
  #View(especif_grm_sem_cov_repli_escore_categ)
  base_final_esp_grm_sem_cov_categ<- rbind(base_final_esp_grm_sem_cov_categ,
                                           especif_grm_sem_cov_repli_escore_categ)
  #View(base_final_esp_grm_sem_cov_categ)
  
  
  ### F1-SCORE
  
  fp_grm_sem_cov_categ<- nrow(inner_join(grm_sign_sem_cov_repli_escore_categ, genes_n_sign, by = c("Genes")))
  tp_grm_sem_cov_categ<- nrow(inner_join(grm_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_grm_sem_cov_categ<- nrow(anti_join(simula_gaw17, grm_sign_sem_cov_repli_escore_categ, by = c("Genes","Cromossomo")))
  
  f1_grm_sem_cov_categ$f1<- (2*tp_grm_sem_cov_categ)/((2*tp_grm_sem_cov_categ) + fn_grm_sem_cov_categ + fp_grm_sem_cov_categ)
  
  
  base_final_f1_grm_sem_cov_categ<- rbind(base_final_f1_grm_sem_cov_categ,f1_grm_sem_cov_categ)
  
  
  
  
  
  sensib_grm_sem_cov_repli_escore_categ <- inner_join(grm_sign_sem_cov_repli_escore_categ,
                                                      simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_grm_sem_cov_repli_escore_categ$S<- c(rep(nrow(sensib_grm_sem_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_grm_sem_cov_repli_escore_categ)))
  #View(sensib_grm_sem_cov_repli_escore_categ)
  
  
  base_final_sens_grm_sem_cov_categ<- rbind(base_final_sens_grm_sem_cov_categ,sensib_grm_sem_cov_repli_escore_categ)
  #View(base_final_sens_grm_sem_cov_categ)
  
  
  
  ############# GPCM ESCORE CATEG
  
  gpcm_sign_sem_cov_repli_escore_categ<- gene_raro_sem_cov_repli[,c(1,5,7,8)]
  #View(gpcm_sign_sem_cov_repli_escore_categ)
  str(gpcm_sign_sem_cov_repli_escore_categ)
  
  gpcm_sign_sem_cov_repli_escore_categ$`Valor-p Arvore GPCM` <- as.numeric(gpcm_sign_sem_cov_repli_escore_categ$`Valor-p Arvore GPCM`)
  
  gpcm_sign_sem_cov_repli_escore_categ <- gpcm_sign_sem_cov_repli_escore_categ %>%
    filter(`Valor-p Arvore GPCM` < `Nível sign gpcm`) %>%
    arrange(`Valor-p Arvore GPCM`)
  gpcm_sign_sem_cov_repli_escore_categ$replica<-c(rep(repli,nrow(gpcm_sign_sem_cov_repli_escore_categ)))
  #View(gpcm_sign_sem_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(gpcm_sign_sem_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_sign_gpcm_sem_cov_categ<- rbind(num_genes_sign_gpcm_sem_cov_categ,base)
  
  
  especif_gpcm_sem_cov_repli_escore_categ <- anti_join(gpcm_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  #View(especif_gpcm_sem_cov_repli_escore_categ)
  especif_gpcm_sem_cov_repli_escore_categ$E<- (3196-nrow(especif_gpcm_sem_cov_repli_escore_categ))/(3196)
  
  base_final_esp_gpcm_sem_cov_categ<- rbind(base_final_esp_gpcm_sem_cov_categ,
                                            especif_gpcm_sem_cov_repli_escore_categ)
  
  
  ### F1-SCORE
  
  fp_gpcm_sem_cov_categ<- nrow(inner_join(gpcm_sign_sem_cov_repli_escore_categ, genes_n_sign, by = c("Genes")))
  tp_gpcm_sem_cov_categ<- nrow(inner_join(gpcm_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_gpcm_sem_cov_categ<- nrow(anti_join(simula_gaw17, gpcm_sign_sem_cov_repli_escore_categ, by = c("Genes","Cromossomo")))
  
  f1_gpcm_sem_cov_categ$f1<- (2*tp_gpcm_sem_cov_categ)/((2*tp_gpcm_sem_cov_categ) + fn_gpcm_sem_cov_categ + fp_gpcm_sem_cov_categ)
  
  
  base_final_f1_gpcm_sem_cov_categ<- rbind(base_final_f1_gpcm_sem_cov_categ,f1_gpcm_sem_cov_categ)
  
  
  
  #View(especif_gpcm_sem_cov_repli_escore_categ)
  sensib_gpcm_sem_cov_repli_escore_categ <- inner_join(gpcm_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_gpcm_sem_cov_repli_escore_categ$S<- c(rep(nrow(sensib_gpcm_sem_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_gpcm_sem_cov_repli_escore_categ)))
  
  
  base_final_sens_gpcm_sem_cov_categ<- rbind(base_final_sens_gpcm_sem_cov_categ,
                                             sensib_gpcm_sem_cov_repli_escore_categ)
  #View(base_final_sens_gpcm_sem_cov_categ)
}


#View(base_final_sens_grm_sem_cov_cont)
#View(base_final_esp_grm_sem_cov_cont)
#View(base_final_f1_grm_sem_cov_cont)

#View(base_final_sens_gpcm_sem_cov_cont)
#View(base_final_esp_gpcm_sem_cov_cont)
#View(base_final_f1_gpcm_sem_cov_cont)

#View(base_final_sens_grm_sem_cov_categ)
#View(base_final_esp_grm_sem_cov_categ)
#View(base_final_f1_grm_sem_cov_categ)

#View(base_final_sens_gpcm_sem_cov_categ)
#View(base_final_esp_gpcm_sem_cov_categ)
#View(base_final_f1_gpcm_sem_cov_categ)

#View(num_genes_sign_grm_sem_cov_cont)
#View(num_genes_sign_gpcm_sem_cov_cont)
#View(num_genes_sign_grm_sem_cov_categ)
#View(num_genes_sign_gpcm_sem_cov_categ)


#####################################
options(scipen = 999)

bases <- c(
  "base_final_sens_grm_sem_cov_cont",
  "base_final_sens_gpcm_sem_cov_cont",
  "base_final_sens_grm_sem_cov_categ",
  "base_final_sens_gpcm_sem_cov_categ",
  "base_final_esp_grm_sem_cov_cont",
  "base_final_esp_gpcm_sem_cov_cont",
  "base_final_esp_grm_sem_cov_categ",
  "base_final_esp_gpcm_sem_cov_categ",
  "num_genes_sign_grm_sem_cov_cont",
  "num_genes_sign_gpcm_sem_cov_cont",
  "num_genes_sign_grm_sem_cov_categ",
  "num_genes_sign_gpcm_sem_cov_categ"
)

bases_2 <- c(
  "base_final_f1_grm_sem_cov_cont",
  "base_final_f1_gpcm_sem_cov_cont",
  "base_final_f1_grm_sem_cov_categ",
  "base_final_f1_gpcm_sem_cov_categ"
)

for (nome in bases_2) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0(nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}



##################### AGRUP COM COV

num_genes_agrup_sign_grm_com_cov_cont<- NULL
base_final_esp_grm_agrup_com_cov_cont <- NULL
base_final_sens_grm_agrup_com_cov_cont<- NULL

num_genes_agrup_sign_gpcm_com_cov_cont<- NULL
base_final_esp_gpcm_agrup_com_cov_cont <- NULL
base_final_sens_gpcm_agrup_com_cov_cont<- NULL

num_genes_agrup_sign_grm_com_cov_categ<- NULL
base_final_esp_grm_agrup_com_cov_categ <- NULL
base_final_sens_grm_agrup_com_cov_categ<- NULL

num_genes_agrup_sign_gpcm_com_cov_categ<- NULL
base_final_esp_gpcm_agrup_com_cov_categ <- NULL
base_final_sens_gpcm_agrup_com_cov_categ<- NULL


for (repli in (c(1:10))){
  ############ BANCO DE DADOS REPLICA
  gene_raro_com_cov_agrup_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/resultados_replicas_full/gene_raro_com_cov_agrup/gene_raro_com_cov_agrup_repli",repli,".txt"))
  #View(gene_raro_com_cov_repli)
  
  
  #View(genes_n_sign)
  ############# GRM ESCORE CONT
  grm_agrup_sign_com_cov_repli_escore_cont<- gene_raro_com_cov_agrup_repli[,c(1,2,6,8)]
  names(grm_agrup_sign_com_cov_repli_escore_cont)[names(grm_agrup_sign_com_cov_repli_escore_cont) == "Genes AGRUP"] <- "Genes"
  #View(grm_sign_com_cov_repli_escore_cont)
  #str(grm_sign_com_cov_repli_escore_cont)
  
  grm_agrup_sign_com_cov_repli_escore_cont$`Valor-p GRM AGRUP COV` <-
    as.numeric(grm_agrup_sign_com_cov_repli_escore_cont$`Valor-p GRM AGRUP COV`)
  
  
  grm_agrup_sign_com_cov_repli_escore_cont <- grm_agrup_sign_com_cov_repli_escore_cont %>%
    filter(`Valor-p GRM AGRUP COV` < `Nivel de signi grm`) %>%
    arrange(`Valor-p GRM AGRUP COV`)
  grm_agrup_sign_com_cov_repli_escore_cont$replica<-c(rep(repli,nrow(grm_agrup_sign_com_cov_repli_escore_cont)))
  
  base<-cbind(repli[1],nrow(grm_agrup_sign_com_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_grm_com_cov_cont<- rbind(num_genes_agrup_sign_grm_com_cov_cont,base)
  
  especif_grm_agrup_com_cov_repli_escore_cont <- inner_join(grm_agrup_sign_com_cov_repli_escore_cont, genes_n_sign, by = c("Genes"))
  especif_grm_agrup_com_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_grm_agrup_com_cov_repli_escore_cont))/(3196)),
                                                   nrow(especif_grm_agrup_com_cov_repli_escore_cont)))
  
  #View(especif_grm_com_cov_repli_escore_cont)
  base_final_esp_grm_agrup_com_cov_cont<- rbind(base_final_esp_grm_agrup_com_cov_cont,
                                          especif_grm_agrup_com_cov_repli_escore_cont)
  #View(base_final_esp_grm_com_cov_cont)
  
  
  
  sensib_agrup_grm_com_cov_repli_escore_cont <- inner_join(grm_agrup_sign_com_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_grm_com_cov_repli_escore_cont$S<- c(rep(nrow(sensib_agrup_grm_com_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_agrup_grm_com_cov_repli_escore_cont)))
  #View(sensib_grm_com_cov_repli_escore_cont)
  
  
  base_final_sens_grm_agrup_com_cov_cont<- rbind(base_final_sens_grm_agrup_com_cov_cont,sensib_agrup_grm_com_cov_repli_escore_cont)
  #View(base_final_sens_grm_agrup_com_cov_cont)
  
  ############# GPCM ESCORE CONT
  
  gpcm_agrup_sign_com_cov_repli_escore_cont<- gene_raro_com_cov_agrup_repli[,c(1,3,7,8)]
  names(gpcm_agrup_sign_com_cov_repli_escore_cont)[names(gpcm_agrup_sign_com_cov_repli_escore_cont) == "Genes AGRUP"] <- "Genes"
  
  #View(gpcm_sign_com_cov_repli_escore_cont)
  str(gpcm_agrup_sign_com_cov_repli_escore_cont)
  
  gpcm_agrup_sign_com_cov_repli_escore_cont$`Valor-p GPCM AGRUP COV`<-
    as.numeric(gpcm_agrup_sign_com_cov_repli_escore_cont$`Valor-p GPCM AGRUP COV`)
  
  gpcm_agrup_sign_com_cov_repli_escore_cont <- gpcm_agrup_sign_com_cov_repli_escore_cont %>%
    filter(`Valor-p GPCM AGRUP COV` < `Nível de signi gpcm`) %>%
    arrange(`Valor-p GPCM AGRUP COV`)
  gpcm_agrup_sign_com_cov_repli_escore_cont$replica<-c(rep(repli,nrow(gpcm_agrup_sign_com_cov_repli_escore_cont)))
  #View(gpcm_sign_com_cov_repli_escore_cont)
  
  base<-cbind(repli[1],nrow(gpcm_agrup_sign_com_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_gpcm_com_cov_cont<- rbind(num_genes_agrup_sign_gpcm_com_cov_cont,base)
  
  especif_gpcm_agrup_com_cov_repli_escore_cont <- inner_join(gpcm_agrup_sign_com_cov_repli_escore_cont, genes_n_sign, by = c("Genes"))
  especif_gpcm_agrup_com_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_gpcm_agrup_com_cov_repli_escore_cont))/(3196)),
                                                        nrow(especif_gpcm_agrup_com_cov_repli_escore_cont)))
  
  #View(especif_gpcm_com_cov_repli_escore_cont)
  base_final_esp_gpcm_agrup_com_cov_cont<- rbind(base_final_esp_gpcm_agrup_com_cov_cont,
                                                especif_gpcm_agrup_com_cov_repli_escore_cont)
  #View(base_final_esp_gpcm_com_cov_cont)
  
  
  
  sensib_agrup_gpcm_com_cov_repli_escore_cont <- inner_join(gpcm_agrup_sign_com_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_gpcm_com_cov_repli_escore_cont$S<- c(rep(nrow(sensib_agrup_gpcm_com_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_agrup_gpcm_com_cov_repli_escore_cont)))
  #View(sensib_gpcm_com_cov_repli_escore_cont)
  
  
  base_final_sens_gpcm_agrup_com_cov_cont<- rbind(base_final_sens_gpcm_agrup_com_cov_cont,sensib_agrup_gpcm_com_cov_repli_escore_cont)
  #View(base_final_sens_gpcm_agrup_com_cov_cont)  
  
  ############# GRM ESCORE CATEG
  grm_agrup_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_agrup_repli[,c(1,4,6,8)]
  names(grm_agrup_sign_com_cov_repli_escore_categ)[names(grm_agrup_sign_com_cov_repli_escore_categ) == "Genes AGRUP"] <- "Genes"
  
  grm_agrup_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM AGRUP COV`<-
    as.numeric(grm_agrup_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM AGRUP COV`)
  
  
  grm_agrup_sign_com_cov_repli_escore_categ <- grm_agrup_sign_com_cov_repli_escore_categ %>%
    filter(`Valor-p Arvore GRM AGRUP COV` < `Nivel de signi grm`) %>%
    arrange(`Valor-p Arvore GRM AGRUP COV`)
  #View(grm_sign_com_cov_repli_escore_categ)
  
  
  grm_agrup_sign_com_cov_repli_escore_categ$replica<-c(rep(repli,nrow(grm_agrup_sign_com_cov_repli_escore_categ)))
  #View(grm_sign_com_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(grm_agrup_sign_com_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_grm_com_cov_categ<- rbind(num_genes_agrup_sign_grm_com_cov_categ,base)
  
  especif_grm_agrup_com_cov_repli_escore_categ <- inner_join(grm_agrup_sign_com_cov_repli_escore_categ, genes_n_sign, by = c("Genes"))
  especif_grm_agrup_com_cov_repli_escore_categ$E<- c(rep(((3196-nrow(especif_grm_agrup_com_cov_repli_escore_categ))/(3196)),
                                                        nrow(especif_grm_agrup_com_cov_repli_escore_categ)))
  
  #View(especif_grm_com_cov_repli_escore_categ)
  base_final_esp_grm_agrup_com_cov_categ<- rbind(base_final_esp_grm_agrup_com_cov_categ,
                                                especif_grm_agrup_com_cov_repli_escore_categ)
  #View(base_final_esp_grm_com_cov_categ)
  
  
  
  sensib_agrup_grm_com_cov_repli_escore_categ <- inner_join(grm_agrup_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_grm_com_cov_repli_escore_categ$S<- c(rep(nrow(sensib_agrup_grm_com_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_agrup_grm_com_cov_repli_escore_categ)))
  #View(sensib_grm_com_cov_repli_escore_categ)
  
  
  base_final_sens_grm_agrup_com_cov_categ<- rbind(base_final_sens_grm_agrup_com_cov_categ,sensib_agrup_grm_com_cov_repli_escore_categ)
  #View(base_final_sens_grm_agrup_com_cov_categ)
  
  
  
  
  ############# GPCM ESCORE CATEG
  
  gpcm_agrup_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_agrup_repli[,c(1,5,7,8)]
  names(gpcm_agrup_sign_com_cov_repli_escore_categ)[names(gpcm_agrup_sign_com_cov_repli_escore_categ) == "Genes AGRUP"] <- "Genes"
  
  #View(gpcm_sign_com_cov_repli_escore_categ)
  str(gpcm_agrup_sign_com_cov_repli_escore_categ)
  
  gpcm_agrup_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM AGRUP COV` <-
    as.numeric(gpcm_agrup_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM AGRUP COV`)
  
  gpcm_agrup_sign_com_cov_repli_escore_categ <- gpcm_agrup_sign_com_cov_repli_escore_categ %>%
    filter(`Valor-p Arvore GPCM AGRUP COV` < `Nível de signi gpcm`) %>%
    arrange(`Valor-p Arvore GPCM AGRUP COV`)
  gpcm_agrup_sign_com_cov_repli_escore_categ$replica<-c(rep(repli,nrow(gpcm_agrup_sign_com_cov_repli_escore_categ)))
  #View(gpcm_sign_com_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(gpcm_agrup_sign_com_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_gpcm_com_cov_categ<- rbind(num_genes_agrup_sign_gpcm_com_cov_categ,base)
  
  especif_gpcm_agrup_com_cov_repli_escore_categ <- inner_join(gpcm_agrup_sign_com_cov_repli_escore_categ, genes_n_sign, by = c("Genes"))
  especif_gpcm_agrup_com_cov_repli_escore_categ$E<- c(rep(((3196-nrow(especif_gpcm_agrup_com_cov_repli_escore_categ))/(3196)),
                                                         nrow(especif_gpcm_agrup_com_cov_repli_escore_categ)))
  
  #View(especif_gpcm_com_cov_repli_escore_categ)
  base_final_esp_gpcm_agrup_com_cov_categ<- rbind(base_final_esp_gpcm_agrup_com_cov_categ,
                                                 especif_gpcm_agrup_com_cov_repli_escore_categ)
  #View(base_final_esp_gpcm_com_cov_categ)
  
  
  
  sensib_agrup_gpcm_com_cov_repli_escore_categ <- inner_join(gpcm_agrup_sign_com_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_gpcm_com_cov_repli_escore_categ$S<- c(rep(nrow(sensib_agrup_gpcm_com_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_agrup_gpcm_com_cov_repli_escore_categ)))
  #View(sensib_gpcm_com_cov_repli_escore_categ)
  
  
  base_final_sens_gpcm_agrup_com_cov_categ<- rbind(base_final_sens_gpcm_agrup_com_cov_categ,sensib_agrup_gpcm_com_cov_repli_escore_categ)
  #View(base_final_sens_gpcm_agrup_com_cov_categ)
}


View(base_final_sens_grm_agrup_com_cov_cont)
#View(base_final_esp_grm_agrup_com_cov_cont)

View(base_final_sens_gpcm_agrup_com_cov_cont)
#View(base_final_esp_gpcm_agrup_com_cov_cont)


View(base_final_sens_grm_agrup_com_cov_categ)
#View(base_final_esp_grm_agrup_com_cov_categ)

View(base_final_sens_gpcm_agrup_com_cov_categ)
#View(base_final_esp_gpcm_agrup_com_cov_categ)

View(num_genes_agrup_sign_grm_com_cov_cont)
View(num_genes_agrup_sign_gpcm_com_cov_cont)
View(num_genes_agrup_sign_grm_com_cov_categ)
View(num_genes_agrup_sign_gpcm_com_cov_categ)

#####################################
options(scipen = 999)

bases <- c(
  "base_final_sens_grm_agrup_com_cov_cont",
  "base_final_sens_gpcm_agrup_com_cov_cont",
  "base_final_sens_grm_agrup_com_cov_categ",
  "base_final_sens_gpcm_agrup_com_cov_categ",
  "num_genes_agrup_sign_grm_com_cov_cont",
  "num_genes_agrup_sign_gpcm_com_cov_cont",
  "num_genes_agrup_sign_grm_com_cov_categ",
  "num_genes_agrup_sign_gpcm_com_cov_categ"
)

for (nome in bases) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0(nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}



##################### AGRUP SEM COV
num_genes_agrup_sign_grm_sem_cov_cont<- NULL
base_final_esp_grm_agrup_sem_cov_cont <- NULL
base_final_sens_grm_agrup_sem_cov_cont<- NULL

num_genes_agrup_sign_gpcm_sem_cov_cont<- NULL
base_final_esp_gpcm_agrup_sem_cov_cont <- NULL
base_final_sens_gpcm_agrup_sem_cov_cont<- NULL

num_genes_agrup_sign_grm_sem_cov_categ<- NULL
base_final_esp_grm_agrup_sem_cov_categ <- NULL
base_final_sens_grm_agrup_sem_cov_categ<- NULL

num_genes_agrup_sign_gpcm_sem_cov_categ<- NULL
base_final_esp_gpcm_agrup_sem_cov_categ <- NULL
base_final_sens_gpcm_agrup_sem_cov_categ<- NULL


for (repli in (c(1:10))){
  ############ BANCO DE DADOS REPLICA
  gene_raro_sem_cov_agrup_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/resultados_replicas_full/gene_raro_sem_cov_agrup/gene_raro_sem_cov_agrup_repli",repli,".txt"))
  names(gene_raro_sem_cov_agrup_repli)[names(gene_raro_sem_cov_agrup_repli) == "Genes AGRUP"] <- "Genes"
  
  #View(gene_raro_com_cov_repli)
  
  #View(genes_n_sign)
  ############# GRM ESCORE CONT
  grm_agrup_sign_sem_cov_repli_escore_cont<- gene_raro_sem_cov_agrup_repli[,c(1,2,6,8)]
    #View(grm_sign_sem_cov_repli_escore_cont)
  #str(grm_sign_sem_cov_repli_escore_cont)
  
  grm_agrup_sign_sem_cov_repli_escore_cont$`Valor-p GRM ML AGRUP` <-
    as.numeric(grm_agrup_sign_sem_cov_repli_escore_cont$`Valor-p GRM ML AGRUP`)
  
  
  grm_agrup_sign_sem_cov_repli_escore_cont <- grm_agrup_sign_sem_cov_repli_escore_cont %>%
    filter(`Valor-p GRM ML AGRUP` < `Nível sign GRM`) %>%
    arrange(`Valor-p GRM ML AGRUP`)
  grm_agrup_sign_sem_cov_repli_escore_cont$replica<-c(rep(repli,nrow(grm_agrup_sign_sem_cov_repli_escore_cont)))
  
  base<-cbind(repli[1],nrow(grm_agrup_sign_sem_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_grm_sem_cov_cont<- rbind(num_genes_agrup_sign_grm_sem_cov_cont,base)
  
  especif_grm_agrup_sem_cov_repli_escore_cont <- inner_join(grm_agrup_sign_sem_cov_repli_escore_cont, genes_n_sign, by = c("Genes"))
  especif_grm_agrup_sem_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_grm_agrup_sem_cov_repli_escore_cont))/(3196)),
                                                        nrow(especif_grm_agrup_sem_cov_repli_escore_cont)))
  
  #View(especif_grm_sem_cov_repli_escore_cont)
  base_final_esp_grm_agrup_sem_cov_cont<- rbind(base_final_esp_grm_agrup_sem_cov_cont,
                                                especif_grm_agrup_sem_cov_repli_escore_cont)
  #View(base_final_esp_grm_sem_cov_cont)
  
  
  
  sensib_agrup_grm_sem_cov_repli_escore_cont <- inner_join(grm_agrup_sign_sem_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_grm_sem_cov_repli_escore_cont$S<- c(rep(nrow(sensib_agrup_grm_sem_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_agrup_grm_sem_cov_repli_escore_cont)))
  #View(sensib_grm_sem_cov_repli_escore_cont)
  
  
  base_final_sens_grm_agrup_sem_cov_cont<- rbind(base_final_sens_grm_agrup_sem_cov_cont,sensib_agrup_grm_sem_cov_repli_escore_cont)
  #View(base_final_sens_grm_agrup_sem_cov_cont)
  
  ############# GPCM ESCORE CONT
  
  gpcm_agrup_sign_sem_cov_repli_escore_cont<- gene_raro_sem_cov_agrup_repli[,c(1,3,7,8)]
  #View(gpcm_sign_sem_cov_repli_escore_cont)
  str(gpcm_agrup_sign_sem_cov_repli_escore_cont)
  
  gpcm_agrup_sign_sem_cov_repli_escore_cont$`Valor-p GPCM ML AGRUP`<-
    as.numeric(gpcm_agrup_sign_sem_cov_repli_escore_cont$`Valor-p GPCM ML AGRUP`)
  
  gpcm_agrup_sign_sem_cov_repli_escore_cont <- gpcm_agrup_sign_sem_cov_repli_escore_cont %>%
    filter(`Valor-p GPCM ML AGRUP` < `Nível sign GPCM`) %>%
    arrange(`Valor-p GPCM ML AGRUP`)
  gpcm_agrup_sign_sem_cov_repli_escore_cont$replica<-c(rep(repli,nrow(gpcm_agrup_sign_sem_cov_repli_escore_cont)))
  #View(gpcm_sign_sem_cov_repli_escore_cont)
  
  base<-cbind(repli[1],nrow(gpcm_agrup_sign_sem_cov_repli_escore_cont))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_gpcm_sem_cov_cont<- rbind(num_genes_agrup_sign_gpcm_sem_cov_cont,base)
  
  especif_gpcm_agrup_sem_cov_repli_escore_cont <- inner_join(gpcm_agrup_sign_sem_cov_repli_escore_cont, genes_n_sign, by = c("Genes"))
  especif_gpcm_agrup_sem_cov_repli_escore_cont$E<- c(rep(((3196-nrow(especif_gpcm_agrup_sem_cov_repli_escore_cont))/(3196)),
                                                         nrow(especif_gpcm_agrup_sem_cov_repli_escore_cont)))
  
  #View(especif_gpcm_sem_cov_repli_escore_cont)
  base_final_esp_gpcm_agrup_sem_cov_cont<- rbind(base_final_esp_gpcm_agrup_sem_cov_cont,
                                                 especif_gpcm_agrup_sem_cov_repli_escore_cont)
  #View(base_final_esp_gpcm_sem_cov_cont)
  
  
  
  sensib_agrup_gpcm_sem_cov_repli_escore_cont <- inner_join(gpcm_agrup_sign_sem_cov_repli_escore_cont, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_gpcm_sem_cov_repli_escore_cont$S<- c(rep(nrow(sensib_agrup_gpcm_sem_cov_repli_escore_cont)/nrow(simula_gaw17),nrow(sensib_agrup_gpcm_sem_cov_repli_escore_cont)))
  #View(sensib_gpcm_sem_cov_repli_escore_cont)
  
  
  base_final_sens_gpcm_agrup_sem_cov_cont<- rbind(base_final_sens_gpcm_agrup_sem_cov_cont,sensib_agrup_gpcm_sem_cov_repli_escore_cont)
  #View(base_final_sens_gpcm_agrup_sem_cov_cont)  
  
  ############# GRM ESCORE CATEG
  grm_agrup_sign_sem_cov_repli_escore_categ<- gene_raro_sem_cov_agrup_repli[,c(1,4,6,8)]
  grm_agrup_sign_sem_cov_repli_escore_categ$`Valor-p GRM AGRUP Arvore`<-
    as.numeric(grm_agrup_sign_sem_cov_repli_escore_categ$`Valor-p GRM AGRUP Arvore`)
  
  
  grm_agrup_sign_sem_cov_repli_escore_categ <- grm_agrup_sign_sem_cov_repli_escore_categ %>%
    filter(`Valor-p GRM AGRUP Arvore` < `Nível sign GRM`) %>%
    arrange(`Valor-p GRM AGRUP Arvore`)
  #View(grm_sign_sem_cov_repli_escore_categ)
  
  
  grm_agrup_sign_sem_cov_repli_escore_categ$replica<-c(rep(repli,nrow(grm_agrup_sign_sem_cov_repli_escore_categ)))
  #View(grm_sign_sem_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(grm_agrup_sign_sem_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_grm_sem_cov_categ<- rbind(num_genes_agrup_sign_grm_sem_cov_categ,base)
  
  especif_grm_agrup_sem_cov_repli_escore_categ <- inner_join(grm_agrup_sign_sem_cov_repli_escore_categ, genes_n_sign, by = c("Genes"))
  especif_grm_agrup_sem_cov_repli_escore_categ$E<- c(rep(((3196-nrow(especif_grm_agrup_sem_cov_repli_escore_categ))/(3196)),
                                                         nrow(especif_grm_agrup_sem_cov_repli_escore_categ)))
  
  #View(especif_grm_sem_cov_repli_escore_categ)
  base_final_esp_grm_agrup_sem_cov_categ<- rbind(base_final_esp_grm_agrup_sem_cov_categ,
                                                 especif_grm_agrup_sem_cov_repli_escore_categ)
  #View(base_final_esp_grm_sem_cov_categ)
  
  
  
  sensib_agrup_grm_sem_cov_repli_escore_categ <- inner_join(grm_agrup_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_grm_sem_cov_repli_escore_categ$S<- c(rep(nrow(sensib_agrup_grm_sem_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_agrup_grm_sem_cov_repli_escore_categ)))
  #View(sensib_grm_sem_cov_repli_escore_categ)
  
  
  base_final_sens_grm_agrup_sem_cov_categ<- rbind(base_final_sens_grm_agrup_sem_cov_categ,sensib_agrup_grm_sem_cov_repli_escore_categ)
  #View(base_final_sens_grm_agrup_sem_cov_categ)
  
  
  
  
  ############# GPCM ESCORE CATEG
  
  gpcm_agrup_sign_sem_cov_repli_escore_categ<- gene_raro_sem_cov_agrup_repli[,c(1,5,7,8)]
  #View(gpcm_sign_sem_cov_repli_escore_categ)
  str(gpcm_agrup_sign_sem_cov_repli_escore_categ)
  
  gpcm_agrup_sign_sem_cov_repli_escore_categ$`Valor-p GPCM AGRUP Arvore` <-
    as.numeric(gpcm_agrup_sign_sem_cov_repli_escore_categ$`Valor-p GPCM AGRUP Arvore`)
  
  gpcm_agrup_sign_sem_cov_repli_escore_categ <- gpcm_agrup_sign_sem_cov_repli_escore_categ %>%
    filter(`Valor-p GPCM AGRUP Arvore` < `Nível sign GPCM`) %>%
    arrange(`Valor-p GPCM AGRUP Arvore`)
  gpcm_agrup_sign_sem_cov_repli_escore_categ$replica<-c(rep(repli,nrow(gpcm_agrup_sign_sem_cov_repli_escore_categ)))
  #View(gpcm_sign_sem_cov_repli_escore_categ)
  
  base<-cbind(repli[1],nrow(gpcm_agrup_sign_sem_cov_repli_escore_categ))
  colnames(base)<- c('Replica','N_genes_sign')
  
  num_genes_agrup_sign_gpcm_sem_cov_categ<- rbind(num_genes_agrup_sign_gpcm_sem_cov_categ,base)
  
  especif_gpcm_agrup_sem_cov_repli_escore_categ <- inner_join(gpcm_agrup_sign_sem_cov_repli_escore_categ, genes_n_sign, by = c("Genes"))
  especif_gpcm_agrup_sem_cov_repli_escore_categ$E<- c(rep(((3196-nrow(especif_gpcm_agrup_sem_cov_repli_escore_categ))/(3196)),
                                                          nrow(especif_gpcm_agrup_sem_cov_repli_escore_categ)))
  
  #View(especif_gpcm_sem_cov_repli_escore_categ)
  base_final_esp_gpcm_agrup_sem_cov_categ<- rbind(base_final_esp_gpcm_agrup_sem_cov_categ,
                                                  especif_gpcm_agrup_sem_cov_repli_escore_categ)
  #View(base_final_esp_gpcm_sem_cov_categ)
  
  
  
  sensib_agrup_gpcm_sem_cov_repli_escore_categ <- inner_join(gpcm_agrup_sign_sem_cov_repli_escore_categ, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_agrup_gpcm_sem_cov_repli_escore_categ$S<- c(rep(nrow(sensib_agrup_gpcm_sem_cov_repli_escore_categ)/nrow(simula_gaw17),nrow(sensib_agrup_gpcm_sem_cov_repli_escore_categ)))
  #View(sensib_gpcm_sem_cov_repli_escore_categ)
  
  
  base_final_sens_gpcm_agrup_sem_cov_categ<- rbind(base_final_sens_gpcm_agrup_sem_cov_categ,sensib_agrup_gpcm_sem_cov_repli_escore_categ)
  #View(base_final_sens_gpcm_agrup_sem_cov_categ)
}


View(base_final_sens_grm_agrup_sem_cov_cont)
#View(base_final_esp_grm_agrup_sem_cov_cont)

View(base_final_sens_gpcm_agrup_sem_cov_cont)
#View(base_final_esp_gpcm_agrup_sem_cov_cont)


View(base_final_sens_grm_agrup_sem_cov_categ)
#View(base_final_esp_grm_agrup_sem_cov_categ)

View(base_final_sens_gpcm_agrup_sem_cov_categ)
#View(base_final_esp_gpcm_agrup_sem_cov_categ)

View(num_genes_agrup_sign_grm_sem_cov_cont)
View(num_genes_agrup_sign_gpcm_sem_cov_cont)
View(num_genes_agrup_sign_grm_sem_cov_categ)
View(num_genes_agrup_sign_gpcm_sem_cov_categ)

#####################################
options(scipen = 999)

bases <- c(
  "base_final_sens_grm_agrup_sem_cov_cont",
  "base_final_sens_gpcm_agrup_sem_cov_cont",
  "base_final_sens_grm_agrup_sem_cov_categ",
  "base_final_sens_gpcm_agrup_sem_cov_categ",
  "num_genes_agrup_sign_grm_sem_cov_cont",
  "num_genes_agrup_sign_gpcm_sem_cov_cont",
  "num_genes_agrup_sign_grm_sem_cov_categ",
  "num_genes_agrup_sign_gpcm_sem_cov_categ"
)

for (nome in bases) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0(nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}


############## METODOS TRADICIONAIS SEM COV


base_final_sens_burden_s_cov<- NULL
base_final_esp_burden_s_cov<- NULL
base_final_f1_burden_s_cov<- NULL
f1_burden_s_cov<- NULL

base_final_sens_skat_s_cov<- NULL
base_final_esp_skat_s_cov<- NULL
base_final_f1_skat_s_cov<- NULL
f1_skat_s_cov<- NULL

base_final_sens_skat_o_s_cov<- NULL
base_final_esp_skat_o_s_cov<- NULL
base_final_f1_skat_o_s_cov<- NULL
f1_skat_o_s_cov<- NULL

num_genes_sign_burden_s_cov<- NULL

num_genes_sign_skat_s_cov<- NULL

num_genes_sign_skat_o_s_cov<- NULL


for (repli in (c(1:10))){
  ############ BANCO DE DADOS BURDEN REPLICA
  burden_sem_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/result_met_trad_s_cov/resul_modelo_s_cov_burden_repli",repli,".txt"))
  
  #View(gene_raro_sem_cov_repli)
    
  
  
  burden_sem_cov_sign_repli <- burden_sem_cov_sign_repli %>%
    arrange(`Valor-p Burden`)
  
  base<-cbind(repli[1],nrow(burden_sem_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_burden_s_cov<- rbind(num_genes_sign_burden_s_cov,base)
  
  
  especif_burden_sem_cov_repli <- inner_join(burden_sem_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_burden_sem_cov_repli$E<- c(rep(((3196-nrow(especif_burden_sem_cov_repli))/(3196)),
                                                  nrow(especif_burden_sem_cov_repli)))
  
  #View(especif_burden_sem_cov_repli)
  base_final_esp_burden_s_cov<- rbind(base_final_esp_burden_s_cov,especif_burden_sem_cov_repli)
  #View(base_final_esp_burden_s_cov)
  
  ### F1-SCORE
  
  fp_burden_s_cov<- nrow(inner_join(burden_sem_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_burden_s_cov<- nrow(inner_join(burden_sem_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_burden_s_cov<- nrow(anti_join(simula_gaw17, burden_sem_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_burden_s_cov$f1<- (2*tp_burden_s_cov)/((2*tp_burden_s_cov) + fn_burden_s_cov + fp_burden_s_cov)
  
  
  base_final_f1_burden_s_cov<- rbind(base_final_f1_burden_s_cov,f1_burden_s_cov)
  
  
  
  
  sensib_burden_sem_cov_repli <- inner_join(burden_sem_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_burden_sem_cov_repli$S<- c(rep(nrow(sensib_burden_sem_cov_repli)/nrow(simula_gaw17),nrow(sensib_burden_sem_cov_repli)))
  #View(sensib_burden_sem_cov_repli)
  
  
  base_final_sens_burden_s_cov<- rbind(base_final_sens_burden_s_cov,sensib_burden_sem_cov_repli)
  #View(base_final_sens_burden_s_cov)
  
  ############# SKAT
  
  skat_sem_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/result_met_trad_s_cov/resul_modelo_s_cov_skat_repli",repli,".txt"))
  
  #View(skat_sem_cov_sign_repli)
  
  
  
  skat_sem_cov_sign_repli <- skat_sem_cov_sign_repli %>%
    arrange(`Valor-p SKAT`)
  
  base<-cbind(repli[1],nrow(skat_sem_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_skat_s_cov<- rbind(num_genes_sign_skat_s_cov,base)
  
  
  especif_skat_sem_cov_repli <- inner_join(skat_sem_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_skat_sem_cov_repli$E<- c(rep(((3196-nrow(especif_skat_sem_cov_repli))/(3196)),
                                         nrow(especif_skat_sem_cov_repli)))
  
  #View(especif_skat_sem_cov_repli)
  base_final_esp_skat_s_cov<- rbind(base_final_esp_skat_s_cov,especif_skat_sem_cov_repli)
  #View(base_final_esp_skat_s_cov)
  
  
  ### F1-SCORE
  
  fp_skat_s_cov<- nrow(inner_join(skat_sem_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_skat_s_cov<- nrow(inner_join(skat_sem_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_skat_s_cov<- nrow(anti_join(simula_gaw17, skat_sem_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_skat_s_cov$f1<- (2*tp_skat_s_cov)/((2*tp_skat_s_cov) + fn_skat_s_cov + fp_skat_s_cov)
  
  
  base_final_f1_skat_s_cov<- rbind(base_final_f1_skat_s_cov,f1_skat_s_cov)
  
  
  sensib_skat_sem_cov_repli <- inner_join(skat_sem_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_skat_sem_cov_repli$S<- c(rep(nrow(sensib_skat_sem_cov_repli)/nrow(simula_gaw17),nrow(sensib_skat_sem_cov_repli)))
  #View(sensib_skat_sem_cov_repli)
  
  
  base_final_sens_skat_s_cov<- rbind(base_final_sens_skat_s_cov,sensib_skat_sem_cov_repli)
  #View(base_final_sens_skat_s_cov)
  
  
  ############# SKAT-O
  skat_o_sem_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/result_met_trad_s_cov/resul_modelo_s_cov_skat_o_repli",repli,".txt"))
  
  #View(skat_o_sem_cov_sign_repli)
  
  
  
  skat_o_sem_cov_sign_repli <- skat_o_sem_cov_sign_repli %>%
    arrange(`Valor-p SKAT-O`)
  
  base<-cbind(repli[1],nrow(skat_o_sem_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_skat_o_s_cov<- rbind(num_genes_sign_skat_o_s_cov,base)
  
  
  especif_skat_o_sem_cov_repli <- inner_join(skat_o_sem_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_skat_o_sem_cov_repli$E<- c(rep(((3196-nrow(especif_skat_o_sem_cov_repli))/(3196)),
                                       nrow(especif_skat_o_sem_cov_repli)))
  
  #View(especif_skat_o_sem_cov_repli)
  base_final_esp_skat_o_s_cov<- rbind(base_final_esp_skat_o_s_cov,especif_skat_o_sem_cov_repli)
  #View(base_final_esp_skat_o_s_cov)
  
  ### F1-SCORE
  
  fp_skat_o_s_cov<- nrow(inner_join(skat_o_sem_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_skat_o_s_cov<- nrow(inner_join(skat_o_sem_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_skat_o_s_cov<- nrow(anti_join(simula_gaw17, skat_o_sem_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_skat_o_s_cov$f1<- (2*tp_skat_o_s_cov)/((2*tp_skat_o_s_cov) + fn_skat_o_s_cov + fp_skat_o_s_cov)
  
  
  base_final_f1_skat_o_s_cov<- rbind(base_final_f1_skat_o_s_cov,f1_skat_o_s_cov)
  
  
  sensib_skat_o_sem_cov_repli <- inner_join(skat_o_sem_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_skat_o_sem_cov_repli$S<- c(rep(nrow(sensib_skat_o_sem_cov_repli)/nrow(simula_gaw17),nrow(sensib_skat_o_sem_cov_repli)))
  #View(sensib_skat_o_sem_cov_repli)
  
  
  base_final_sens_skat_o_s_cov<- rbind(base_final_sens_skat_o_s_cov,sensib_skat_o_sem_cov_repli)
  #View(base_final_sens_skat_o_s_cov)
}



#View(base_final_sens_burden_s_cov)
#View(base_final_esp_burden_s_cov)
#View(base_final_f1_burden_s_cov)

#View(base_final_sens_skat_s_cov)
#View(base_final_esp_skat_s_cov)
#View(base_final_f1_skat_s_cov)

#View(base_final_sens_skat_o_s_cov)
#View(base_final_esp_skat_o_s_cov)
#View(base_final_f1_skat_o_s_cov)

#View(num_genes_sign_burden_s_cov)

#View(num_genes_sign_skat_s_cov)

#View(num_genes_sign_skat_o_s_cov)



#####################################
options(scipen = 999)

bases <- c(
  "base_final_sens_burden_s_cov",
  "base_final_sens_skat_s_cov",
  "base_final_sens_skat_o_s_cov",
  "base_final_esp_burden_s_cov",
  "base_final_esp_skat_s_cov",
  "base_final_esp_skat_o_s_cov",
  "num_genes_sign_burden_s_cov",
  "num_genes_sign_skat_s_cov",
  "num_genes_sign_skat_o_s_cov"
)

bases_2 <- c(
  "base_final_f1_burden_s_cov",
  "base_final_f1_skat_s_cov",
  "base_final_f1_skat_o_s_cov"
)

for (nome in bases_2) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0(nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}




############## METODOS TRADICIONAIS COM COV


base_final_sens_burden_c_cov<- NULL
base_final_esp_burden_c_cov<- NULL
base_final_f1_burden_c_cov<- NULL
f1_burden_c_cov<- NULL

base_final_sens_skat_c_cov<- NULL
base_final_esp_skat_c_cov<- NULL
base_final_f1_skat_c_cov<- NULL
f1_skat_c_cov<- NULL

base_final_sens_skat_o_c_cov<- NULL
base_final_esp_skat_o_c_cov<- NULL
base_final_f1_skat_o_c_cov<- NULL
f1_skat_o_c_cov<- NULL

num_genes_sign_burden_c_cov<- NULL

num_genes_sign_skat_c_cov<- NULL

num_genes_sign_skat_o_c_cov<- NULL


for (repli in (c(1:10))){
  ############ BANCO DE DADOS BURDEN REPLICA
  burden_com_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/result_met_trad_c_cov/resul_modelo_c_cov_burden_repli",repli,".txt"))
  
  #View(burden_com_cov_sign_repli)
  
  
  
  burden_com_cov_sign_repli <- burden_com_cov_sign_repli %>%
    arrange(`Valor-p Burden com Cov`)
  
  base<-cbind(repli[1],nrow(burden_com_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_burden_c_cov<- rbind(num_genes_sign_burden_c_cov,base)
  
  
  especif_burden_com_cov_repli <- inner_join(burden_com_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_burden_com_cov_repli$E<- c(rep(((3196-nrow(especif_burden_com_cov_repli))/(3196)),
                                         nrow(especif_burden_com_cov_repli)))
  
  #View(especif_burden_com_cov_repli)
  base_final_esp_burden_c_cov<- rbind(base_final_esp_burden_c_cov,especif_burden_com_cov_repli)
  #View(base_final_esp_burden_c_cov)
  
  
  ### F1-SCORE
  
  fp_burden_c_cov<- nrow(inner_join(burden_com_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_burden_c_cov<- nrow(inner_join(burden_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_burden_c_cov<- nrow(anti_join(simula_gaw17, burden_com_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_burden_c_cov$f1<- (2*tp_burden_c_cov)/((2*tp_burden_c_cov) + fn_burden_c_cov + fp_burden_c_cov)
  
  
  base_final_f1_burden_c_cov<- rbind(base_final_f1_burden_c_cov,f1_burden_c_cov)
  
  
  sensib_burden_com_cov_repli <- inner_join(burden_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_burden_com_cov_repli$S<- c(rep(nrow(sensib_burden_com_cov_repli)/nrow(simula_gaw17),nrow(sensib_burden_com_cov_repli)))
  #View(sensib_burden_com_cov_repli)
  
  
  base_final_sens_burden_c_cov<- rbind(base_final_sens_burden_c_cov,sensib_burden_com_cov_repli)
  #View(base_final_sens_burden_c_cov)
  
  ############# SKAT
  
  skat_com_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/result_met_trad_c_cov/resul_modelo_c_cov_skat_repli",repli,".txt"))
  
  #View(skat_com_cov_sign_repli)
  
  
  
  skat_com_cov_sign_repli <- skat_com_cov_sign_repli %>%
    arrange(`Valor-p SKAT com Cov`)
  
  base<-cbind(repli[1],nrow(skat_com_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_skat_c_cov<- rbind(num_genes_sign_skat_c_cov,base)
  
  
  especif_skat_com_cov_repli <- inner_join(skat_com_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_skat_com_cov_repli$E<- c(rep(((3196-nrow(especif_skat_com_cov_repli))/(3196)),
                                       nrow(especif_skat_com_cov_repli)))
  
  #View(especif_skat_com_cov_repli)
  base_final_esp_skat_c_cov<- rbind(base_final_esp_skat_c_cov,especif_skat_com_cov_repli)
  #View(base_final_esp_skat_c_cov)
  
  
  ### F1-SCORE
  
  fp_skat_c_cov<- nrow(inner_join(skat_com_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_skat_c_cov<- nrow(inner_join(skat_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_skat_c_cov<- nrow(anti_join(simula_gaw17, skat_com_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_skat_c_cov$f1<- (2*tp_skat_c_cov)/((2*tp_skat_c_cov) + fn_skat_c_cov + fp_skat_c_cov)
  
  
  base_final_f1_skat_c_cov<- rbind(base_final_f1_skat_c_cov,f1_skat_c_cov)
  
  
  
  sensib_skat_com_cov_repli <- inner_join(skat_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_skat_com_cov_repli$S<- c(rep(nrow(sensib_skat_com_cov_repli)/nrow(simula_gaw17),nrow(sensib_skat_com_cov_repli)))
  #View(sensib_skat_com_cov_repli)
  
  
  base_final_sens_skat_c_cov<- rbind(base_final_sens_skat_c_cov,sensib_skat_com_cov_repli)
  #View(base_final_sens_skat_c_cov)
  
  
  ############# SKAT-O
  skat_o_com_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/result_met_trad_c_cov/resul_modelo_c_cov_skat_o_repli",repli,".txt"))
  
  #View(skat_o_com_cov_sign_repli)
  
  
  
  skat_o_com_cov_sign_repli <- skat_o_com_cov_sign_repli %>%
    arrange(`Valor-p SKAT-O com Cov`)
  
  base<-cbind(repli[1],nrow(skat_o_com_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_skat_o_c_cov<- rbind(num_genes_sign_skat_o_c_cov,base)
  
  
  especif_skat_o_com_cov_repli <- inner_join(skat_o_com_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_skat_o_com_cov_repli$E<- c(rep(((3196-nrow(especif_skat_o_com_cov_repli))/(3196)),
                                         nrow(especif_skat_o_com_cov_repli)))
  
  #View(especif_skat_o_com_cov_repli)
  base_final_esp_skat_o_c_cov<- rbind(base_final_esp_skat_o_c_cov,especif_skat_o_com_cov_repli)
  #View(base_final_esp_skat_o_c_cov)
  
  ### F1-SCORE
  
  fp_skat_o_c_cov<- nrow(inner_join(skat_o_com_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_skat_o_c_cov<- nrow(inner_join(skat_o_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_skat_o_c_cov<- nrow(anti_join(simula_gaw17, skat_o_com_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_skat_o_c_cov$f1<- (2*tp_skat_o_c_cov)/((2*tp_skat_o_c_cov) + fn_skat_o_c_cov + fp_skat_o_c_cov)
  
  
  base_final_f1_skat_o_c_cov<- rbind(base_final_f1_skat_o_c_cov,f1_skat_o_c_cov)
  
  
  
  sensib_skat_o_com_cov_repli <- inner_join(skat_o_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_skat_o_com_cov_repli$S<- c(rep(nrow(sensib_skat_o_com_cov_repli)/nrow(simula_gaw17),nrow(sensib_skat_o_com_cov_repli)))
  #View(sensib_skat_o_com_cov_repli)
  
  
  base_final_sens_skat_o_c_cov<- rbind(base_final_sens_skat_o_c_cov,sensib_skat_o_com_cov_repli)
  #View(base_final_sens_skat_o_c_cov)
}



#View(base_final_sens_burden_c_cov)
#View(base_final_esp_burden_c_cov)
#View(base_final_f1_burden_c_cov)

#View(base_final_sens_skat_c_cov)
#View(base_final_esp_skat_c_cov)
#View(base_final_f1_skat_c_cov)

#View(base_final_sens_skat_o_c_cov)
#View(base_final_esp_skat_o_c_cov)
#View(base_final_f1_skat_o_c_cov)

#View(num_genes_sign_burden_c_cov)

#View(num_genes_sign_skat_c_cov)

#View(num_genes_sign_skat_o_c_cov)



#####################################
options(scipen = 999)

bases <- c(
  "base_final_sens_burden_c_cov",
  "base_final_sens_skat_c_cov",
  "base_final_sens_skat_o_c_cov",
  "base_final_esp_burden_c_cov",
  "base_final_esp_skat_c_cov",
  "base_final_esp_skat_o_c_cov",
  "num_genes_sign_burden_c_cov",
  "num_genes_sign_skat_c_cov",
  "num_genes_sign_skat_o_c_cov"
)

bases2 <- c(
  "base_final_f1_burden_c_cov",
  "base_final_f1_skat_c_cov",
  "base_final_f1_skat_o_c_cov"
)

for (nome in bases2) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0(nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}



### VERIFICANDO SE GENES IDENTIFICADOS ERRADOS ESTÃO PRÓXIMOS DOS CERTOS -  PARA OS METODOS TRAD SEM COV
metodos_trad_s_cov<- c(
  "Burden_s_cov",
  "Skat_s_cov",
  "Skat_o_s_cov"
)

for (metodo in metodos_trad_s_cov){
  simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
  #View(simula_gaw17)
  gene_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/gene_info.txt")
  
  colnames(gene_info)<- c("Genes","Cromossomo","Posicao","Gene_start","Gene_end","Gene_length","Orientation")
  
  ### gene falso positivo - metodos
  top_10_metodo <- read_excel(paste0("Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/resumo_met_trad_s_cov/",metodo,"/top10_",metodo,".xlsx"))
  #View(top_10_metodo)
  
  top_10_metodo <- anti_join(
    top_10_metodo,
    simula_gaw17 %>% select(Genes),
    by = "Genes"
  )
  top_10_metodo$Vdd_signif <- 0
  
  top_10_metodo <- left_join(
    top_10_metodo,
    gene_info %>% select(Genes, Cromossomo, Gene_start, Gene_end),
    by = "Genes"
  )
  ## PEGA SÓ QUEM TÁ EM CROMOSSOMO QUE TEM GENE VDD POSITIVO
  top10_metodo_v3<- top_10_metodo[top_10_metodo$Cromossomo %in% simula_gaw17$Cromossomo,]
  #View(top10_metodo_v3)
  simula_gaw17<- simula_gaw17[simula_gaw17$Cromossomo %in% top_10_metodo$Cromossomo,]
  #View(simula_gaw17)
  simula_gaw17 <- left_join(
    simula_gaw17,
    gene_info %>% select(Genes, Gene_start, Gene_end),
    by = c("Genes")
  )
  final<- left_join(
    simula_gaw17,
    top10_metodo_v3,
    by = c("Cromossomo")
  )
  final$distancia <- ifelse(
    !is.na(final$Gene_start.x) & !is.na(final$Gene_start.y),
    ifelse(
      final$Gene_start.x < final$Gene_start.y,
      final$Gene_end.x - final$Gene_start.y,
      final$Gene_end.y - final$Gene_start.x
    ),
    NA
  )/1e8
  
  final$proximo <- ifelse(abs(final$distancia)<0.5,1,0)
  #View(final)
  final_v2 <- final %>%
    rename_with(~ gsub("\\.x$", "_vdd_posit", .x)) %>%
    rename_with(~ gsub("\\.y$", "_ident", .x))
  #View(final_v2)
  
  write.table(
    final_v2,
    file = paste0("verif_prox_genes_errados_",metodo, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}


### VERIFICANDO SE GENES IDENTIFICADOS ERRADOS ESTÃO PRÓXIMOS DOS CERTOS -  PARA OS METODOS TRAD COM COV

metodos_trad_c_cov<- c(
  "Burden_c_cov",
  "Skat_c_cov",
  "Skat_o_c_cov"
)


for (metodo in metodos_trad_c_cov){
  simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
  #View(simula_gaw17)
  gene_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/gene_info.txt")
  
  colnames(gene_info)<- c("Genes","Cromossomo","Posicao","Gene_start","Gene_end","Gene_length","Orientation")
  
  ### gene falso positivo - metodos
  top_10_metodo <- read_excel(paste0("Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/resumo_met_trad_c_cov/",metodo,"/top10_",metodo,".xlsx"))
  #View(top_10_metodo)
  
  top_10_metodo <- anti_join(
    top_10_metodo,
    simula_gaw17 %>% select(Genes),
    by = "Genes"
  )
  top_10_metodo$Vdd_signif <- 0
  
  top_10_metodo <- left_join(
    top_10_metodo,
    gene_info %>% select(Genes, Cromossomo, Gene_start, Gene_end),
    by = "Genes"
  )
  ## PEGA SÓ QUEM TÁ EM CROMOSSOMO QUE TEM GENE VDD POSITIVO
  top10_metodo_v3<- top_10_metodo[top_10_metodo$Cromossomo %in% simula_gaw17$Cromossomo,]
  #View(top10_metodo_v3)
  simula_gaw17<- simula_gaw17[simula_gaw17$Cromossomo %in% top_10_metodo$Cromossomo,]
  #View(simula_gaw17)
  simula_gaw17 <- left_join(
    simula_gaw17,
    gene_info %>% select(Genes, Gene_start, Gene_end),
    by = c("Genes")
  )
  final<- left_join(
    simula_gaw17,
    top10_metodo_v3,
    by = c("Cromossomo")
  )
  final$distancia <- ifelse(
    !is.na(final$Gene_start.x) & !is.na(final$Gene_start.y),
    ifelse(
      final$Gene_start.x < final$Gene_start.y,
      final$Gene_end.x - final$Gene_start.y,
      final$Gene_end.y - final$Gene_start.x
    ),
    NA
  )/1e8
  
  final$proximo <- ifelse(abs(final$distancia)<0.5,1,0)
  #View(final)
  final_v2 <- final %>%
    rename_with(~ gsub("\\.x$", "_vdd_posit", .x)) %>%
    rename_with(~ gsub("\\.y$", "_ident", .x))
  #View(final_v2)
  
  write.table(
    final_v2,
    file = paste0("verif_prox_genes_errados_",metodo, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}




### VERIFICANDO SE GENES IDENTIFICADOS ERRADOS ESTÃO PRÓXIMOS DOS CERTOS -  PARA OS METODOS TRI SEM COV
metodos_tri<- c(
  "GRM_cont",
  "GRM_categ",
  "GPCM_cont",
  "GPCM_categ"
)

for (metodo in metodos_tri){
  simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
  #View(simula_gaw17)
  gene_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/gene_info.txt")
  
  colnames(gene_info)<- c("Genes","Cromossomo","Posicao","Gene_start","Gene_end","Gene_length","Orientation")
  

  ### gene falso positivo - metodos
  top_10_metodo <- read_excel(paste0("Mestrado/Defesa/resultados/resultados_replicas_full/20260419_result_replicas_s_cov/",metodo,"/top10_",metodo,"_s_cov.xlsx"))
  #View(top_10_metodo)
  
  top_10_metodo <- anti_join(
    top_10_metodo,
    simula_gaw17 %>% select(Genes),
    by = "Genes"
  )
  top_10_metodo$Vdd_signif <- 0
  
  top_10_metodo <- left_join(
    top_10_metodo,
    gene_info %>% select(Genes, Cromossomo, Gene_start, Gene_end),
    by = "Genes"
  )
  ## PEGA SÓ QUEM TÁ EM CROMOSSOMO QUE TEM GENE VDD POSITIVO
  top10_metodo_v3<- top_10_metodo[top_10_metodo$Cromossomo %in% simula_gaw17$Cromossomo,]
  #View(top10_metodo_v3)
  simula_gaw17<- simula_gaw17[simula_gaw17$Cromossomo %in% top_10_metodo$Cromossomo,]
  #View(simula_gaw17)
  simula_gaw17 <- left_join(
    simula_gaw17,
    gene_info %>% select(Genes, Gene_start, Gene_end),
    by = c("Genes")
  )
  final<- left_join(
    simula_gaw17,
    top10_metodo_v3,
    by = c("Cromossomo")
  )
  final$distancia <- ifelse(
    !is.na(final$Gene_start.x) & !is.na(final$Gene_start.y),
    ifelse(
      final$Gene_start.x < final$Gene_start.y,
      final$Gene_end.x - final$Gene_start.y,
      final$Gene_end.y - final$Gene_start.x
    ),
    NA
  )/1e8
  
  final$proximo <- ifelse(abs(final$distancia)<0.5,1,0)
  #View(final)
  final_v2 <- final %>%
    rename_with(~ gsub("\\.x$", "_vdd_posit", .x)) %>%
    rename_with(~ gsub("\\.y$", "_ident", .x))
  #View(final_v2)
  
  write.table(
    final_v2,
    file = paste0("verif_prox_genes_errados_",metodo, "_s_cov.csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}



### VERIFICANDO SE GENES IDENTIFICADOS ERRADOS ESTÃO PRÓXIMOS DOS CERTOS -  PARA OS METODOS TRI COM COV
for (metodo in metodos_tri){
  simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
  #View(simula_gaw17)
  gene_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/gene_info.txt")
  
  colnames(gene_info)<- c("Genes","Cromossomo","Posicao","Gene_start","Gene_end","Gene_length","Orientation")
  
  
  ### gene falso positivo - metodos
  top_10_metodo <- read_excel(paste0("Mestrado/Defesa/resultados/resultados_replicas_full/20260419_result_replicas_c_cov/",metodo,"/top10_",metodo,"_c_cov.xlsx"))
  #View(top_10_metodo)
  
  top_10_metodo <- anti_join(
    top_10_metodo,
    simula_gaw17 %>% select(Genes),
    by = "Genes"
  )
  top_10_metodo$Vdd_signif <- 0
  
  top_10_metodo <- left_join(
    top_10_metodo,
    gene_info %>% select(Genes, Cromossomo, Gene_start, Gene_end),
    by = "Genes"
  )
  ## PEGA SÓ QUEM TÁ EM CROMOSSOMO QUE TEM GENE VDD POSITIVO
  top10_metodo_v3<- top_10_metodo[top_10_metodo$Cromossomo %in% simula_gaw17$Cromossomo,]
  #View(top10_metodo_v3)
  simula_gaw17<- simula_gaw17[simula_gaw17$Cromossomo %in% top_10_metodo$Cromossomo,]
  #View(simula_gaw17)
  simula_gaw17 <- left_join(
    simula_gaw17,
    gene_info %>% select(Genes, Gene_start, Gene_end),
    by = c("Genes")
  )
  final<- left_join(
    simula_gaw17,
    top10_metodo_v3,
    by = c("Cromossomo")
  )
  final$distancia <- ifelse(
    !is.na(final$Gene_start.x) & !is.na(final$Gene_start.y),
    ifelse(
      final$Gene_start.x < final$Gene_start.y,
      final$Gene_end.x - final$Gene_start.y,
      final$Gene_end.y - final$Gene_start.x
    ),
    NA
  )/1e8
  
  final$proximo <- ifelse(abs(final$distancia)<0.5,1,0)
  #View(final)
  final_v2 <- final %>%
    rename_with(~ gsub("\\.x$", "_vdd_posit", .x)) %>%
    rename_with(~ gsub("\\.y$", "_ident", .x))
  #View(final_v2)
  
  write.table(
    final_v2,
    file = paste0("verif_prox_genes_errados_",metodo, "_c_cov.csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}



# BOXPLOT SENSIBILIDADE S COV

sensibilid_s_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/sensibilidade_por_metodo_s_cov.xlsx")
#View(sensibilid_s_cov)

boxplot(
  sensibilid_s_cov,
  xaxt = "n",   # remove nomes padrão
  col = "#FF6699",
  ylim=c(0.0,0.7)
)

text(
  x = 1:ncol(sensibilid_s_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(sensibilid_s_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)


#apply(sensibilid_s_cov, 2, mean)
#apply(sensibilid_s_cov, 2, sd)

# BOXPLOT ESPECIFICIDADE S COV
especif_s_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/especificidade_por_metodo_s_cov.xlsx")
#View(sensibilid_c_cov)

boxplot(
  especif_s_cov,
  xaxt = "n",   # remove nomes padrão
  col = '#FF7538',
  ylim=c(0.85,1)
  #main='Especificidade por método sem covariável'
)

text(
  x = 1:ncol(especif_s_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(especif_s_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
#apply(especif_s_cov, 2, mean)
#apply(especif_s_cov, 2, sd)


# NUMERO DE GENES IDENTIFIC MODELO S COV

n_genes_s_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/n_genes_signif_por_metodo_s_cov.xlsx")
#View(sensibilid_c_cov)

boxplot(
  n_genes_s_cov,
  xaxt = "n",   # remove nomes padrão
  col = '#CC66FF',
  ylim=c(0,200)
  #main='Quantidade de genes identificados por método sem covariável'
)

text(
  x = 1:ncol(n_genes_s_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(n_genes_s_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
#apply(n_genes_s_cov, 2, mean)
#apply(n_genes_s_cov, 2, sd)


# BOXPLOT F1 SCORE

f1_s_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/f1_escore_por_metodo_s_cov.xlsx")
#View(sensibilid_s_cov)

boxplot(
  f1_s_cov,
  xaxt = "n",   # remove nomes padrão
  col = "#2A9D8F",
  ylim=c(0.0,0.7)
)

text(
  x = 1:ncol(f1_s_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(f1_s_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)



# BOXPLOT SENSIBILIDADE C COV
sensibilid_c_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/sensibilidade_por_metodo_c_cov.xlsx")
#View(sensibilid_c_cov)

boxplot(
  sensibilid_c_cov,
  xaxt = "n",   # remove nomes padrão
  col = "#FF6699",
  ylim=c(0.0,0.7)
  #main='Sensibilidade por método com covariável'
)

text(
  x = 1:ncol(sensibilid_c_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(sensibilid_c_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
summary(sensibilid_c_cov)

# BOXPLOT ESPECIFICIDADE C COV
especif_c_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/especificidade_por_metodo_c_cov.xlsx")
#View(sensibilid_c_cov)

boxplot(
  especif_c_cov,
  xaxt = "n",   # remove nomes padrão
  col = '#FF7538',
  ylim=c(0.85,1)
  #main='Especificidade por método com covariável'
)

text(
  x = 1:ncol(especif_c_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(especif_c_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
#summary(especif_c_cov)

# NUMERO DE GENES IDENTIFIC MODELO C COV

n_genes_c_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/n_genes_signif_por_metodo_c_cov.xlsx")
#View(sensibilid_c_cov)

boxplot(
  n_genes_c_cov,
  xaxt = "n",   # remove nomes padrão
  col = '#CC66FF',
  ylim=c(0,250)#main='Quantidade de genes identificados por método sem covariável'
)

text(
  x = 1:ncol(n_genes_c_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(n_genes_c_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
#summary(n_genes_c_cov)



# BOXPLOT F1 SCORE

f1_c_cov<- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260429_resumo_result_replicas/f1_escore_por_metodo_c_cov.xlsx")
#View(sensibilid_s_cov)

boxplot(
  f1_c_cov,
  xaxt = "n",   # remove nomes padrão
  col = "#2A9D8F",
  ylim=c(0.0,0.7)
)

text(
  x = 1:ncol(f1_c_cov),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(f1_c_cov),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
