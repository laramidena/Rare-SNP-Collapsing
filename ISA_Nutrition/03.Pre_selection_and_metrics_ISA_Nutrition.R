library(readr)
library(dplyr)

############## METODOS TRADICIONAIS COM COV

num_genes_sign_burden_c_cov<- NULL

num_genes_sign_skat_c_cov<- NULL

num_genes_sign_skat_o_c_cov<- NULL


resul_modelo_c_cov_burden_dados_ISA <- 

  ############ BANCO DE DADOS BURDEN REPLICA
burden_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_burden_dados_ISA.txt")

View(resul_modelo_c_cov_burden_dados_ISA)


burden_com_cov_sign_repli <- burden_com_cov_sign_repli %>%
  arrange(`Valor-p Burden com Cov`)

base<-cbind('Burden',nrow(burden_com_cov_sign_repli))
colnames(base)<- c('Metodo','N_genes_sign')

############# SKAT

skat_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_dados_ISA.txt")


View(skat_com_cov_sign_repli)



skat_com_cov_sign_repli <- skat_com_cov_sign_repli %>%
  arrange(`Valor-p SKAT com Cov`)

base_2<-cbind('SKAT',nrow(skat_com_cov_sign_repli))
colnames(base_2)<- c('Metodo','N_genes_sign')



############# SKAT-O
skat_o_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_o_dados_ISA.txt")

View(skat_o_com_cov_sign_repli)



skat_o_com_cov_sign_repli <- skat_o_com_cov_sign_repli %>%
  arrange(`Valor-p SKAT-O com Cov`)

base_3<-cbind('SKAT-O',nrow(skat_o_com_cov_sign_repli))
colnames(base_3)<- c('Metodo','N_genes_sign')

final<- rbind(base,base_2,base_3)
num_genes_sign_met_trad_dados_ISA<- final


options(scipen = 999)

bases <- c(
  "num_genes_sign_met_trad_dados_ISA"
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




############# RESULTADOS TRI - INDIVIDUAIS


num_genes_sign_grm_com_cov_categ<- NULL
num_genes_sign_gpcm_com_cov_categ<- NULL


gene_raro_com_cov_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_tri_dados_ISA/gene_raro_com_cov_TRI_dados_ISA_v2.txt")
#View(gene_raro_com_cov_repli)




#View(base_final_sens_gpcm_com_cov_cont)


############# GRM ESCORE CATEG -- 

grm_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_repli[,c(1,4,6,8)]
View(grm_sign_com_cov_repli_escore_categ)

grm_sign_com_cov_repli_escore_categ <-
  grm_sign_com_cov_repli_escore_categ[
    grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM C COV` != "sem corte",
  ]


grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM C COV`<- as.numeric(grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM C COV`)

grm_sign_com_cov_repli_escore_categ <- grm_sign_com_cov_repli_escore_categ %>%
  filter(`Valor-p Arvore GRM C COV` < `Nivel de sign grm`) %>%
  arrange(`Valor-p Arvore GRM C COV`)




base<-cbind('GRM categorizado com cov',nrow(grm_sign_com_cov_repli_escore_categ))
colnames(base)<- c('Metodo','N_genes_sign')
View(base)

genes_sign_grm_com_cov_categ_ISA_final<- grm_sign_com_cov_repli_escore_categ

View(genes_sign_grm_com_cov_categ_ISA_final)


############# GPCM ESCORE CATEG

gpcm_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_repli[,c(1,5,7,8)]
#View(gpcm_sign_com_cov_repli_escore_categ)


gpcm_sign_com_cov_repli_escore_categ <-
  gpcm_sign_com_cov_repli_escore_categ[
    gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM C COV` != "sem corte",
  ]



gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM C COV` <- as.numeric(gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM C COV`)

gpcm_sign_com_cov_repli_escore_categ <- gpcm_sign_com_cov_repli_escore_categ %>%
  filter(`Valor-p Arvore GPCM C COV` < `Nível de sign gpcm`) %>%
  arrange(`Valor-p Arvore GPCM C COV`)




base2<-cbind('GPCM categorizado com cov',nrow(gpcm_sign_com_cov_repli_escore_categ))
colnames(base2)<- c('Metodo','N_genes_sign')





genes_sign_gpcm_com_cov_categ_ISA_final <- gpcm_sign_com_cov_repli_escore_categ
View(genes_sign_gpcm_com_cov_categ_ISA_final)

num_genes_sign_tri_dados_ISA<- rbind(base,base2)
View(num_genes_sign_tri_dados_ISA)

##### EXPORTANDO
options(scipen = 999)

bases <- c(
  "num_genes_sign_tri_dados_ISA",
  "genes_sign_grm_com_cov_categ_ISA_final",
  "genes_sign_gpcm_com_cov_categ_ISA_final"
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





############# RESULTADOS TRI - AGRUPADOS


num_genes_sign_grm_com_cov_categ<- NULL
num_genes_sign_gpcm_com_cov_categ<- NULL


gene_raro_com_cov_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_tri_dados_ISA/gene_raro_com_cov_agrup_dados_ISA.txt")
#View(gene_raro_com_cov_repli)




#View(base_final_sens_gpcm_com_cov_cont)


############# GRM ESCORE CATEG -- 

grm_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_repli[,c(1,4,6,8)]
View(grm_sign_com_cov_repli_escore_categ)

grm_sign_com_cov_repli_escore_categ <-
  grm_sign_com_cov_repli_escore_categ[
    grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM AGRUP COV` != "sem corte",
  ]


grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM AGRUP COV`<- as.numeric(grm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GRM AGRUP COV`)

grm_sign_com_cov_repli_escore_categ <- grm_sign_com_cov_repli_escore_categ %>%
  filter(`Valor-p Arvore GRM AGRUP COV` < `Nivel de signi grm`) %>%
  arrange(`Valor-p Arvore GRM AGRUP COV`)




base<-cbind('GRM categorizado com cov AGRUP',nrow(grm_sign_com_cov_repli_escore_categ))
colnames(base)<- c('Metodo agrup','N_genes_sign')
View(base)

genes_sign_grm_com_cov_categ_ISA_agrup_final<- grm_sign_com_cov_repli_escore_categ

View(genes_sign_grm_com_cov_categ_ISA_agrup_final)


############# GPCM ESCORE CATEG

gpcm_sign_com_cov_repli_escore_categ<- gene_raro_com_cov_repli[,c(1,5,7,8)]
#View(gpcm_sign_com_cov_repli_escore_categ)


gpcm_sign_com_cov_repli_escore_categ <-
  gpcm_sign_com_cov_repli_escore_categ[
    gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM AGRUP COV` != "sem corte",
  ]



gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM AGRUP COV` <- as.numeric(gpcm_sign_com_cov_repli_escore_categ$`Valor-p Arvore GPCM AGRUP COV`)

gpcm_sign_com_cov_repli_escore_categ <- gpcm_sign_com_cov_repli_escore_categ %>%
  filter(`Valor-p Arvore GPCM AGRUP COV` < `Nível de signi gpcm`) %>%
  arrange(`Valor-p Arvore GPCM AGRUP COV`)




base2<-cbind('GPCM categorizado com cov AGRUP',nrow(gpcm_sign_com_cov_repli_escore_categ))
colnames(base2)<- c('Metodo agrup','N_genes_sign')





genes_sign_gpcm_com_cov_categ_ISA_agrup_final <- gpcm_sign_com_cov_repli_escore_categ
View(genes_sign_gpcm_com_cov_categ_ISA_agrup_final)

num_genes_sign_tri_dados_ISA_agrup<- rbind(base,base2)
View(num_genes_sign_tri_dados_ISA_agrup)

##### EXPORTANDO
options(scipen = 999)

bases <- c(
  "num_genes_sign_tri_dados_ISA_agrup",
  "genes_sign_grm_com_cov_categ_ISA_agrup_final",
  "genes_sign_gpcm_com_cov_categ_ISA_agrup_final"
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



################# AQUI É A PARTE DE PRÉ-SELEÇÃO COM METODO TRADICIONAL COM ALPHA=5%
################ SEGUIDO DE APLICAÇÃO DA TRI

library(dplyr)
library(readxl)
library(readr)




############## METODOS TRADICIONAIS COM COV



num_genes_sign_burden_c_cov_005<- NULL

num_genes_sign_skat_c_cov_005<- NULL

num_genes_sign_skat_o_c_cov_005<- NULL


############ BANCO DE DADOS BURDEN
burden_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_burden_dados_ISA_005.txt")


#View(burden_com_cov_sign_repli)



burden_com_cov_sign_repli <- burden_com_cov_sign_repli %>%
  arrange(`Valor-p Burden com Cov`)

base<-cbind('Burden 5%',nrow(burden_com_cov_sign_repli))
colnames(base)<- c('Metodo pre-selecao','N_genes_sign')
View(base)

base_final_ident_burden_c_cov_005<-burden_com_cov_sign_repli

num_genes_sign_burden_c_cov_005<- rbind(num_genes_sign_burden_c_cov_005,base)
  
  
  
############# SKAT

skat_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_dados_ISA_005.txt")

#View(skat_com_cov_sign_repli)



skat_com_cov_sign_repli <- skat_com_cov_sign_repli %>%
  arrange(`Valor-p SKAT com Cov`)

base2<-cbind('SKAT 5%',nrow(skat_com_cov_sign_repli))
colnames(base2)<- c('Metodo pre-selecao','N_genes_sign')

base_final_ident_skat_c_cov_005<-skat_com_cov_sign_repli

num_genes_sign_skat_c_cov_005<- rbind(num_genes_sign_skat_c_cov_005,base2)



############# SKAT-O
skat_o_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_o_dados_ISA_005.txt")

#View(skat_o_com_cov_sign_repli)



skat_o_com_cov_sign_repli <- skat_o_com_cov_sign_repli %>%
  arrange(`Valor-p SKAT-O com Cov`)

base3<-cbind('SKAT-O 5%',nrow(skat_o_com_cov_sign_repli))
colnames(base3)<- c('Metodo pre-selecao','N_genes_sign')
  
  
base_final_ident_skat_o_c_cov_005<- skat_o_com_cov_sign_repli

num_genes_sign_skat_o_c_cov_005<- rbind(num_genes_sign_skat_o_c_cov_005,base3)



num_genes_sign_met_trad_cov_005 <- rbind(num_genes_sign_burden_c_cov_005,
                                         num_genes_sign_skat_c_cov_005,
                                         num_genes_sign_skat_o_c_cov_005)



#View(num_genes_sign_met_trad_cov_005)

#View(base_final_ident_skat_o_c_cov_005)


#####################################
options(scipen = 999)

bases <- c(
  "num_genes_sign_met_trad_cov_005"
)

#for (nome in bases) {

#  df <- get(nome)

#  write.table(
#    df,
#    file = paste0(nome, ".csv"),
#    sep = ";",
#    dec = ".",
#    row.names = FALSE,
#    col.names = TRUE
#  )
#}

######### VERIFICAR SE TÁ TUDO CERTO



######################## PUXA OS RESULTADOS DA TRI CATEG COM COV

############ GRM

library(readxl)
genes_identif_GRM_categ_c_cov <- read_delim("Mestrado/Defesa/resultados_dados_ISA/result_tri_dados_ISA/genes_sign_grm_com_cov_categ_ISA_final.csv", 
                                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

#View(genes_identif_GRM_categ_c_cov)


#View(base_final_ident_burden_c_cov_005)

### PEGA OS GENES PRÉ-SELECIONADOS VIA TRI E VERIFICA COMO ELES FICAM COM
## METODOLOGIAS TRADICIONAIS COM ALPHA 5% (PARA REDUZIR OS FALSO POSITIVOS IDENTIFICADOS NA TRI
## COM ALPHA 5% CORRIGIDO POR BONFERRONI)

## PRÉ SELEÇÃO COM A TRI
pre_sel_grm_categ_burden <- left_join(genes_identif_GRM_categ_c_cov,
                                      base_final_ident_burden_c_cov_005,
                                      by = c("Genes","Cromossomo"))

View(pre_sel_grm_categ_burden)

pre_sel_grm_categ_skat<- left_join(pre_sel_grm_categ_burden,
                                   base_final_ident_skat_c_cov_005,
                                   by = c("Genes","Cromossomo"))
#View(base_final_ident_skat_o_c_cov_005)

pre_sel_grm_categ_final<- left_join(pre_sel_grm_categ_skat,
                                    base_final_ident_skat_o_c_cov_005[,-c(3)],
                                    by = c("Genes","Cromossomo"))


pre_sel_grm_categ_final <- pre_sel_grm_categ_final %>%
  arrange(Cromossomo)
View(pre_sel_grm_categ_final)

### CORRIGE O ALPHA
base <- pre_sel_grm_categ_final %>%
  add_count(Cromossomo, name = "qtd_comb_cr_repli")
#View(base)
base$alpha_bf <-  0.05/base$qtd_comb_cr_repli

## BASE FINAL SIGNIF BURDEN PRÉ SELECIONADO POR GRM CATEG COM COV
burden_signif_GRM_categ_ISA <- base %>%
  filter(`Valor-p Burden com Cov` < alpha_bf) %>%
  arrange(`Valor-p Burden com Cov`) %>%
  select(
    Genes,
    Cromossomo,
    `Valor-p Burden com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

#View(burden_signif_GRM_categ_ISA)

## BASE FINAL SIGNIF SKAT PRÉ SELECIONADO POR GRM CATEG COM COV
skat_signif_GRM_categ_ISA <- base %>%
  filter(`Valor-p SKAT com Cov` < alpha_bf) %>%
  arrange(`Valor-p SKAT com Cov`) %>%
  select(
    Genes,
    Cromossomo,
    `Valor-p SKAT com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

#View(skat_signif_GRM_categ_ISA)

## BASE FINAL SIGNIF SKAT-O PRÉ SELECIONADO POR GRM CATEG COM COV
skat_o_signif_GRM_categ_ISA <- base %>%
  filter(`Valor-p SKAT-O com Cov` < alpha_bf) %>%
  arrange(`Valor-p SKAT-O com Cov`) %>%
  select(
    Genes,
    Cromossomo,
    `Valor-p SKAT-O com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

#View(skat_o_signif_GRM_categ_ISA)

#### ORDENANDO PELO VALOR-P



#### EXPORTA GENES SELECIONADOS PELOS MÉT TRADICIONAIS APÓS PRÉ SELEÇÃO PELA TRI CATEG
options(scipen = 999)

bases <- c(
  "burden_signif_GRM_categ_ISA",
  "skat_signif_GRM_categ_ISA",
  "skat_o_signif_GRM_categ_ISA"
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






################ GPCM


library(readxl)
genes_identif_GPCM_categ_c_cov <- read_delim("Mestrado/Defesa/resultados_dados_ISA/result_tri_dados_ISA/genes_sign_gpcm_com_cov_categ_ISA_final.csv", 
                                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(genes_identif_GPCM_categ_c_cov)


#View(base_final_ident_burden_c_cov_005)

## PRÉ SELEÇÃO COM A TRI
pre_sel_gpcm_categ_burden <- left_join(genes_identif_GPCM_categ_c_cov,
                                       base_final_ident_burden_c_cov_005,
                                       by = c("Genes","Cromossomo"))

pre_sel_gpcm_categ_skat<- left_join(pre_sel_gpcm_categ_burden,
                                    base_final_ident_skat_c_cov_005,
                                    by = c("Genes","Cromossomo"))
#View(base_final_ident_skat_o_c_cov_005)

pre_sel_gpcm_categ_final<- left_join(pre_sel_gpcm_categ_skat,
                                     base_final_ident_skat_o_c_cov_005[,-c(3)],
                                     by = c("Genes","Cromossomo"))


pre_sel_gpcm_categ_final <- pre_sel_gpcm_categ_final %>%
  arrange(Cromossomo)
#View(pre_sel_gpcm_categ_final)

### CORRIGE O ALPHA
base <- pre_sel_gpcm_categ_final %>%
  add_count(Cromossomo, name = "qtd_comb_cr_repli")
base$alpha_bf <-  0.05/base$qtd_comb_cr_repli

#View(base)

## BASE FINAL SIGNIF BURDEN PRÉ SELECIONADO POR GPCM CATEG COM COV
burden_signif_GPCM_categ_ISA <- base %>%
  filter(`Valor-p Burden com Cov` < alpha_bf) %>%
  arrange(`Valor-p Burden com Cov`) %>%
  select(
    Genes,
    Cromossomo,
    `Valor-p Burden com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

## BASE FINAL SIGNIF SKAT PRÉ SELECIONADO POR GPCM CATEG COM COV
skat_signif_GPCM_categ_ISA <- base %>%
  filter(`Valor-p SKAT com Cov` < alpha_bf) %>%
  arrange(`Valor-p SKAT com Cov`) %>%
  select(
    Genes,
    Cromossomo,
    `Valor-p SKAT com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

## BASE FINAL SIGNIF SKAT-O PRÉ SELECIONADO POR GPCM CATEG COM COV
skat_o_signif_GPCM_categ_ISA <- base %>%
  filter(`Valor-p SKAT-O com Cov` < alpha_bf) %>%
  arrange(`Valor-p SKAT-O com Cov`) %>%
  select(
    Genes,
    Cromossomo,
    `Valor-p SKAT-O com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

#View(skat_signif_GPCM_categ)

#### EXPORTA GENES SELECIONADOS PELOS MÉT TRADICIONAIS APÓS PRÉ SELEÇÃO PELA TRI CATEG
options(scipen = 999)

bases <- c(
  "burden_signif_GPCM_categ_ISA",
  "skat_signif_GPCM_categ_ISA",
  "skat_o_signif_GPCM_categ_ISA"
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

n1<- cbind("Burden - GRM",nrow(burden_signif_GRM_categ_ISA))
n2<- cbind("SKAT - GRM",nrow(skat_signif_GRM_categ_ISA))
n3<- cbind("SKAT-O - GRM",nrow(skat_o_signif_GRM_categ_ISA))
num_genes_GRM <- rbind(n1,n2,n3)




n1<- cbind("Burden - GPCM",nrow(burden_signif_GPCM_categ_ISA))
n2<- cbind("SKAT - GPCM",nrow(skat_signif_GPCM_categ_ISA))
n3<- cbind("SKAT-O - GPCM",nrow(skat_o_signif_GPCM_categ_ISA))
num_genes_GPCM <- rbind(n1,n2,n3)

num_genes_pre_sel_final <- rbind(num_genes_GRM,num_genes_GPCM)
colnames(num_genes_pre_sel_final)<- c('Metodo combinado','N_genes_sign')


bases <- c(
  "num_genes_pre_sel_final"
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

############ VENDO QUAIS FORAM SELECIONADOS POR TODOS OS MÉTODOS

View(burden_signif_GPCM_categ_ISA)
View(skat_signif_GPCM_categ_ISA)
View(skat_o_signif_GPCM_categ_ISA)
View(genes_sign_grm_com_cov_categ_ISA_final)
View(genes_sign_gpcm_com_cov_categ_ISA_final)

View(burden_signif_GRM_categ_ISA)
View(skat_signif_GRM_categ_ISA)
View(skat_o_signif_GRM_categ_ISA)



met_1 <- full_join(
  genes_sign_grm_com_cov_categ_ISA_final %>%
    select(Genes, Cromossomo, `Valor-p Arvore GRM C COV`),
  
  genes_sign_gpcm_com_cov_categ_ISA_final %>%
    select(Genes, Cromossomo, `Valor-p Arvore GPCM C COV`),
  
  by = c("Genes", "Cromossomo")
)

met_2<- full_join(
  met_1,
  skat_o_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, 
           `Valor-p SKAT-O GPCM` = `Valor-p SKAT-O com Cov`),
  
  by = c("Genes", "Cromossomo")
)


met_3<- full_join(
  met_2 ,skat_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, 
           `Valor-p SKAT GPCM` = `Valor-p SKAT com Cov`),
  by = c("Genes", "Cromossomo")
)
  


met_4<- full_join(met_3,burden_signif_GPCM_categ_ISA%>%
                  select(Genes, Cromossomo, 
                           `Valor-p Burden GPCM` = `Valor-p Burden com Cov`),
                  by = c("Genes","Cromossomo"))

View(met_4)


met_5<- full_join(met_4,skat_o_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p SKAT-O GRM` = `Valor-p SKAT-O com Cov`),
                  by = c("Genes","Cromossomo"))

met_6<- full_join(met_5,skat_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p SKAT GRM` = `Valor-p SKAT com Cov`),
                  by = c("Genes","Cromossomo"))


met_7<- full_join(met_6,burden_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p Burden GRM` = `Valor-p Burden com Cov`),
                  by = c("Genes","Cromossomo"))

View(met_7)



skat_o_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_o_dados_ISA.txt")



met_final<- full_join(met_7,skat_o_com_cov_sign_repli%>%
                        select(Genes, Cromossomo,`Valor-p SKAT-O com Cov`),
                              by = c("Genes","Cromossomo"))

View(met_final)


compara_result_trad_tri_pre_sel_full_join<- met_final


bases <- c(
  "compara_result_trad_tri_pre_sel_full_join"
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



######## INNER




tri <- inner_join(
  genes_sign_grm_com_cov_categ_ISA_final %>%
    select(Genes, Cromossomo, `Valor-p Arvore GRM C COV`),
  
  genes_sign_gpcm_com_cov_categ_ISA_final %>%
    select(Genes, Cromossomo, `Valor-p Arvore GPCM C COV`),
  
  by = c("Genes", "Cromossomo")
)

burden_pre_sel <- inner_join(
  burden_signif_GRM_categ_ISA %>%
    select(Genes, Cromossomo, `Valor-p Burden GRM` = `Valor-p Burden com Cov`),
  
  burden_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, `Valor-p Burden GPCM` = `Valor-p Burden com Cov`),
  
  by = c("Genes", "Cromossomo")
)

skat_pre_sel <- inner_join(
  skat_signif_GRM_categ_ISA %>%
    select(Genes, Cromossomo, `Valor-p SKAT GRM` = `Valor-p SKAT com Cov`),
  
  skat_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, `Valor-p SKAT GPCM` = `Valor-p SKAT com Cov`),
  
  by = c("Genes", "Cromossomo")
)

skat_o_pre_sel<- inner_join(
  skat_o_signif_GRM_categ_ISA %>%
    select(Genes, Cromossomo, `Valor-p SKAT-O GRM` = `Valor-p SKAT-O com Cov`),
  
  skat_o_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, `Valor-p SKAT-O GPCM` = `Valor-p SKAT-O com Cov`),
  
  by = c("Genes", "Cromossomo")
)
  
inner_1<- inner_join(tri,
                     skat_o_pre_sel,
                    by = c("Genes","Cromossomo"))

inner_2<- inner_join(inner_1,
                     skat_pre_sel,
                     by = c("Genes","Cromossomo"))



inner_final <- inner_join(inner_2,
                          burden_pre_sel,
                          by = c("Genes", "Cromossomo")
                          ) %>%
                arrange(`Valor-p Arvore GRM C COV`)

View(inner_final)



compara_result_trad_tri_pre_sel_inner_join<- inner_final



bases <- c(
  "compara_result_trad_tri_pre_sel_inner_join"
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


##### TOP 10 GENES MAIS SIGNIF GPCM E LEFT COM AS OUTRAS METODOLOGIAS

met_1 <- left_join(
  genes_sign_gpcm_com_cov_categ_ISA_final %>%
    arrange(`Valor-p Arvore GPCM C COV`) %>%
    head(10) %>%
    select(Genes, Cromossomo, `Valor-p Arvore GPCM C COV`),
  
  genes_sign_grm_com_cov_categ_ISA_final %>%
    select(Genes, Cromossomo, `Valor-p Arvore GRM C COV`),
  
  by = c("Genes", "Cromossomo")
)

### GRM PRE SEL

met_2<- left_join(met_1,burden_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p Burden GRM` = `Valor-p Burden com Cov`),
                  by = c("Genes","Cromossomo"))

met_3<- left_join(met_2,skat_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p SKAT GRM` = `Valor-p SKAT com Cov`),
                  by = c("Genes","Cromossomo"))

met_4<- left_join(met_3,skat_o_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p SKAT-O GRM` = `Valor-p SKAT-O com Cov`),
                  by = c("Genes","Cromossomo"))



### GPCM PRE SEL

met_5<- left_join(met_4,burden_signif_GPCM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p Burden GPCM` = `Valor-p Burden com Cov`),
                  by = c("Genes","Cromossomo"))


met_6<- left_join(
  met_5 ,skat_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, 
           `Valor-p SKAT GPCM` = `Valor-p SKAT com Cov`),
  by = c("Genes", "Cromossomo")
)

met_7<- left_join(
  met_6,
  skat_o_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, 
           `Valor-p SKAT-O GPCM` = `Valor-p SKAT-O com Cov`),
  
  by = c("Genes", "Cromossomo")
)




skat_o_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_o_dados_ISA.txt")



met_final<- left_join(met_7,skat_o_com_cov_sign_repli%>%
                        select(Genes, Cromossomo,`Valor-p SKAT-O com Cov`),
                      by = c("Genes","Cromossomo"))

View(met_final)


top10_result_trad_tri_pre_sel_dados_ISA<- met_final


bases <- c(
  "top10_result_trad_tri_pre_sel_dados_ISA"
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




##### TOP 10 GENES MAIS SIGNIF GPCM E FULL JOIN COM AS OUTRAS METODOLOGIAS

top10 <- left_join(
  genes_sign_gpcm_com_cov_categ_ISA_final %>%
    arrange(`Valor-p Arvore GPCM C COV`) %>%
    head(10) %>%
    select(Genes, Cromossomo, `Valor-p Arvore GPCM C COV`),
  
  genes_sign_grm_com_cov_categ_ISA_final %>%
    select(Genes, Cromossomo, `Valor-p Arvore GRM C COV`),
  
  by = c("Genes", "Cromossomo")
)

### GRM PRE SEL

met_2<- full_join(met_1,burden_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p Burden GRM` = `Valor-p Burden com Cov`),
                  by = c("Genes","Cromossomo"))

met_3<- full_join(met_2,skat_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p SKAT GRM` = `Valor-p SKAT com Cov`),
                  by = c("Genes","Cromossomo"))

met_4<- full_join(met_3,skat_o_signif_GRM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p SKAT-O GRM` = `Valor-p SKAT-O com Cov`),
                  by = c("Genes","Cromossomo"))



### GPCM PRE SEL

met_5<- full_join(met_4,burden_signif_GPCM_categ_ISA%>%
                    select(Genes, Cromossomo, 
                           `Valor-p Burden GPCM` = `Valor-p Burden com Cov`),
                  by = c("Genes","Cromossomo"))


met_6<- full_join(
  met_5 ,skat_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, 
           `Valor-p SKAT GPCM` = `Valor-p SKAT com Cov`),
  by = c("Genes", "Cromossomo")
)

met_7<- full_join(
  met_6,
  skat_o_signif_GPCM_categ_ISA %>%
    select(Genes, Cromossomo, 
           `Valor-p SKAT-O GPCM` = `Valor-p SKAT-O com Cov`),
  
  by = c("Genes", "Cromossomo")
)




skat_o_com_cov_sign_repli <- read_csv("Mestrado/Defesa/resultados_dados_ISA/result_met_trad_dados_ISA/resul_modelo_c_cov_skat_o_dados_ISA.txt")



met_antepenul<- full_join(met_7,skat_o_com_cov_sign_repli%>%
                        select(Genes, Cromossomo,`Valor-p SKAT-O com Cov`),
                      by = c("Genes","Cromossomo"))


met_penul<- left_join(met_antepenul,genes_sign_gpcm_com_cov_categ_ISA_final%>%
                            select(Genes, Cromossomo,`Valor-p Arvore GPCM C COV`),
                          by = c("Genes","Cromossomo"))

met_ultl<- left_join(met_penul,genes_sign_grm_com_cov_categ_ISA_final%>%
                        select(Genes, Cromossomo,`Valor-p Arvore GRM C COV`),
                      by = c("Genes","Cromossomo"))

View(met_ultl)


junta_result_trad_tri_pre_sel_dados_ISA_certo<- met_ultl


bases <- c(
  "junta_result_trad_tri_pre_sel_dados_ISA_certo"
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




