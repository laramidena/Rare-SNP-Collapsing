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




############## METODOS TRADICIONAIS COM COV


base_final_sens_burden_c_cov_005<- NULL
base_final_esp_burden_c_cov_005<- NULL
base_final_f1_burden_c_cov_005<- NULL
f1_burden_c_cov<- NULL

base_final_sens_skat_c_cov_005<- NULL
base_final_esp_skat_c_cov_005<- NULL
base_final_f1_skat_c_cov_005<- NULL
f1_skat_c_cov<- NULL

base_final_sens_skat_o_c_cov_005<- NULL
base_final_esp_skat_o_c_cov_005<- NULL
base_final_f1_skat_o_c_cov_005<- NULL
f1_skat_o_c_cov<- NULL

num_genes_sign_burden_c_cov_005<- NULL

num_genes_sign_skat_c_cov_005<- NULL

num_genes_sign_skat_o_c_cov_005<- NULL

base_final_ident_burden_c_cov_005<- NULL
base_final_ident_skat_c_cov_005<- NULL
base_final_ident_skat_o_c_cov_005 <- NULL


for (repli in (c(1:10))){
  ############ BANCO DE DADOS BURDEN REPLICA
  burden_com_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/20260714_result_met_trad_c_cov_alpha_005/resul_modelo_c_cov_burden_repli",repli,".txt"))
  
  
  #View(burden_com_cov_sign_repli)
  
  
  
  burden_com_cov_sign_repli <- burden_com_cov_sign_repli %>%
    arrange(`Valor-p Burden com Cov`)
  
  base<-cbind(repli[1],nrow(burden_com_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_burden_c_cov_005<- rbind(num_genes_sign_burden_c_cov_005,base)
  
  
  especif_burden_com_cov_repli <- inner_join(burden_com_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_burden_com_cov_repli$E<- c(rep(((3196-nrow(especif_burden_com_cov_repli))/(3196)),
                                         nrow(especif_burden_com_cov_repli)))
  
  #View(especif_burden_com_cov_repli)
  base_final_esp_burden_c_cov_005<- rbind(base_final_esp_burden_c_cov_005,especif_burden_com_cov_repli)
  #View(base_final_esp_burden_c_cov_005)
  
  
  ### F1-SCORE
  
  fp_burden_c_cov<- nrow(inner_join(burden_com_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_burden_c_cov<- nrow(inner_join(burden_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_burden_c_cov<- nrow(anti_join(simula_gaw17, burden_com_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_burden_c_cov$f1<- (2*tp_burden_c_cov)/((2*tp_burden_c_cov) + fn_burden_c_cov + fp_burden_c_cov)
  
  
  base_final_f1_burden_c_cov_005<- rbind(base_final_f1_burden_c_cov_005,f1_burden_c_cov)
  
  
  sensib_burden_com_cov_repli <- inner_join(burden_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_burden_com_cov_repli$S<- c(rep(nrow(sensib_burden_com_cov_repli)/nrow(simula_gaw17),nrow(sensib_burden_com_cov_repli)))
  #View(sensib_burden_com_cov_repli)
  
  
  base_final_sens_burden_c_cov_005<- rbind(base_final_sens_burden_c_cov_005,sensib_burden_com_cov_repli)
  #View(base_final_sens_burden_c_cov_005)
  
  base_final_ident_burden_c_cov_005<- rbind(base_final_sens_burden_c_cov_005[,-ncol(base_final_sens_burden_c_cov_005)],
                                            base_final_esp_burden_c_cov_005[,-ncol(base_final_esp_burden_c_cov_005)])
  
  ############# SKAT
  
  skat_com_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/20260714_result_met_trad_c_cov_alpha_005/resul_modelo_c_cov_skat_repli",repli,".txt"))
  
  #View(skat_com_cov_sign_repli)
  
  
  
  skat_com_cov_sign_repli <- skat_com_cov_sign_repli %>%
    arrange(`Valor-p SKAT com Cov`)
  
  base<-cbind(repli[1],nrow(skat_com_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_skat_c_cov_005<- rbind(num_genes_sign_skat_c_cov_005,base)
  
  
  especif_skat_com_cov_repli <- inner_join(skat_com_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_skat_com_cov_repli$E<- c(rep(((3196-nrow(especif_skat_com_cov_repli))/(3196)),
                                       nrow(especif_skat_com_cov_repli)))
  
  #View(especif_skat_com_cov_repli)
  base_final_esp_skat_c_cov_005<- rbind(base_final_esp_skat_c_cov_005,especif_skat_com_cov_repli)
  #View(base_final_esp_skat_c_cov_005)
  
  
  ### F1-SCORE
  
  fp_skat_c_cov<- nrow(inner_join(skat_com_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_skat_c_cov<- nrow(inner_join(skat_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_skat_c_cov<- nrow(anti_join(simula_gaw17, skat_com_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_skat_c_cov$f1<- (2*tp_skat_c_cov)/((2*tp_skat_c_cov) + fn_skat_c_cov + fp_skat_c_cov)
  
  
  base_final_f1_skat_c_cov_005<- rbind(base_final_f1_skat_c_cov_005,f1_skat_c_cov)
  
  
  
  sensib_skat_com_cov_repli <- inner_join(skat_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_skat_com_cov_repli$S<- c(rep(nrow(sensib_skat_com_cov_repli)/nrow(simula_gaw17),nrow(sensib_skat_com_cov_repli)))
  #View(sensib_skat_com_cov_repli)
  
  
  base_final_sens_skat_c_cov_005<- rbind(base_final_sens_skat_c_cov_005,sensib_skat_com_cov_repli)
  #View(base_final_sens_skat_c_cov_005)
  
  base_final_ident_skat_c_cov_005<- rbind(base_final_sens_skat_c_cov_005[,-ncol(base_final_sens_skat_c_cov_005)],
                                          base_final_esp_skat_c_cov_005[,-ncol(base_final_esp_skat_c_cov_005)])
  
  ############# SKAT-O
  skat_o_com_cov_sign_repli <- read_csv(paste0("C:/Users/laram/OneDrive/Documentos/Mestrado/Defesa/resultados/resultados_replicas_full/result_met_trad/20260714_result_met_trad_c_cov_alpha_005/resul_modelo_c_cov_skat_o_repli",repli,".txt"))
  
  #View(skat_o_com_cov_sign_repli)
  
  
  
  skat_o_com_cov_sign_repli <- skat_o_com_cov_sign_repli %>%
    arrange(`Valor-p SKAT-O com Cov`)
  
  base<-cbind(repli[1],nrow(skat_o_com_cov_sign_repli))
  colnames(base)<- c('Replica','N_genes_sign')
  
  
  num_genes_sign_skat_o_c_cov_005<- rbind(num_genes_sign_skat_o_c_cov_005,base)
  
  
  especif_skat_o_com_cov_repli <- inner_join(skat_o_com_cov_sign_repli, genes_n_sign, by = c("Genes"))
  especif_skat_o_com_cov_repli$E<- c(rep(((3196-nrow(especif_skat_o_com_cov_repli))/(3196)),
                                         nrow(especif_skat_o_com_cov_repli)))
  
  #View(especif_skat_o_com_cov_repli)
  base_final_esp_skat_o_c_cov_005<- rbind(base_final_esp_skat_o_c_cov_005,especif_skat_o_com_cov_repli)
  #View(base_final_esp_skat_o_c_cov_005)
  
  ### F1-SCORE
  
  fp_skat_o_c_cov<- nrow(inner_join(skat_o_com_cov_sign_repli, genes_n_sign, by = c("Genes")))
  tp_skat_o_c_cov<- nrow(inner_join(skat_o_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo")))
  fn_skat_o_c_cov<- nrow(anti_join(simula_gaw17, skat_o_com_cov_sign_repli, by = c("Genes","Cromossomo")))
  
  f1_skat_o_c_cov$f1<- (2*tp_skat_o_c_cov)/((2*tp_skat_o_c_cov) + fn_skat_o_c_cov + fp_skat_o_c_cov)
  
  
  base_final_f1_skat_o_c_cov_005<- rbind(base_final_f1_skat_o_c_cov_005,f1_skat_o_c_cov)
  
  
  
  sensib_skat_o_com_cov_repli <- inner_join(skat_o_com_cov_sign_repli, simula_gaw17, by = c("Genes","Cromossomo"))
  sensib_skat_o_com_cov_repli$S<- c(rep(nrow(sensib_skat_o_com_cov_repli)/nrow(simula_gaw17),nrow(sensib_skat_o_com_cov_repli)))
  #View(sensib_skat_o_com_cov_repli)
  
  
  base_final_sens_skat_o_c_cov_005<- rbind(base_final_sens_skat_o_c_cov_005,sensib_skat_o_com_cov_repli)
  #View(base_final_sens_skat_o_c_cov_005)
  
  base_final_ident_skat_o_c_cov_005<- rbind(base_final_sens_skat_o_c_cov_005[,-ncol(base_final_sens_skat_o_c_cov_005)],
                                            base_final_esp_skat_o_c_cov_005[,-ncol(base_final_esp_skat_o_c_cov_005)])
}



#View(base_final_sens_burden_c_cov_005)
#View(base_final_esp_burden_c_cov_005)
#View(base_final_f1_burden_c_cov_005)

#View(base_final_sens_skat_c_cov_005)
#View(base_final_esp_skat_c_cov_005)
#View(base_final_f1_skat_c_cov_005)

#View(base_final_sens_skat_o_c_cov_005)
#View(base_final_esp_skat_o_c_cov_005)
#View(base_final_f1_skat_o_c_cov_005)

#View(num_genes_sign_burden_c_cov_005)

#View(num_genes_sign_skat_c_cov_005)

#View(num_genes_sign_skat_o_c_cov_005)

View(base_final_ident_burden_c_cov_005)


#####################################
options(scipen = 999)

bases <- c(
  "base_final_sens_burden_c_cov_005",
  "base_final_sens_skat_c_cov_005",
  "base_final_sens_skat_o_c_cov_005",
  "base_final_esp_burden_c_cov_005",
  "base_final_esp_skat_c_cov_005",
  "base_final_esp_skat_o_c_cov_005",
  "num_genes_sign_burden_c_cov_005",
  "num_genes_sign_skat_c_cov_005",
  "num_genes_sign_skat_o_c_cov_005",
  "base_final_f1_burden_c_cov_005",
  "base_final_f1_skat_c_cov_005",
  "base_final_f1_skat_o_c_cov_005",
  "base_final_ident_burden_c_cov_005",
  "base_final_ident_skat_c_cov_005",
  "base_final_ident_skat_o_c_cov_005"
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


######################## PUXA OS RESULTADOS DA TRI CATEG COM COV

############ GRM

library(readxl)
genes_identif_GRM_categ_c_cov <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260419_result_replicas_c_cov/GRM_categ/genes_identif_GRM_categ_c_cov.xlsx")
View(genes_identif_GRM_categ_c_cov)


#View(base_final_ident_burden_c_cov_005)

## PRÉ SELEÇÃO COM A TRI
pre_sel_grm_categ_burden <- left_join(genes_identif_GRM_categ_c_cov,
                                      base_final_ident_burden_c_cov_005,
                                      by = c("Genes","Cromossomo","replica" = "Replica"))

pre_sel_grm_categ_skat<- left_join(pre_sel_grm_categ_burden,
                                     base_final_ident_skat_c_cov_005,
                                     by = c("Genes","Cromossomo","replica" = "Replica"))
View(base_final_ident_skat_o_c_cov_005)

pre_sel_grm_categ_final<- left_join(pre_sel_grm_categ_skat,
                                     base_final_ident_skat_o_c_cov_005[,-c(3)],
                                     by = c("Genes","Cromossomo","replica" = "Replica"))


pre_sel_grm_categ_final <- pre_sel_grm_categ_final %>%
  arrange(replica,Cromossomo)
View(pre_sel_grm_categ_final)

### CORRIGE O ALPHA
base <- pre_sel_grm_categ_final %>%
  add_count(Cromossomo, replica, name = "qtd_comb_cr_repli")
base$alpha_bf <-  0.05/base$qtd_comb_cr_repli

## BASE FINAL SIGNIF BURDEN PRÉ SELECIONADO POR GRM CATEG COM COV
burden_signif_GRM_categ <- base %>%
  filter(`Valor-p Burden com Cov` < alpha_bf) %>%
  select(
    Genes,
    Cromossomo,
    replica,
    `Valor-p Burden com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

## BASE FINAL SIGNIF SKAT PRÉ SELECIONADO POR GRM CATEG COM COV
skat_signif_GRM_categ <- base %>%
  filter(`Valor-p SKAT com Cov` < alpha_bf) %>%
  select(
    Genes,
    Cromossomo,
    replica,
    `Valor-p SKAT com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

## BASE FINAL SIGNIF SKAT-O PRÉ SELECIONADO POR GRM CATEG COM COV
skat_o_signif_GRM_categ <- base %>%
  filter(`Valor-p SKAT-O com Cov` < alpha_bf) %>%
  select(
    Genes,
    Cromossomo,
    replica,
    `Valor-p SKAT-O com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

View(skat_signif_GRM_categ)

#### EXPORTA GENES SELECIONADOS PELOS MÉT TRADICIONAIS APÓS PRÉ SELEÇÃO PELA TRI CATEG
options(scipen = 999)

bases <- c(
  "burden_signif_GRM_categ",
  "skat_signif_GRM_categ",
  "skat_o_signif_GRM_categ"
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



#### VDD POSITIVO
simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
#View(simula_gaw17)



snp_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/snp_info")
#View(snp_info)
snp_info<-snp_info[,4]
genes <- unique(snp_info)
colnames(genes)<- "Genes"
#View(genes)


genes_n_sign<- anti_join(genes, simula_gaw17, by = c("Genes"))
View(genes_n_sign)




#### FUNÇÃO SENSIBILIDADE, ESPECIFICIDADE, F1-SCORE E NUM DE GENES GRM

library(dplyr)

calcula_metricas <- function(resultados,
                             simula_gaw17,
                             genes_n_sign,
                             n_negativos = 3196){
  
  n_vp <- n_distinct(simula_gaw17$Genes)
  
  ## Número de genes identificados
  num_genes <- resultados %>%
    group_by(replica) %>%
    summarise(
      Genes = n(),
      .groups = "drop"
    )
  
  ## Especificidade
  especificidade <- resultados %>%
    inner_join(genes_n_sign, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      FP = n_distinct(Genes),
      Especificidade = (n_negativos - FP) / n_negativos,
      .groups = "drop"
    )
  
  ## FP para F1
  fp <- resultados %>%
    inner_join(genes_n_sign, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      FP = n(),
      .groups = "drop"
    )
  
  ## TP para F1
  tp <- resultados %>%
    inner_join(simula_gaw17, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      TP = n(),
      .groups = "drop"
    )
  
  ## FN para F1
  fn <- resultados %>%
    group_by(replica) %>%
    summarise(
      FN = nrow(
        anti_join(
          simula_gaw17,
          cur_data(),
          by = c("Genes", "Cromossomo")
        )
      ),
      .groups = "drop"
    )
  
  ## F1-score
  f1 <- tp %>%
    left_join(fp, by = "replica") %>%
    left_join(fn, by = "replica") %>%
    mutate(
      F1 = (2 * TP) / (2 * TP + FP + FN)
    ) %>%
    select(replica, F1)
  
  ## Sensibilidade
  sensibilidade <- resultados %>%
    inner_join(simula_gaw17, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      TP = n_distinct(Genes),
      Sensibilidade = TP / n_vp,
      .groups = "drop"
    )
  
  ## Junta tudo
  num_genes %>%
    left_join(sensibilidade %>% select(replica, Sensibilidade), by = "replica") %>%
    left_join(especificidade %>% select(replica, Especificidade), by = "replica") %>%
    left_join(f1, by = "replica")
}

#-------------------------------------------------------
# Aplica para os três métodos
#-------------------------------------------------------

metricas_burden <- calcula_metricas(
  burden_signif_GRM_categ,
  simula_gaw17,
  genes_n_sign
)


metricas_skat <- calcula_metricas(
  skat_signif_GRM_categ,
  simula_gaw17,
  genes_n_sign
)

metricas_skato <- calcula_metricas(
  skat_o_signif_GRM_categ,
  simula_gaw17,
  genes_n_sign
)

#-------------------------------------------------------
# Bases finais
#-------------------------------------------------------

sensibilidade_GRM <- metricas_burden %>%
  select(replica, burden = Sensibilidade) %>%
  left_join(metricas_skat %>%
              select(replica, skat = Sensibilidade),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = Sensibilidade),
            by = "replica")

especificidade_GRM <- metricas_burden %>%
  select(replica, burden = Especificidade) %>%
  left_join(metricas_skat %>%
              select(replica, skat = Especificidade),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = Especificidade),
            by = "replica")

num_genes_GRM <- metricas_burden %>%
  select(replica, burden = Genes) %>%
  left_join(metricas_skat %>%
              select(replica, skat = Genes),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = Genes),
            by = "replica")

f1_score_GRM <- metricas_burden %>%
  select(replica, burden = F1) %>%
  left_join(metricas_skat %>%
              select(replica, skat = F1),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = F1),
            by = "replica")


View(sensibilidade_GRM)
View(especificidade_GRM)
View(num_genes_GRM)
View(f1_score_GRM)


options(scipen = 999)

bases <- c(
  "sensibilidade_GRM",
  "especificidade_GRM",
  "num_genes_GRM",
  "f1_score_GRM"
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
genes_identif_GPCM_categ_c_cov <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/20260419_result_replicas_c_cov/GPCM_categ/genes_identif_GPCM_categ_c_cov.xlsx")
View(genes_identif_GPCM_categ_c_cov)


#View(base_final_ident_burden_c_cov_005)

## PRÉ SELEÇÃO COM A TRI
pre_sel_gpcm_categ_burden <- left_join(genes_identif_GPCM_categ_c_cov,
                                      base_final_ident_burden_c_cov_005,
                                      by = c("Genes","Cromossomo","replica" = "Replica"))

pre_sel_gpcm_categ_skat<- left_join(pre_sel_gpcm_categ_burden,
                                   base_final_ident_skat_c_cov_005,
                                   by = c("Genes","Cromossomo","replica" = "Replica"))
#View(base_final_ident_skat_o_c_cov_005)

pre_sel_gpcm_categ_final<- left_join(pre_sel_gpcm_categ_skat,
                                    base_final_ident_skat_o_c_cov_005[,-c(3)],
                                    by = c("Genes","Cromossomo","replica" = "Replica"))


pre_sel_gpcm_categ_final <- pre_sel_gpcm_categ_final %>%
  arrange(replica,Cromossomo)
#View(pre_sel_gpcm_categ_final)

### CORRIGE O ALPHA
base <- pre_sel_gpcm_categ_final %>%
  add_count(Cromossomo, replica, name = "qtd_comb_cr_repli")
base$alpha_bf <-  0.05/base$qtd_comb_cr_repli

#View(base)

## BASE FINAL SIGNIF BURDEN PRÉ SELECIONADO POR GPCM CATEG COM COV
burden_signif_GPCM_categ <- base %>%
  filter(`Valor-p Burden com Cov` < alpha_bf) %>%
  select(
    Genes,
    Cromossomo,
    replica,
    `Valor-p Burden com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

## BASE FINAL SIGNIF SKAT PRÉ SELECIONADO POR GPCM CATEG COM COV
skat_signif_GPCM_categ <- base %>%
  filter(`Valor-p SKAT com Cov` < alpha_bf) %>%
  select(
    Genes,
    Cromossomo,
    replica,
    `Valor-p SKAT com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

## BASE FINAL SIGNIF SKAT-O PRÉ SELECIONADO POR GPCM CATEG COM COV
skat_o_signif_GPCM_categ <- base %>%
  filter(`Valor-p SKAT-O com Cov` < alpha_bf) %>%
  select(
    Genes,
    Cromossomo,
    replica,
    `Valor-p SKAT-O com Cov`,
    qtd_comb_cr_repli,
    alpha_bf
  )

#View(skat_signif_GPCM_categ)

#### EXPORTA GENES SELECIONADOS PELOS MÉT TRADICIONAIS APÓS PRÉ SELEÇÃO PELA TRI CATEG
options(scipen = 999)

bases <- c(
  "burden_signif_GPCM_categ",
  "skat_signif_GPCM_categ",
  "skat_o_signif_GPCM_categ"
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



#### VDD POSITIVO
simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
#View(simula_gaw17)



snp_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/snp_info")
#View(snp_info)
snp_info<-snp_info[,4]
genes <- unique(snp_info)
colnames(genes)<- "Genes"
#View(genes)


genes_n_sign<- anti_join(genes, simula_gaw17, by = c("Genes"))
#View(genes_n_sign)




#### FUNÇÃO SENSIBILIDADE, ESPECIFICIDADE, F1-SCORE E NUM DE GENES GPCM

library(dplyr)

calcula_metricas <- function(resultados,
                             simula_gaw17,
                             genes_n_sign,
                             n_negativos = 3196){
  
  n_vp <- n_distinct(simula_gaw17$Genes)
  
  ## Número de genes identificados
  num_genes <- resultados %>%
    group_by(replica) %>%
    summarise(
      Genes = n(),
      .groups = "drop"
    )
  
  ## Especificidade
  especificidade <- resultados %>%
    inner_join(genes_n_sign, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      FP = n_distinct(Genes),
      Especificidade = (n_negativos - FP) / n_negativos,
      .groups = "drop"
    )
  
  ## FP para F1
  fp <- resultados %>%
    inner_join(genes_n_sign, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      FP = n(),
      .groups = "drop"
    )
  
  ## TP para F1
  tp <- resultados %>%
    inner_join(simula_gaw17, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      TP = n(),
      .groups = "drop"
    )
  
  ## FN para F1
  fn <- resultados %>%
    group_by(replica) %>%
    summarise(
      FN = nrow(
        anti_join(
          simula_gaw17,
          cur_data(),
          by = c("Genes", "Cromossomo")
        )
      ),
      .groups = "drop"
    )
  
  ## F1-score
  f1 <- tp %>%
    left_join(fp, by = "replica") %>%
    left_join(fn, by = "replica") %>%
    mutate(
      F1 = (2 * TP) / (2 * TP + FP + FN)
    ) %>%
    select(replica, F1)
  
  ## Sensibilidade
  sensibilidade <- resultados %>%
    inner_join(simula_gaw17, by = "Genes") %>%
    group_by(replica) %>%
    summarise(
      TP = n_distinct(Genes),
      Sensibilidade = TP / n_vp,
      .groups = "drop"
    )
  
  ## Junta tudo
  num_genes %>%
    left_join(sensibilidade %>% select(replica, Sensibilidade), by = "replica") %>%
    left_join(especificidade %>% select(replica, Especificidade), by = "replica") %>%
    left_join(f1, by = "replica")
}

#-------------------------------------------------------
# Aplica para os três métodos
#-------------------------------------------------------

metricas_burden <- calcula_metricas(
  burden_signif_GPCM_categ,
  simula_gaw17,
  genes_n_sign
)


metricas_skat <- calcula_metricas(
  skat_signif_GPCM_categ,
  simula_gaw17,
  genes_n_sign
)

metricas_skato <- calcula_metricas(
  skat_o_signif_GPCM_categ,
  simula_gaw17,
  genes_n_sign
)

#-------------------------------------------------------
# Bases finais
#-------------------------------------------------------

sensibilidade_GPCM <- metricas_burden %>%
  select(replica, burden = Sensibilidade) %>%
  left_join(metricas_skat %>%
              select(replica, skat = Sensibilidade),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = Sensibilidade),
            by = "replica")

especificidade_GPCM <- metricas_burden %>%
  select(replica, burden = Especificidade) %>%
  left_join(metricas_skat %>%
              select(replica, skat = Especificidade),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = Especificidade),
            by = "replica")

num_genes_GPCM <- metricas_burden %>%
  select(replica, burden = Genes) %>%
  left_join(metricas_skat %>%
              select(replica, skat = Genes),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = Genes),
            by = "replica")

f1_score_GPCM <- metricas_burden %>%
  select(replica, burden = F1) %>%
  left_join(metricas_skat %>%
              select(replica, skat = F1),
            by = "replica") %>%
  left_join(metricas_skato %>%
              select(replica, skat_o = F1),
            by = "replica")




options(scipen = 999)

bases <- c(
  "sensibilidade_GPCM",
  "especificidade_GPCM",
  "num_genes_GPCM",
  "f1_score_GPCM"
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

############## ARRUMAR AQUI

### VERIFICANDO SE GENES IDENTIFICADOS ERRADOS ESTÃO PRÓXIMOS DOS CERTOS -  PARA PRÉ SELEÇÃO VIA METODOS TRI GRM 
metodos_tri<- c(
  "burden_GRM",
  "skat_GRM",
  "skat_o_GRM"
)

library(dplyr)


for (metodo in metodos_tri){
  simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
  #View(simula_gaw17)
  gene_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/gene_info.txt")
  
  colnames(gene_info)<- c("Genes","Cromossomo","Posicao","Gene_start","Gene_end","Gene_length","Orientation")
  
  
  ### gene falso positivo - metodos
  top_10_metodo <- read_excel(paste0("Mestrado/Defesa/resultados/resultados_replicas_full/20260728_result_replicas_2etapas/GRM/top10_",metodo,".xlsx"))
  
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
    file = paste0("verif_prox_genes_errados_pre_sel_",metodo, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}







### VERIFICANDO SE GENES IDENTIFICADOS ERRADOS ESTÃO PRÓXIMOS DOS CERTOS -  PARA PRÉ SELEÇÃO VIA METODOS TRI GPCM 
metodos_tri<- c(
  "burden_GPCM",
  "skat_GPCM",
  "skat_o_GPCM"
)



for (metodo in metodos_tri){
  simula_gaw17 <- read_excel("Mestrado/Defesa/resultados/resultados_replicas_full/simula_gaw17.xlsx")
  #View(simula_gaw17)
  gene_info <- read_csv("Mestrado/Defesa/resultados/resultados_replicas_full/gene_info.txt")
  
  colnames(gene_info)<- c("Genes","Cromossomo","Posicao","Gene_start","Gene_end","Gene_length","Orientation")
  
  
  ### gene falso positivo - metodos
  top_10_metodo <- read_excel(paste0("Mestrado/Defesa/resultados/resultados_replicas_full/20260728_result_replicas_2etapas/GPCM/top10_",metodo,".xlsx"))
  
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
    file = paste0("verif_prox_genes_errados_pre_sel_",metodo, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}


#################### JUNTA GRM E GPCM COMO PRÉ SELEÇÃO DA TRI PARA VER S, E, F1-SCORE E NUM_GENES_IDENTIF

sensibilidade <- sensibilidade_GRM %>%
  rename(
    `Teste de Burden\nGRM` = burden,
    `SKAT\nGRM` = skat,
    `SKAT-O\nGRM` = skat_o
  ) %>%
  inner_join(
    sensibilidade_GPCM %>%
      rename(
        `Teste de Burden\nGPCM` = burden,
        `SKAT\nGPCM` = skat,
        `SKAT-O\nGPCM` = skat_o
      ),
    by = "replica"
  )

sensibilidade<- sensibilidade[,-c(1)]

#### GRÁFICO DA SENSIBILIDADE


boxplot(
  sensibilidade,
  xaxt = "n",   # remove nomes padrão
  col = "#FF6699",
  ylim=c(0.0,0.7)
  #main='Sensibilidade por método com covariável'
)

text(
  x = 1:ncol(sensibilidade),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(sensibilidade),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)
#summary(sensibilidade)


especificidade <- especificidade_GRM %>%
  rename(
    `Teste de Burden\nGRM` = burden,
    `SKAT\nGRM` = skat,
    `SKAT-O\nGRM` = skat_o
  ) %>%
  inner_join(
    especificidade_GPCM %>%
      rename(
        `Teste de Burden\nGPCM` = burden,
        `SKAT\nGPCM` = skat,
        `SKAT-O\nGPCM` = skat_o
      ),
    by = "replica"
  )

especificidade<- especificidade[,-c(1)]


boxplot(
  especificidade,
  xaxt = "n",   # remove nomes padrão
  col = '#FF7538',
  ylim=c(0.85,1)
  #main='Especificidade por método com covariável'
)

text(
  x = 1:ncol(especificidade),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(especificidade),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)


#### NUM DE GENES IDENTIFICADOS

num_genes <- num_genes_GRM %>%
  rename(
    `Teste de Burden\nGRM` = burden,
    `SKAT\nGRM` = skat,
    `SKAT-O\nGRM` = skat_o
  ) %>%
  inner_join(
    num_genes_GPCM %>%
      rename(
        `Teste de Burden\nGPCM` = burden,
        `SKAT\nGPCM` = skat,
        `SKAT-O\nGPCM` = skat_o
      ),
    by = "replica"
  )

num_genes<- num_genes[,-c(1)]

boxplot(
  num_genes,
  xaxt = "n",   # remove nomes padrão
  col = '#CC66FF',
  ylim=c(0,250)#main='Quantidade de genes identificados por método sem covariável'
)

text(
  x = 1:ncol(num_genes),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(num_genes),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)


### F1 SCORE
f1_score <- f1_score_GRM %>%
  rename(
    `Teste de Burden\nGRM` = burden,
    `SKAT\nGRM` = skat,
    `SKAT-O\nGRM` = skat_o
  ) %>%
  inner_join(
    f1_score_GPCM %>%
      rename(
        `Teste de Burden\nGPCM` = burden,
        `SKAT\nGPCM` = skat,
        `SKAT-O\nGPCM` = skat_o
      ),
    by = "replica"
  )

f1_score<- f1_score[,-c(1)]


boxplot(
  f1_score,
  xaxt = "n",   # remove nomes padrão
  col = "#2A9D8F",
  ylim=c(0.0,0.7)
)

text(
  x = 1:ncol(f1_score),
  y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]),
  labels = colnames(f1_score),
  srt = 45,     # rotação diagonal
  adj = 1,
  xpd = TRUE,
  cex = 0.85
)

#### VENDO QUAIS OS GENES VERDADEIROS FORAM CORRETAMENTE IDENTIFICADOS E QUAIS NÃO FORAM IDENTIFICADOS
library(tidyr)
obter_tp_fn <- function(resultados, simula_gaw17){
  
  tp <- resultados %>%
    inner_join(simula_gaw17, by = c("Genes", "Cromossomo")) %>%
    distinct(replica, Genes, Cromossomo) %>%
    count(Genes, name = "TP")
  
  fn <- resultados %>%
    group_by(replica) %>%
    group_modify(~ anti_join(simula_gaw17, .x,
                             by = c("Genes", "Cromossomo"))) %>%
    ungroup() %>%
    distinct(replica, Genes, Cromossomo) %>%
    count(Genes, name = "FN")
  
  full_join(tp, fn, by = "Genes") %>%
    mutate(
      TP = replace_na(TP, 0),
      FN = replace_na(FN, 0)
    ) %>%
    arrange(desc(TP))
}

freq_burden_GRM <- obter_tp_fn(burden_signif_GRM_categ, simula_gaw17)
#View(freq_burden_GRM)
freq_skat_GRM  <- obter_tp_fn(skat_signif_GRM_categ, simula_gaw17)
View(freq_skat_GRM)
freq_skat_o_GRM  <- obter_tp_fn(skat_o_signif_GRM_categ, simula_gaw17)
View(freq_skat_o_GRM)

freq_burden_GPCM <- obter_tp_fn(burden_signif_GPCM_categ, simula_gaw17)
View(freq_burden_GPCM)
freq_skat_GPCM  <- obter_tp_fn(skat_signif_GPCM_categ, simula_gaw17)
View(freq_skat_GPCM)
freq_skat_o_GPCM  <- obter_tp_fn(skat_o_signif_GPCM_categ, simula_gaw17)



options(scipen = 999)

bases <- c(
  "freq_burden_GRM",
  "freq_skat_GRM",
  "freq_skat_o_GRM",
  "freq_burden_GPCM",
  "freq_skat_GPCM",
  "freq_skat_o_GPCM"
)

for (nome in bases) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0("genes_vdd_signif_",nome, ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}


### JUNTAR TODOS OS GENES CORRETOS COM TP DE TODOS OS MÉTODOS E TANTO GRM QNTO GPCM

resultado <- Reduce(
  function(x, y) merge(x, y, by = "Genes", all = TRUE),
  lapply(bases, function(nome) {
    base <- get(nome)
    dados <- base[, c("Genes", "TP")]
    names(dados)[2] <- nome
    dados
  })
)
View(resultado)


bases <- c(
  "resultado"
)

for (nome in bases) {
  
  df <- get(nome)
  
  write.table(
    df,
    file = paste0("resumo_",nome,"_genes_vdd_signif", ".csv"),
    sep = ";",
    dec = ".",
    row.names = FALSE,
    col.names = TRUE
  )
}
