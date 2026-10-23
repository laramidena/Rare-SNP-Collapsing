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


fenotipos<-read.csv(file="~/Dados_ISA/banco_minicurso_841.csv",header=TRUE) ## essa base será utilizada apenas para fazermos match de indivíduos entre bases, mais abaixo no match
#View(fenotipos)
#### MATRIZ DE CORRELAÇÃO ENTRE OS INDIVÍDUOS

grm_v2<-read.table("~/Dados_ISA/isa.grm841x841.domicilio.grm",sep="\t")
#View(grm_v2)

## MATRIZ DE CORRELAÇÃO (PARA EXCLUIR IND MUITO CORRELACIONADOS)
relatedness_threshold <- 0.354
related_ids <- grm_v2[which((grm_v2$V4 > relatedness_threshold) & (grm_v2$V1!=grm_v2$V2)),] #ind relacionados (mesma familia)
diagonal_ids <- grm_v2[which((grm_v2$V4 > relatedness_threshold) & (grm_v2$V1==grm_v2$V2)),] #todos ind (diag princ)
unrelated_ids <- diagonal_ids[!(diagonal_ids$V2 %in% related_ids$V2), ] #base sem ind relacionados

## COMP PRINC DE ANCESTRALIDADE (ordena as duas bases pelo ID e faz merge delas)
comp<-data.frame(read_dta("~/Dados_ISA/bancoISA_GWAS_n841_2024.dta")[,c(1,3,5,19,24,26,28,30,32,33)]) # dessa base vem os fenótipos
comp2<-data.frame(read_dta("~/Dados_ISA/CPancestralidadeISA_2024.dta")[,c(1,4,5)]) # dessa base vem os componentes principais de ancestralidade
comp<-comp[order(comp[,1]),]
comp2<-comp2[order(comp2[,1]),]
#View(comp_tot) ## merge das duas bases sem o ID agora (tem os componentes principais de ancestralidade - CP1 e CP2, sexo, idade e IMC)
#

#View(comp_tot) #base final com fenótipos e com comp princ de ancestralidade e idade, sexo, etc


## PARA PUXAR O FENÓTIPO NA ORDEM CERTA
feno<-3

#### SNPS COMUNS CR 8
snp_comum<-paste("~/Dados_ISA/SNPs_chr",8,"_total_imput_fam.txt",sep="")
SNPs_comuns<-read.csv(file=snp_comum, header=TRUE, sep=" ")
#View(SNPs_comuns)
names(SNPs_comuns) <- sub("_.*$", "", names(SNPs_comuns))



relevantes<- c("rs12679789","rs17601051","rs72719740","rs940079","rs79159586")


enableJIT(3)
for (cr in c(1:22)){
  comp_tot<-cbind(comp[,-1],comp2[,-1])
  ## FAZ O MATCH DOS INDIVÍDUOS
  match<-NULL
  for (i in 1:nrow(fenotipos)){
    if (!is.na(fenotipos[i,13]) & !is.na(fenotipos[i,19])) {match[i]<-which(round(fenotipos[i,13],2)==round(comp_tot[,4],2) & round(fenotipos[i,19],4)==round(comp_tot[,7],4) & fenotipos[i,3]==comp_tot[,2])} else {if (!is.na(fenotipos[i,13])) match[i]<-which(round(fenotipos[i,13],2)==round(comp_tot[,4],2) & (is.na(comp_tot[,7])) & fenotipos[i,3]==comp_tot[,2]) else match[i]<-which(round(fenotipos[i,19],4)==round(comp_tot[,7],4) & (is.na(comp_tot[,4])) & fenotipos[i,3]==comp_tot[,2])}}
  #	
  
  comp_tot<-as.matrix(comp_tot[match,])
  
  
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
  
  # SUBSTITUIR AQUI
  
  base_snp<-paste("~/Dados_ISA/SNPs_chr",cr,"_total_imput_fam.txt",sep="")
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
  lista_genes <- readRDS("~/Dados_ISA/lista_chr_SNPs_genes_0726.rds")
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
  
  
  #### TENTANDO AGRUPAR OS SNPS Q ESTÃO NO MESMO GENE
  genes<-as.character(data.frame(table(tabela_gene_cr13[,2]))[,1]) #nomes dos genes
  
  
  
  ######## AGORA É O MODELO COM COVARIÁVEIS SELECIONADAS PELO LASSO
  
  
  #### SELECIONANDO SNPS COMUNS PELO LASSO
  set.seed(2511)
  
  
  #### TESTANDO SIGNIFICANCIA DOS GENES DE SNPS RAROS USANDO COMO COVARIÁVEIS OS SNPS COMUNS E VARIÁVEIS SELECIONADAS PELO LASSO
  
  corr_ML_grm <- NULL
  corr_ML_gpcm<- NULL
  valor_p_grm<- NULL
  valor_p_gpcm<- NULL
  matriz_grm<- matrix(0,nrow=nrow(dados_13),ncol=length(genes))
  matriz_gpcm<- matrix(0,nrow=nrow(dados_13),ncol=length(genes))
  ponto_corte_grm<- NULL
  ponto_corte_gpcm<- NULL
  valor_p_arvore_grm<- NULL
  valor_p_arvore_gpcm<- NULL
  a_grm_cov<-NULL
  b1_grm_cov<-NULL
  a_gpcm_cov<-NULL
  b1_gpcm_cov<-NULL
  
  library(mirt)
  set.seed(2511)
  enableJIT(3)
  for (i in 1:length(genes)){ #gene 798 tinha dado problema no cr 2
    snps<-which(tabela_gene_cr13[,2]==genes[i]) #aqui sai a posição na lista genes (tá em ordem alfabética) dos snps q estão no gene pedido da tabela_gene_cr13
    #novo:
    snps_val <- SNPs_raros_cr13[,snps,drop=FALSE]
    #cat(
    #  "Gene:", genes[i],
    #  "| Nº SNPs:", length(snps),
    #  "| ncol:", ncol(snps_val),
    #  "| colnames:", !is.null(colnames(snps_val)),
    #  "\n"
    #)
    #snps_val <- SNPs_raros_cr13[,snps]
    nomes_col_gene <- colnames(SNPs_raros_cr13[which(tabela_gene_cr13[,2]==genes[i])])
    resultado<- try({
      mod_grm <- mirt(snps_val, 1, itemtype = "graded", TOL = 0.0001,
                      technical = list(NCYCLES = 100000))
      mod_gpcm <- mirt(snps_val, 1, itemtype = "gpcm", TOL = 0.0001,
                       technical = list(NCYCLES = 100000))
      coef_grm<- coef(mod_grm, IRTpars = TRUE, simplify = TRUE)
      coef_gpcm<- coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)
    })
    if (inherits(resultado, "try-error")) {
      valor_p_grm[i] <- NA
      valor_p_gpcm[i] <- NA
      ponto_corte_grm[i]<- NA
      ponto_corte_gpcm[i]<- NA
      valor_p_arvore_grm[i]<- NA
      valor_p_arvore_gpcm[i]<- NA
      next
    }
    #if ((cr==22) & (i==50)) {
    #  valor_p_grm[i] <- 'problema'
    #  valor_p_gpcm[i] <- 'problema'
    #  ponto_corte_grm[i]<- 'problema'
    #  ponto_corte_gpcm[i]<- 'problema'
    # valor_p_arvore_grm[i]<- 'problema'
    #  valor_p_arvore_gpcm[i]<- 'problema'
    #  next
    #} 
    else {
      a_grm_cov<-c(a_grm_cov,coef(mod_grm, IRTpars = TRUE, simplify = TRUE)$items[,1])
      b1_grm_cov<-c(b1_grm_cov,coef(mod_grm, IRTpars = TRUE, simplify = TRUE)$items[,2])
      a_gpcm_cov<-c(a_gpcm_cov,coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)$items[,1])
      b1_gpcm_cov<-c(b1_gpcm_cov,coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)$items[,2])
      
      ### COLOQUEI UM TRY
      fs_grm<-try(fscores(mod_grm,full.scores=FALSE, method = 'ML'))
      #fs_grm<-fscores(mod_grm,full.scores=FALSE, method = 'ML')
      ## adicionei isso aqui
      if (inherits(fs_grm, "try-error")) {
        valor_p_grm[i] <- NA
        valor_p_gpcm[i] <- NA
        ponto_corte_grm[i]<- NA
        ponto_corte_gpcm[i]<- NA
        valor_p_arvore_grm[i]<- NA
        valor_p_arvore_gpcm[i]<- NA
        next
      }
      ### a parti disso aqui já tinha
      resp_mirt_df_grm<- as.data.frame(fs_grm)
      resp_mirt_df_grm$ID <- do.call(paste0, resp_mirt_df_grm[ ,nomes_col_gene])
      snps_val$ID <- do.call(paste0, snps_val[ ,nomes_col_gene])
      snps_val$Q1 <- Q2
      left_cr13_grm <- dplyr::left_join(snps_val, dplyr::select(resp_mirt_df_grm, c(ID, F1)), by = "ID")
      
      
      #### SUBSTITUINDO (NA VISÃO INDIVIDUO) QUEM TEM SCORE F1 INF E -INF PELO MAX+1 E MÍN -1
      
      
      min_finito_grm <- min(left_cr13_grm$F1[is.finite(left_cr13_grm$F1)])-1
      max_finito_grm <- max(left_cr13_grm$F1[is.finite(left_cr13_grm$F1)])+1
      left_cr13_grm$F1[left_cr13_grm$F1 == -Inf] <- min_finito_grm
      left_cr13_grm$F1[left_cr13_grm$F1 == Inf] <- max_finito_grm
      
      corr_ML_grm[i]<- cor(left_cr13_grm$Q1,left_cr13_grm$F1)
      
      if (all(is.infinite(left_cr13_grm$F1))){
        valor_p_grm[i] <- NA
        valor_p_gpcm[i] <- NA
        ponto_corte_grm[i]<- NA
        ponto_corte_gpcm[i]<- NA
        valor_p_arvore_grm[i]<- NA
        valor_p_arvore_gpcm[i]<- NA
        next
      }
      
      ajuste_grm<- lm(left_cr13_grm$Q1~(left_cr13_grm$F1 + var_comum_sign))
      #ajuste_grm<- lm(left_cr13_grm$Q1~as.factor(left_cr13_grm$F1_cat)+ var_comum_sign)
      
      matriz_grm[,i]<-left_cr13_grm$F1
      valor_p_grm[i]<- anova(ajuste_grm)[[5]][1]
      
      ### categórico grm
      
      arvore_reg_grm <- rpart(left_cr13_grm$Q1~(left_cr13_grm$F1),  method = "anova", minbucket=10,maxdepth=1)
      
      
      if (!is.null(arvore_reg_grm$splits)) {
        ponto_corte_grm[i] <- arvore_reg_grm$splits[,"index"][1]
        left_cr13_grm$F1_cat <- ifelse(left_cr13_grm$F1 < as.numeric(ponto_corte_grm[i]), "negativo","positivo")
        left_cr13_grm$F1_cat <- factor(left_cr13_grm$F1_cat, levels = c("negativo", "positivo"))
        ajuste_arvore_grm <- lm(left_cr13_grm$Q1~(left_cr13_grm$F1_cat) + var_comum_sign)
        valor_p_arvore_grm[i]<- anova(ajuste_arvore_grm)[[5]][1]
        
      } else {
        ponto_corte_grm[i] <- 'sem corte'  # ou qualquer outro valor padrão
        valor_p_arvore_grm[i] <- 'sem corte'
      }
      
      
      #### AGORA É GPCM
      fs_gpcm<-fscores(mod_gpcm,full.scores=FALSE, method = 'ML')
      resp_mirt_df_gpcm<- as.data.frame(fs_gpcm)
      resp_mirt_df_gpcm$ID <- do.call(paste0, resp_mirt_df_gpcm[ ,nomes_col_gene])
      snps_val$ID <- do.call(paste0, snps_val[ ,nomes_col_gene])
      snps_val$Q1 <- Q2
      left_cr13_gpcm <- dplyr::left_join(snps_val, dplyr::select(resp_mirt_df_gpcm, c(ID, F1)), by = "ID")
      
      
      
      
      ## substitui -inf no F1 pelo [menor valor (sem ser o -inf) -1] e
      # +inf no F1 por [maior valor (sem ser +inf) +1]
      min_finito_gpcm <- min(left_cr13_gpcm$F1[is.finite(left_cr13_gpcm$F1)])-1
      max_finito_gpcm <- max(left_cr13_gpcm$F1[is.finite(left_cr13_gpcm$F1)])+1
      left_cr13_gpcm$F1[left_cr13_gpcm$F1 == -Inf] <- min_finito_gpcm
      left_cr13_gpcm$F1[left_cr13_gpcm$F1 == Inf] <- max_finito_gpcm
      
      
      
      #corr_gpcm_WLE[i]<- cor(left_cr13_gpcm$Q1,left_cr13_gpcm$F1)
      corr_ML_gpcm[i]<- cor(left_cr13_gpcm$Q1,left_cr13_gpcm$F1)
      
      #ajuste_gpcm<- lm(left_cr13_gpcm$Q1~as.factor(left_cr13_gpcm$F1_cat) + var_comum_sign)
      ajuste_gpcm<- lm(left_cr13_gpcm$Q1~(left_cr13_gpcm$F1 + var_comum_sign))
      
      matriz_gpcm[,i]<-left_cr13_gpcm$F1
      valor_p_gpcm[i]<- anova(ajuste_gpcm)[[5]][1]
      
      ### categórico gpcm
      
      arvore_reg_gpcm <- rpart(left_cr13_gpcm$Q1~(left_cr13_gpcm$F1),  method = "anova", minbucket=10,maxdepth=1)
      
      #ponto_corte_gpcm[i]<- arvore_reg_gpcm$splits[,"index"][1]
      
      if (!is.null(arvore_reg_gpcm$splits)) {
        ponto_corte_gpcm[i] <- arvore_reg_gpcm$splits[,"index"][1]
        left_cr13_gpcm$F1_cat <- ifelse(left_cr13_gpcm$F1 < as.numeric(ponto_corte_gpcm[i]), "negativo","positivo")
        left_cr13_gpcm$F1_cat <- factor(left_cr13_gpcm$F1_cat, levels = c("negativo", "positivo"))
        ajuste_arvore_gpcm <- lm(left_cr13_gpcm$Q1~as.factor(left_cr13_gpcm$F1_cat) + var_comum_sign)
        valor_p_arvore_gpcm[i]<- anova(ajuste_arvore_gpcm)[[5]][1]
        
      } else {
        ponto_corte_gpcm[i] <- 'sem corte'  # ou qualquer outro valor padrão
        valor_p_arvore_gpcm[i] <- 'sem corte'
      }
    }
  }
  #valor_p_arvore_grm
  #valor_p_arvore_gpcm
  #valor_p_grm
  #valor_p_gpcm
  
  
  
  
  
  
  
  sign_gene_cr13<-cbind(genes, valor_p_grm,valor_p_gpcm, valor_p_arvore_grm, valor_p_arvore_gpcm, cr)
  colnames(sign_gene_cr13)<- c("Genes", "Valor-p GRM C COV","Valor-p GPCM C COV", "Valor-p Arvore GRM C COV", "Valor-p Arvore GPCM C COV", "Cromossomo")
  
  sign_gene_cr13_cov_df <- data.frame(sign_gene_cr13)
  #View(sign_gene_cr13)
  
  ## NOVO:
  
  alpha_grm<- rep((0.05/(length(genes))),length(genes))
  alpha_gpcm<- rep((0.05/(length(genes))),length(genes))
  
  
  sign_gene_cr13<-cbind(genes, valor_p_grm,valor_p_gpcm, valor_p_arvore_grm, valor_p_arvore_gpcm, alpha_grm, alpha_gpcm, cr)
  colnames(sign_gene_cr13)<- c("Genes", "Valor-p GRM C COV","Valor-p GPCM C COV", "Valor-p Arvore GRM C COV", "Valor-p Arvore GPCM C COV", "Nivel de sign grm", "Nível de sign gpcm", "Cromossomo")
  #View(sign_gene_cr13)
  ##### UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM GENE RARO É SIGNIFICATIVO
  ind_gene_raro_sign_bf_com_cov<- which(as.numeric(sign_gene_cr13[,2])<(0.05/(length(genes)))|
                                          as.numeric(sign_gene_cr13[,3])<(0.05/(length(genes)))|
                                          as.numeric(sign_gene_cr13[,4])<(0.05/(length(genes)))|
                                          as.numeric(sign_gene_cr13[,5])<(0.05/(length(genes))))
  
  
  gene_raro_sign_bf_com_cov<- as.matrix(sign_gene_cr13[ind_gene_raro_sign_bf_com_cov,])
  #View(gene_raro_sign_bf_com_cov)
  
  
  
  
  cat('Cromossomo',cr,"-",'rodou com sucesso',"\n",file=paste("flag_successo_dados_ISA_met_TRI_v2",sep=""),append=T)
  
  
  
  ################ EXPORTANDO
  
  #View(gene_raro_sign_bf_com_cov) #resultado gene a gene
  #View(gene_raro_sign_bf_com_cov_agrup_final) #resultado com genes agrupados para SNPs q não rodavam
  
  
  
  if (dim(as.matrix(gene_raro_sign_bf_com_cov))[2]==1){
    gene_raro_sign_bf_com_cov<- matrix(gene_raro_sign_bf_com_cov,nrow=1)
  }
  #cat('',gene_raro_sign_bf_com_cov,file=paste("gene_raro_com_cov_repli",repli,".txt",sep=""),append=T)
  write.table(
    gene_raro_sign_bf_com_cov,
    file = paste0("gene_raro_com_cov_TRI_dados_ISA_v2.txt"),
    sep = ",",
    row.names = FALSE,
    col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
    append = (cr != 1)
  )
  
  
  
  
  #### exportando os coef a e b 
  tabelas_a = list(summary(a_grm_cov),
                   summary(a_gpcm_cov))
  
  nomes_tabelas_a = c("Summary a grm com cov",
                      "Summary a gpcm com cov")
  
  cat('Cromossomo',cr,"\n",file=paste("summary_coef_a_dados_ISA_v2.txt",sep=""),append=T)
  for (i in seq_along(tabelas_a)) {
    cat(nomes_tabelas_a[i], "\n", file=paste("summary_coef_a_dados_ISA_v2.txt",sep=""),append=T)
    capture.output(print(tabelas_a[[i]], row.names = FALSE), file=paste("summary_coef_a_dados_ISA_v2.txt",sep=""),append=T)
    cat("\n", file = paste("summary_coef_a_dados_ISA_v2.txt",sep=""),append=T)
  }
  # Linha em branco extra entre cromossomos
  cat("\n", file = paste("summary_coef_a_dados_ISA_v2.txt",sep=""),append=T)
  
  
  #### exportando os coef a e b 
  tabelas_b = list(summary(b1_grm_cov),
                   summary(b1_gpcm_cov))
  
  nomes_tabelas_b = c("Summary b1 grm com cov",
                      "Summary b1 gpcm com cov")
  
  cat('Cromossomo',cr,"\n",file=paste("summary_coef_b1_dados_ISA_v2.txt",sep=""),append=T)
  for (i in seq_along(tabelas_b)) {
    cat(nomes_tabelas_b[i], "\n", file=paste("summary_coef_b1_dados_ISA_v2.txt",sep=""),append=T)
    capture.output(print(tabelas_b[[i]], row.names = FALSE), file=paste("summary_coef_b1_dados_ISA_v2.txt",sep=""),append=T)
    cat("\n", file=paste("summary_coef_b1_dados_ISA_v2.txt",sep=""),append=T)
  }
  # Linha em branco extra entre cromossomos
  cat("\n", file=paste("summary_coef_b1_dados_ISA_v2.txt",sep=""),append=T)
  
  
  ########################################################### À PARTIR DAQUI É CÓDIGO NOVO
  
  ######################################## AGORA É AGRUPADO
  qtd_sucessos_grm_cov<-length(which(is.na(valor_p_grm)==FALSE))
  qtd_sucessos_gpcm_cov<- length(which(is.na(valor_p_gpcm)==FALSE))
  
  sign_gene_cr13_cov_df <- data.frame(sign_gene_cr13) 
  
  pos_snps_pouco_var <- which((is.na(sign_gene_cr13_cov_df$Valor.p.GRM.C.COV)==TRUE) | (is.na(sign_gene_cr13_cov_df$Valor.p.GPCM.C.COV)==TRUE)| (sign_gene_cr13_cov_df$Valor.p.GRM.C.COV=='problema') | (sign_gene_cr13_cov_df$Valor.p.GPCM.C.COV=='problema'))
  #pos_snps_pouco_var ### genes com pouca variabilidade
  genes_na <- as.data.frame(sign_gene_cr13_cov_df[pos_snps_pouco_var,1])
  colnames(genes_na)<-c("gene")
  
  snps_filtrados <- tabela_gene_cr13[tabela_gene_cr13$gene %in% genes_na$gene, ]
  ### PEGANDO SÓ OS NOMES DOS SNPS
  snps_dos_genes_na <- snps_filtrados$SNP
  snps_agrup_cov <- SNPs_raros_cr13[, snps_dos_genes_na, drop = FALSE]
  nomes_agrup<- colnames(snps_agrup_cov)
  valor_p_agrup_grm<- NULL
  valor_p_arvore_grm_agrup<- NULL
  ponto_corte_grm_agrup<- NULL
  valor_p_agrup_gpcm<- NULL
  valor_p_arvore_gpcm_agrup<- NULL
  ponto_corte_gpcm_agrup<- NULL
  a_grm_cov_agrup<- NULL
  b1_grm_cov_agrup<- NULL
  a_gpcm_cov_agrup<- NULL
  b1_gpcm_cov_agrup<- NULL
  agrup<- NULL
  enableJIT(3)
  if (ncol(snps_agrup_cov)>20){
    seque<- seq(1,ncol(snps_agrup_cov), by = 10)
    for(i in 1:(length(seque)-1)){
      #for(i in 18:(length(seque)-1)){
      snps_agrupados<- snps_agrup_cov[,seque[i]:(seque[i+1]-1)]
      mod_grm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "graded", TOL = 0.0001,
                                technical = list(NCYCLES = 100000)),silent=TRUE)
      mod_gpcm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "gpcm", TOL = 0.0001,
                                 technical = list(NCYCLES = 100000)),silent=TRUE)
      nomes_agrup<- colnames(snps_agrupados)
      agrup<- c(agrup,snps_filtrados[seque[i]:(seque[i+1]-1),2])
      
      #if ((cr==11) & (i==17)) {
      #  valor_p_agrup_grm[i]<- 'problema'
      #  valor_p_arvore_grm_agrup[i]<- 'problema'
      #  ponto_corte_grm_agrup[i]<- 'problema'
      #  a_grm_cov_agrup[i]<- 'problema'
      #  b1_grm_cov_agrup[i]<- 'problema'
      #  next
      #}
      if (inherits(mod_grm_agrup, "try-error")) {
        valor_p_agrup_grm[i]<- 'funcao falhou'
        valor_p_arvore_grm_agrup[i]<- 'funcao falhou'
        ponto_corte_grm_agrup[i]<- 'funcao falhou'
        a_grm_cov_agrup[i]<- 'funcao falhou'
        b1_grm_cov_agrup[i]<- 'funcao falhou'
        
      }
      else{
        a_grm_cov_agrup<- c(a_grm_cov_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_grm_cov_agrup<-c(b1_grm_cov_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
        fs_grm_agrup<-fscores(mod_grm_agrup,full.scores=FALSE, method = 'ML')
        resp_mirt_df_grm_agrup<- as.data.frame(fs_grm_agrup)
        resp_mirt_df_grm_agrup$ID <- do.call(paste0, resp_mirt_df_grm_agrup[ ,nomes_agrup])
        snps_agrup_cov$ID <- do.call(paste0, snps_agrup_cov[ ,nomes_agrup])
        snps_agrup_cov$Q1 <- Q2
        
        left_cr13_grm_agrup <- dplyr::left_join(snps_agrup_cov, dplyr::select(resp_mirt_df_grm_agrup, c(ID, F1)), by = "ID")
        
        #### SUBSTITUINDO (NA VISÃO INDIVIDUO) QUEM TEM SCORE F1 INF E -INF PELO MAX+1 E MÍN -1
        
        
        min_finito_grm_agrup <- min(left_cr13_grm_agrup$F1[is.finite(left_cr13_grm_agrup$F1)])-1
        max_finito_grm_agrup <- max(left_cr13_grm_agrup$F1[is.finite(left_cr13_grm_agrup$F1)])+1
        left_cr13_grm_agrup$F1[left_cr13_grm_agrup$F1 == -Inf] <- min_finito_grm_agrup
        left_cr13_grm_agrup$F1[left_cr13_grm_agrup$F1 == Inf] <- max_finito_grm_agrup
        
        if ((all(is.infinite(left_cr13_grm_agrup$F1)))| ((is.na(min_finito_grm_agrup)) & (is.infinite(max_finito_grm_agrup))) | ((is.na(max_finito_grm_agrup)) & (is.infinite(min_finito_grm_agrup))) | ((is.infinite(max_finito_grm_agrup)) & (is.infinite(min_finito_grm_agrup)))){
          valor_p_agrup_grm[i]<- NA
          valor_p_arvore_grm_agrup[i] <- NA
          ponto_corte_grm_agrup[i] <- NA
        }
        else{  
          ajuste_grm_agrup<- lm(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1+ var_comum_sign))
          
          valor_p_agrup_grm[i]<- anova(ajuste_grm_agrup)[[5]][1]
          
          
          arvore_reg_grm_agrup<- rpart(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
          
          if (!is.null(arvore_reg_grm_agrup$splits)) {
            ponto_corte_grm_agrup[i] <- arvore_reg_grm_agrup$splits[,"index"][1]
            left_cr13_grm_agrup$F1_cat <- ifelse(left_cr13_grm_agrup$F1 < as.numeric(ponto_corte_grm_agrup[i]), "negativo","positivo")
            left_cr13_grm_agrup$F1_cat <- factor(left_cr13_grm_agrup$F1_cat, levels = c("negativo", "positivo"))
            ajuste_arvore_grm_agrup <- lm(left_cr13_grm_agrup$Q1~as.factor(left_cr13_grm_agrup$F1_cat) + var_comum_sign)
            valor_p_arvore_grm_agrup[i]<- anova(ajuste_arvore_grm_agrup)[[5]][1]
            
          } else {
            ponto_corte_grm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
            valor_p_arvore_grm_agrup[i] <- 'sem corte'
          }
          
        }
      }
      if (inherits(mod_gpcm_agrup, "try-error")) {
        valor_p_agrup_gpcm[i]<- 'funcao falhou'
        valor_p_arvore_gpcm_agrup[i]<- 'funcao falhou'
        ponto_corte_gpcm_agrup[i] <- 'funcao falhou'
        a_gpcm_cov_agrup[i]<- 'funcao falhou'
        b1_gpcm_cov_agrup[i]<- 'funcao falhou'
        
      } else {
        a_gpcm_cov_agrup<- c(a_gpcm_cov_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_gpcm_cov_agrup<-c(b1_gpcm_cov_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
        
        
        #### AGORA É GPCM
        fs_gpcm_agrup<-fscores(mod_gpcm_agrup,full.scores=FALSE, method = 'ML')
        resp_mirt_df_gpcm_agrup<- as.data.frame(fs_gpcm_agrup)
        resp_mirt_df_gpcm_agrup$ID <- do.call(paste0, resp_mirt_df_gpcm_agrup[ ,nomes_agrup])
        snps_agrup_cov$ID <- do.call(paste0, snps_agrup_cov[ ,nomes_agrup])
        snps_agrup_cov$Q1 <- Q2
        left_cr13_gpcm_agrup <- dplyr::left_join(snps_agrup_cov, dplyr::select(resp_mirt_df_gpcm_agrup, c(ID, F1)), by = "ID")
        
        
        
        ## substitui -inf no F1 pelo [menor valor (sem ser o -inf) -1] e
        # +inf no F1 por [maior valor (sem ser +inf) +1]
        min_finito_gpcm_agrup <- min(left_cr13_gpcm_agrup$F1[is.finite(left_cr13_gpcm_agrup$F1)])-1
        max_finito_gpcm_agrup <- max(left_cr13_gpcm_agrup$F1[is.finite(left_cr13_gpcm_agrup$F1)])+1
        left_cr13_gpcm_agrup$F1[left_cr13_gpcm_agrup$F1 == -Inf] <- min_finito_gpcm_agrup
        left_cr13_gpcm_agrup$F1[left_cr13_gpcm_agrup$F1 == Inf] <- max_finito_gpcm_agrup
        
        if ((all(is.infinite(left_cr13_gpcm_agrup$F1)))| ((is.na(min_finito_gpcm_agrup)) & (is.infinite(max_finito_gpcm_agrup))) | ((is.na(max_finito_gpcm_agrup)) & (is.infinite(min_finito_gpcm_agrup))) | ((is.infinite(max_finito_gpcm_agrup)) & (is.infinite(min_finito_gpcm_agrup)))){
          valor_p_agrup_gpcm[i]<- NA
          ponto_corte_gpcm_agrup[i] <- NA
          valor_p_arvore_gpcm_agrup[i] <- NA
        }
        else{  
          ajuste_gpcm_agrup<- lm(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1 + var_comum_sign))
          valor_p_agrup_gpcm[i]<- anova(ajuste_gpcm_agrup)[[5]][1]
          
          
          arvore_reg_gpcm_agrup <- rpart(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
          
          if (!is.null(arvore_reg_gpcm_agrup$splits)) {
            ponto_corte_gpcm_agrup[i] <- arvore_reg_gpcm_agrup$splits[,"index"][1]
            left_cr13_gpcm_agrup$F1_cat <- ifelse(left_cr13_gpcm_agrup$F1 < as.numeric(ponto_corte_gpcm_agrup[i]), "negativo","positivo")
            left_cr13_gpcm_agrup$F1_cat <- factor(left_cr13_gpcm_agrup$F1_cat, levels = c("negativo", "positivo"))
            ajuste_arvore_gpcm_agrup <- lm(left_cr13_gpcm_agrup$Q1~as.factor(left_cr13_gpcm_agrup$F1_cat) + var_comum_sign)
            valor_p_arvore_gpcm_agrup[i]<- anova(ajuste_arvore_gpcm_agrup)[[5]][1]
            
          } else {
            ponto_corte_gpcm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
            valor_p_arvore_gpcm_agrup[i] <- 'sem corte'
          }
        } 
      }  
    }
    i=i+1
    snps_agrupados<- snps_agrup_cov[,seque[i]:(ncol(snps_agrup_cov)-2)]
    mod_grm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "graded", TOL = 0.0001,
                              technical = list(NCYCLES = 100000)),silent=TRUE)
    mod_gpcm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "gpcm", TOL = 0.0001,
                               technical = list(NCYCLES = 100000)),silent=TRUE)
    nomes_agrup<- colnames(snps_agrupados)
    #agrup<- snps_filtrados[seque[i]:seque[i+1]-1,2]
    agrup<- c(agrup,snps_filtrados[seque[i]:(ncol(snps_agrup_cov)-2),2])
    
    
    if (inherits(mod_grm_agrup, "try-error")) {
      valor_p_agrup_grm[i]<- 'funcao falhou'
      valor_p_arvore_grm_agrup[i]<- 'funcao falhou'
      ponto_corte_grm_agrup[i]<- 'funcao falhou'
      a_grm_cov_agrup<- 'funcao falhou'
      b1_grm_cov_agrup<- 'funcao falhou'
    } else{
      a_grm_cov_agrup<- c(a_grm_cov_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
      b1_grm_cov_agrup<-c(b1_grm_cov_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
      fs_grm_agrup<-fscores(mod_grm_agrup,full.scores=FALSE, method = 'ML')
      resp_mirt_df_grm_agrup<- as.data.frame(fs_grm_agrup)
      resp_mirt_df_grm_agrup$ID <- do.call(paste0, resp_mirt_df_grm_agrup[ ,nomes_agrup])
      snps_agrup_cov$ID <- do.call(paste0, snps_agrup_cov[ ,nomes_agrup])
      snps_agrup_cov$Q1 <- Q2
      
      left_cr13_grm_agrup <- dplyr::left_join(snps_agrup_cov, dplyr::select(resp_mirt_df_grm_agrup, c(ID, F1)), by = "ID")
      
      #### SUBSTITUINDO (NA VISÃO INDIVIDUO) QUEM TEM SCORE F1 INF E -INF PELO MAX+1 E MÍN -1
      
      
      min_finito_grm_agrup <- min(left_cr13_grm_agrup$F1[is.finite(left_cr13_grm_agrup$F1)])-1
      max_finito_grm_agrup <- max(left_cr13_grm_agrup$F1[is.finite(left_cr13_grm_agrup$F1)])+1
      left_cr13_grm_agrup$F1[left_cr13_grm_agrup$F1 == -Inf] <- min_finito_grm_agrup
      left_cr13_grm_agrup$F1[left_cr13_grm_agrup$F1 == Inf] <- max_finito_grm_agrup
      
      if ((all(is.infinite(left_cr13_grm_agrup$F1)))| ((is.na(min_finito_grm_agrup)) & (is.infinite(max_finito_grm_agrup))) | ((is.na(max_finito_grm_agrup)) & (is.infinite(min_finito_grm_agrup))) | ((is.infinite(max_finito_grm_agrup)) & (is.infinite(min_finito_grm_agrup)))){
        valor_p_agrup_grm[i]<- NA
        valor_p_arvore_grm_agrup[i] <- NA
        ponto_corte_grm_agrup[i] <- NA
      }
      else{  
        ajuste_grm_agrup<- lm(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1+ var_comum_sign))
        
        valor_p_agrup_grm[i]<- anova(ajuste_grm_agrup)[[5]][1]
        
        
        arvore_reg_grm_agrup<- rpart(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
        
        if (!is.null(arvore_reg_grm_agrup$splits)) {
          ponto_corte_grm_agrup[i] <- arvore_reg_grm_agrup$splits[,"index"][1]
          left_cr13_grm_agrup$F1_cat <- ifelse(left_cr13_grm_agrup$F1 < as.numeric(ponto_corte_grm_agrup[i]), "negativo","positivo")
          left_cr13_grm_agrup$F1_cat <- factor(left_cr13_grm_agrup$F1_cat, levels = c("negativo", "positivo"))
          ajuste_arvore_grm_agrup <- lm(left_cr13_grm_agrup$Q1~as.factor(left_cr13_grm_agrup$F1_cat) + var_comum_sign)
          valor_p_arvore_grm_agrup[i]<- anova(ajuste_arvore_grm_agrup)[[5]][1]
        } else {
          ponto_corte_grm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
          valor_p_arvore_grm_agrup[i] <- 'sem corte'
        }
      }
    }
    if (inherits(mod_gpcm_agrup, "try-error")) {
      valor_p_agrup_gpcm[i]<- 'funcao falhou'
      valor_p_arvore_gpcm_agrup[i]<- 'funcao falhou'
      ponto_corte_gpcm_agrup[i] <- 'funcao falhou'
      a_gpcm_cov_agrup<- 'funcao falhou'
      b1_gpcm_cov_agrup<- 'funcao falhou'
      
    } else {
      a_gpcm_cov_agrup<- c(a_gpcm_cov_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
      b1_gpcm_cov_agrup<-c(b1_gpcm_cov_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
      
      
      #### AGORA É GPCM
      fs_gpcm_agrup<-fscores(mod_gpcm_agrup,full.scores=FALSE, method = 'ML')
      resp_mirt_df_gpcm_agrup<- as.data.frame(fs_gpcm_agrup)
      resp_mirt_df_gpcm_agrup$ID <- do.call(paste0, resp_mirt_df_gpcm_agrup[ ,nomes_agrup])
      snps_agrup_cov$ID <- do.call(paste0, snps_agrup_cov[ ,nomes_agrup])
      snps_agrup_cov$Q1 <- Q2
      left_cr13_gpcm_agrup <- dplyr::left_join(snps_agrup_cov, dplyr::select(resp_mirt_df_gpcm_agrup, c(ID, F1)), by = "ID")
      
      
      
      ## substitui -inf no F1 pelo [menor valor (sem ser o -inf) -1] e
      # +inf no F1 por [maior valor (sem ser +inf) +1]
      min_finito_gpcm_agrup <- min(left_cr13_gpcm_agrup$F1[is.finite(left_cr13_gpcm_agrup$F1)])-1
      max_finito_gpcm_agrup <- max(left_cr13_gpcm_agrup$F1[is.finite(left_cr13_gpcm_agrup$F1)])+1
      left_cr13_gpcm_agrup$F1[left_cr13_gpcm_agrup$F1 == -Inf] <- min_finito_gpcm_agrup
      left_cr13_gpcm_agrup$F1[left_cr13_gpcm_agrup$F1 == Inf] <- max_finito_gpcm_agrup
      
      if ((all(is.infinite(left_cr13_gpcm_agrup$F1)))| ((is.na(min_finito_gpcm_agrup)) & (is.infinite(max_finito_gpcm_agrup))) | ((is.na(max_finito_gpcm_agrup)) & (is.infinite(min_finito_gpcm_agrup))) | ((is.infinite(max_finito_gpcm_agrup)) & (is.infinite(min_finito_gpcm_agrup)))){
        valor_p_agrup_gpcm[i]<- NA
        ponto_corte_gpcm_agrup[i] <- NA
        valor_p_arvore_gpcm_agrup[i] <- NA
      } else{  
        ajuste_gpcm_agrup<- lm(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1 + var_comum_sign))
        valor_p_agrup_gpcm[i]<- anova(ajuste_gpcm_agrup)[[5]][1]
        
        
        arvore_reg_gpcm_agrup <- rpart(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
        
        if (!is.null(arvore_reg_gpcm_agrup$splits)) {
          ponto_corte_gpcm_agrup[i] <- arvore_reg_gpcm_agrup$splits[,"index"][1]
          left_cr13_gpcm_agrup$F1_cat <- ifelse(left_cr13_gpcm_agrup$F1 < as.numeric(ponto_corte_gpcm_agrup[i]), "negativo","positivo")
          left_cr13_gpcm_agrup$F1_cat <- factor(left_cr13_gpcm_agrup$F1_cat, levels = c("negativo", "positivo"))
          ajuste_arvore_gpcm_agrup <- lm(left_cr13_gpcm_agrup$Q1~as.factor(left_cr13_gpcm_agrup$F1_cat) + var_comum_sign)
          valor_p_arvore_gpcm_agrup[i]<- anova(ajuste_arvore_gpcm_agrup)[[5]][1]
          
        } else {
          ponto_corte_gpcm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
          valor_p_arvore_gpcm_agrup[i] <- 'sem corte'
        }
      } 
    }
  }
  valor_p_agrup_gpcm
  valor_p_agrup_gpcm
  valor_p_arvore_grm_agrup
  valor_p_arvore_gpcm_agrup
  
  #snps_agrup_cov <-snps_agrup_cov[,161:175]
  #agrup<- c(agrup,snps_filtrados[161:175,2])
  ###############
  if (ncol(snps_agrup_cov)<=20){
    i=1
    seque=1
    set.seed(2511)
    mod_grm_agrup <- try(mirt(snps_agrup_cov, 1, itemtype = "graded", TOL = 0.0001,
                              technical = list(NCYCLES = 100000)),silent=TRUE)
    mod_gpcm_agrup <- try(mirt(snps_agrup_cov, 1, itemtype = "gpcm", TOL = 0.0001,
                               technical = list(NCYCLES = 100000)),silent=TRUE)
    nomes_agrup<- colnames(snps_agrup_cov)
    agrup<- c(agrup,snps_filtrados[,2])
    
    if (inherits(mod_grm_agrup, "try-error")) {
      valor_p_agrup_grm[i]<- 'funcao falhou'
      valor_p_arvore_grm_agrup[i]<- 'funcao falhou'
      ponto_corte_grm_agrup[i]<- 'funcao falhou'
      a_grm_cov_agrup<- 'funcao falhou'
      b1_grm_cov_agrup<- 'funcao falhou'
      
    }
    else{
      a_grm_cov_agrup<- c(a_grm_cov_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
      b1_grm_cov_agrup<-c(b1_grm_cov_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
      fs_grm_agrup<-fscores(mod_grm_agrup,full.scores=FALSE, method = 'ML')
      resp_mirt_df_grm_agrup<- as.data.frame(fs_grm_agrup)
      resp_mirt_df_grm_agrup$ID <- do.call(paste0, resp_mirt_df_grm_agrup[ ,nomes_agrup])
      snps_agrup_cov$ID <- do.call(paste0, snps_agrup_cov[ ,nomes_agrup])
      snps_agrup_cov$Q1 <- Q2
      
      left_cr13_grm_agrup <- dplyr::left_join(snps_agrup_cov, dplyr::select(resp_mirt_df_grm_agrup, c(ID, F1)), by = "ID")
      
      #### SUBSTITUINDO (NA VISÃO INDIVIDUO) QUEM TEM SCORE F1 INF E -INF PELO MAX+1 E MÍN -1
      
      
      min_finito_grm_agrup <- min(left_cr13_grm_agrup$F1[is.finite(left_cr13_grm_agrup$F1)])-1
      max_finito_grm_agrup <- max(left_cr13_grm_agrup$F1[is.finite(left_cr13_grm_agrup$F1)])+1
      left_cr13_grm_agrup$F1[left_cr13_grm_agrup$F1 == -Inf] <- min_finito_grm_agrup
      left_cr13_grm_agrup$F1[left_cr13_grm_agrup$F1 == Inf] <- max_finito_grm_agrup
      
      if ((all(is.infinite(left_cr13_grm_agrup$F1)))| ((is.na(min_finito_grm_agrup)) & (is.infinite(max_finito_grm_agrup))) | ((is.na(max_finito_grm_agrup)) & (is.infinite(min_finito_grm_agrup))) | ((is.infinite(max_finito_grm_agrup)) & (is.infinite(min_finito_grm_agrup)))){
        valor_p_agrup_grm[i]<- NA
        valor_p_arvore_grm_agrup[i] <- NA
        ponto_corte_grm_agrup[i] <- NA
      }
      else{  
        ajuste_grm_agrup<- lm(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1+ var_comum_sign))
        
        valor_p_agrup_grm[i]<- anova(ajuste_grm_agrup)[[5]][1]
        
        
        arvore_reg_grm_agrup<- rpart(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
        
        if (!is.null(arvore_reg_grm_agrup$splits)) {
          ponto_corte_grm_agrup[i] <- arvore_reg_grm_agrup$splits[,"index"][1]
          left_cr13_grm_agrup$F1_cat <- ifelse(left_cr13_grm_agrup$F1 < as.numeric(ponto_corte_grm_agrup[i]), "negativo","positivo")
          left_cr13_grm_agrup$F1_cat <- factor(left_cr13_grm_agrup$F1_cat, levels = c("negativo", "positivo"))
          ajuste_arvore_grm_agrup <- lm(left_cr13_grm_agrup$Q1~as.factor(left_cr13_grm_agrup$F1_cat) + var_comum_sign)
          valor_p_arvore_grm_agrup[i]<- anova(ajuste_arvore_grm_agrup)[[5]][1]
          
        } else {
          ponto_corte_grm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
          valor_p_arvore_grm_agrup[i] <- 'sem corte'
        }
        
      }
    }
    if (inherits(mod_gpcm_agrup, "try-error")) {
      valor_p_agrup_gpcm[i]<- 'funcao falhou'
      valor_p_arvore_gpcm_agrup[i]<- 'funcao falhou'
      ponto_corte_gpcm_agrup[i] <- 'funcao falhou'
      a_gpcm_cov_agrup<- 'funcao falhou'
      b1_gpcm_cov_agrup<- 'funcao falhou'
      
    } else {
      a_gpcm_cov_agrup<- c(a_gpcm_cov_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
      b1_gpcm_cov_agrup<-c(b1_gpcm_cov_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
      
      
      #### AGORA É GPCM
      fs_gpcm_agrup<-fscores(mod_gpcm_agrup,full.scores=FALSE, method = 'ML')
      resp_mirt_df_gpcm_agrup<- as.data.frame(fs_gpcm_agrup)
      resp_mirt_df_gpcm_agrup$ID <- do.call(paste0, resp_mirt_df_gpcm_agrup[ ,nomes_agrup])
      snps_agrup_cov$ID <- do.call(paste0, snps_agrup_cov[ ,nomes_agrup])
      snps_agrup_cov$Q1 <- Q2
      left_cr13_gpcm_agrup <- dplyr::left_join(snps_agrup_cov, dplyr::select(resp_mirt_df_gpcm_agrup, c(ID, F1)), by = "ID")
      
      
      
      ## substitui -inf no F1 pelo [menor valor (sem ser o -inf) -1] e
      # +inf no F1 por [maior valor (sem ser +inf) +1]
      min_finito_gpcm_agrup <- min(left_cr13_gpcm_agrup$F1[is.finite(left_cr13_gpcm_agrup$F1)])-1
      max_finito_gpcm_agrup <- max(left_cr13_gpcm_agrup$F1[is.finite(left_cr13_gpcm_agrup$F1)])+1
      left_cr13_gpcm_agrup$F1[left_cr13_gpcm_agrup$F1 == -Inf] <- min_finito_gpcm_agrup
      left_cr13_gpcm_agrup$F1[left_cr13_gpcm_agrup$F1 == Inf] <- max_finito_gpcm_agrup
      
      if ((all(is.infinite(left_cr13_gpcm_agrup$F1)))| ((is.na(min_finito_gpcm_agrup)) & (is.infinite(max_finito_gpcm_agrup))) | ((is.na(max_finito_gpcm_agrup)) & (is.infinite(min_finito_gpcm_agrup))) | ((is.infinite(max_finito_gpcm_agrup)) & (is.infinite(min_finito_gpcm_agrup)))){
        valor_p_agrup_gpcm[i]<- NA
        ponto_corte_gpcm_agrup[i] <- NA
        valor_p_arvore_gpcm_agrup[i] <- NA
      }
      else{  
        ajuste_gpcm_agrup<- lm(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1 + var_comum_sign))
        valor_p_agrup_gpcm[i]<- anova(ajuste_gpcm_agrup)[[5]][1]
        
        
        arvore_reg_gpcm_agrup <- rpart(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
        
        if (!is.null(arvore_reg_gpcm_agrup$splits)) {
          ponto_corte_gpcm_agrup[i] <- arvore_reg_gpcm_agrup$splits[,"index"][1]
          left_cr13_gpcm_agrup$F1_cat <- ifelse(left_cr13_gpcm_agrup$F1 < as.numeric(ponto_corte_gpcm_agrup[i]), "negativo","positivo")
          left_cr13_gpcm_agrup$F1_cat <- factor(left_cr13_gpcm_agrup$F1_cat, levels = c("negativo", "positivo"))
          ajuste_arvore_gpcm_agrup <- lm(left_cr13_gpcm_agrup$Q1~as.factor(left_cr13_gpcm_agrup$F1_cat) + var_comum_sign)
          valor_p_arvore_gpcm_agrup[i]<- anova(ajuste_arvore_gpcm_agrup)[[5]][1]
          
        } else {
          ponto_corte_gpcm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
          valor_p_arvore_gpcm_agrup[i] <- 'sem corte'
        }
      }
    }  
  }
  
  valor_p_agrup_grm
  valor_p_arvore_grm_agrup
  valor_p_agrup_gpcm
  valor_p_arvore_gpcm_agrup
  
  
  ############ ARRUMAR A PARTIR DAQUI
  alpha_grm<- rep((0.05/(length(seque)+qtd_sucessos_grm_cov)),length(genes))
  alpha_gpcm<- rep((0.05/(length(seque)+qtd_sucessos_gpcm_cov)),length(genes))
  
  
  #### FALTA ARRUMAR AQUI
  valor_p_agrup_grm_cov_exten<- NULL
  valor_p_agrup_gpcm_cov_exten<- NULL
  valor_p_arvore_grm_agrup_cov_exten<- NULL
  valor_p_arvore_gpcm_agrup_cov_exten<- NULL
  
  
  
  if (ncol(snps_agrup_cov)>20){
    alpha_grm<- rep(0.05/(length(seque) + qtd_sucessos_grm_cov), (ncol(snps_agrup_cov) - 2))
    alpha_gpcm<- rep(0.05/(length(seque) + qtd_sucessos_gpcm_cov), (ncol(snps_agrup_cov) - 2))
    for (i in 1:(length(seque)-1)){
      valor_p_agrup_grm_cov_exten<- c(valor_p_agrup_grm_cov_exten,rep(valor_p_agrup_grm[i],(seque[i+1]-seque[i])))
      valor_p_agrup_gpcm_cov_exten<- c(valor_p_agrup_gpcm_cov_exten,rep(valor_p_agrup_gpcm[i],(seque[i+1]-seque[i])))
      valor_p_arvore_grm_agrup_cov_exten<- c(valor_p_arvore_grm_agrup_cov_exten,rep(valor_p_arvore_grm_agrup[i],(seque[i+1]-seque[i])))
      valor_p_arvore_gpcm_agrup_cov_exten<- c(valor_p_arvore_gpcm_agrup_cov_exten,rep(valor_p_arvore_gpcm_agrup[i],(seque[i+1]-seque[i])))
    }
    i=i+1
    valor_p_agrup_grm_cov_exten<- c(valor_p_agrup_grm_cov_exten,rep(valor_p_agrup_grm[i],((ncol(snps_agrup_cov)-2)-(seque[i]-1))))
    valor_p_agrup_gpcm_cov_exten<- c(valor_p_agrup_gpcm_cov_exten,rep(valor_p_agrup_gpcm[i],((ncol(snps_agrup_cov)-2)-(seque[i]-1))))
    valor_p_arvore_grm_agrup_cov_exten<- c(valor_p_arvore_grm_agrup_cov_exten,rep(valor_p_arvore_grm_agrup[i],((ncol(snps_agrup_cov)-2)-(seque[i]-1))))
    valor_p_arvore_gpcm_agrup_cov_exten<- c(valor_p_arvore_gpcm_agrup_cov_exten,rep(valor_p_arvore_gpcm_agrup[i],((ncol(snps_agrup_cov)-2)-(seque[i]-1))))
  } else{
    alpha_grm<- rep(0.05/(length(seque) + qtd_sucessos_grm_cov), length(seque))
    alpha_gpcm<- rep(0.05/(length(seque) + qtd_sucessos_gpcm_cov), length(seque))
    valor_p_agrup_grm_cov_exten<- c(valor_p_agrup_grm_cov_exten,rep(valor_p_agrup_grm,length(agrup)))
    valor_p_agrup_gpcm_cov_exten<- c(valor_p_agrup_gpcm_cov_exten,rep(valor_p_agrup_gpcm,length(agrup)))
    valor_p_arvore_grm_agrup_cov_exten<- c(valor_p_arvore_grm_agrup_cov_exten,rep(valor_p_arvore_grm_agrup,length(agrup)))
    valor_p_arvore_gpcm_agrup_cov_exten<- c(valor_p_arvore_gpcm_agrup_cov_exten,rep(valor_p_arvore_gpcm_agrup,length(agrup)))
  }
  
  sign_gene_cr13_agrup_cov<-cbind(agrup, valor_p_agrup_grm_cov_exten,valor_p_agrup_gpcm_cov_exten,valor_p_arvore_grm_agrup_cov_exten, valor_p_arvore_gpcm_agrup_cov_exten, alpha_grm, alpha_gpcm, cr)
  colnames(sign_gene_cr13_agrup_cov)<- c("Genes AGRUP", "Valor-p GRM AGRUP COV","Valor-p GPCM AGRUP COV", "Valor-p Arvore GRM AGRUP COV", "Valor-p Arvore GPCM AGRUP COV", "Nivel de signi grm", "Nível de signi gpcm", "Cromossomo")
  
  
  #View(sign_gene_cr13_agrup_cov)
  
  ##### UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM GENE RARO É SIGNIFICATIVO
  ind_gene_raro_sign_bf_com_cov_agrup<- which(as.numeric(sign_gene_cr13_agrup_cov[,2])<(0.05/(length(seque) + qtd_sucessos_grm_cov))|
                                                as.numeric(sign_gene_cr13_agrup_cov[,3])<(0.05/(length(seque) + qtd_sucessos_gpcm_cov))|
                                                as.numeric(sign_gene_cr13_agrup_cov[,4])<(0.05/(length(seque) + qtd_sucessos_grm_cov))|
                                                as.numeric(sign_gene_cr13_agrup_cov[,5])<(0.05/(length(seque) + qtd_sucessos_gpcm_cov)))
  
  gene_raro_sign_bf_com_cov_agrup<- as.matrix(sign_gene_cr13_agrup_cov[ind_gene_raro_sign_bf_com_cov_agrup,])
  #View(gene_raro_sign_bf_com_cov_agrup)
  
  #### REMOVER DUPLICAÇÃO POR GENE E VALOR-P IGUAL
  gene_raro_sign_bf_com_cov_agrup_final <- gene_raro_sign_bf_com_cov_agrup[!duplicated(gene_raro_sign_bf_com_cov_agrup), ]
  
  
  #View(gene_raro_sign_bf_com_cov_agrup_final)
  
  
  
  
  
  
  
  
  ################ EXPORTANDO
  
  if (dim(as.matrix(gene_raro_sign_bf_com_cov_agrup_final))[2]==1){
    gene_raro_sign_bf_com_cov_agrup_final<- matrix(gene_raro_sign_bf_com_cov_agrup_final,nrow=1)
  }
  #cat('',gene_raro_sign_bf_com_cov_agrup_final,file=paste("gene_raro_com_cov_agrup_repli",repli,".txt",sep=""),append=T)
  write.table(
    gene_raro_sign_bf_com_cov_agrup_final,
    file = paste0("gene_raro_com_cov_agrup_dados_ISA.txt"),
    sep = ",",
    row.names = FALSE,
    col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
    append = (cr != 1)
  )
  
  
  
  
  #### AGORA AGRUPADO
  tabelas_a_agrup = list(summary(a_grm_cov_agrup),
                         summary(a_gpcm_cov_agrup))
  
  
  
  nome_tabelas_a_agrup = c("Summary a grm com cov agrup",
                           "Summary a gpcm com cov agrup")
  
  cat('Cromossomo',cr,"\n",file=paste("summary_coef_a_agrup_dados_ISA.txt",sep=""),append=T)
  for (i in seq_along(tabelas_a_agrup)) {
    cat(nome_tabelas_a_agrup[i], "\n", file=paste("summary_coef_a_agrup_dados_ISA.txt",sep=""),append=T)
    capture.output(print(tabelas_a_agrup[[i]], row.names = FALSE), file=paste("summary_coef_a_agrup_dados_ISA.txt",sep=""),append=T)
    cat("\n", file=paste("summary_coef_a_agrup_dados_ISA.txt",sep=""),append=T)
  }
  
  # Linha em branco extra entre cromossomos
  cat("\n", file=paste("summary_coef_a_agrup_dados_ISA.txt",sep=""),append=T)
  
  
  
  
  
  
  tabelas_b_agrup <- list(summary(b1_grm_cov_agrup),
                          summary(b1_gpcm_cov_agrup))
  nomes_tabelas_b_agrup <- c("Summary b1 grm com cov agrup",
                             "Summary b1 gpcm com cov agrup")
  
  cat('Cromossomo',cr,"\n",file=paste("summary_coef_b1_agrup_dados_ISA.txt",sep=""),append=T)
  for (i in seq_along(tabelas_b_agrup)) {
    cat(nomes_tabelas_b_agrup[i], "\n", file=paste("summary_coef_b1_agrup_dados_ISA.txt",sep=""),append=T)
    capture.output(print(tabelas_b_agrup[[i]], row.names = FALSE), file=paste("summary_coef_b1_agrup_dados_ISA.txt",sep=""),append=T)
    cat("\n", file=paste("summary_coef_b1_agrup_dados_ISA.txt",sep=""),append=T)
  }
  # Linha em branco extra entre cromossomos
  
  cat("\n", file=paste("summary_coef_b1_agrup_dados_ISA.txt",sep=""),append=T)
}
1+1
