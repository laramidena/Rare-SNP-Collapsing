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
#### LEITURA DA BASE DE DADOS DO CROMOSSOMO 1 E DO FENÓTIPO ####
for (repli in 1:10){
  for (cr in 1:22){
    snp_info <-  read_csv("~/Dados/snp_info")
    unr_phen <- read_csv(paste0("~/Dados/unr_phen.",repli))
    #View(unr_phen)
    fuma<- unr_phen$SMOKE
    idade<- unr_phen$AGE
    sexo<- unr_phen$SEX-1
    
    c13_snps <- read_csv(paste0("~/Dados/c", cr, "_snps.unr"))
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
    
    #### CLASSIFICACAO DAS BASE DE DADOS EM 0, 1 E 2 ####
    dados_13<-matrix(0,nrow=nrow(c13_snps),ncol = ncol(c13_snps))
    enableJIT(3)
    for (i in 1:nrow(dados_13)){
      for (j in 1: ncol(dados_13)){
        if (c13_snps[i,j]=="A/A"|c13_snps[i,j]=="T/T") dados_13[i,j]=2
        if (c13_snps[i,j]=="G/G"|c13_snps[i,j]=="C/C") dados_13[i,j]=0
        if (c13_snps[i,j]!="A/A" & c13_snps[i,j]!="T/T" &
            c13_snps[i,j]!="G/G" & c13_snps[i,j]!="C/C") dados_13[i,j]=1
      }
    }
    nomes_13 = colnames(c13_snps)
    colnames(dados_13)<- nomes_13
    
    
    freq_alelica_SNPs_cr13 <-matrix(0,nrow=2,ncol = ncol(c13_snps))
    row.names(freq_alelica_SNPs_cr13)<- c('freq alelo A','freq alelo a')
    colnames(freq_alelica_SNPs_cr13)<- colnames(c13_snps)
    
    for (i in 1:ncol(dados_13)){
      freq_alelica_SNPs_cr13[1,i]<- (0*length(dados_13[which(dados_13[,i]=='2'),i])/(2*nrow(dados_13))) + 
        (1*length(dados_13[which(dados_13[,i]=='1'),i])/(2*nrow(dados_13))) +
        (2*length(dados_13[which(dados_13[,i]=='0'),i])/(2*nrow(dados_13)))
      
      freq_alelica_SNPs_cr13[2,i]<- (2*length(dados_13[which(dados_13[,i]=='2'),i])/(2*nrow(dados_13))) + 
        (1*length(dados_13[which(dados_13[,i]=='1'),i])/(2*nrow(dados_13))) +
        (0*length(dados_13[which(dados_13[,i]=='0'),i])/(2*nrow(dados_13)))
      
    }
    
    #separando os SNPs em raros (MAF<5%) e comuns e colocando os nomes deles
    
    indice_SNPs_raros_cr_13<-which(freq_alelica_SNPs_cr13[,1:ncol(dados_13)]<0.05)
    for (i in 1:length(indice_SNPs_raros_cr_13)){
      if (indice_SNPs_raros_cr_13[i]%%2==0) indice_SNPs_raros_cr_13[i] = indice_SNPs_raros_cr_13[i]/2
      else indice_SNPs_raros_cr_13[i] = (indice_SNPs_raros_cr_13[i]+1)/2
    }
    
    SNPs_raros_cr13<- dados_13[,indice_SNPs_raros_cr_13]
    colnames(SNPs_raros_cr13)<- nomes_13[indice_SNPs_raros_cr_13]
    SNPs_comuns_cr13<- dados_13[,-c(indice_SNPs_raros_cr_13)]
    
    ### FAZENDO O MAF DE CADA SNP DO CR 13 E VERIFICANDO COM A BASE ORIGINAL
    maf_13<- NULL
    for (i in 1:ncol(freq_alelica_SNPs_cr13)){
      maf_13[i]<- min(freq_alelica_SNPs_cr13[1,i],freq_alelica_SNPs_cr13[2,i])
    }
    
    
    ### SEPARANDO OS SNPS RAROS DE ACORDO COM O GENE
    gene<- matrix(0,2,ncol(dados_13))
    tab_nomes<-as.matrix(colnames(SNPs_raros_cr13))
    colnames(tab_nomes)<-c("snp_name")
    tab_nomes<-as.data.frame(tab_nomes)
    SNPs_raros_cr13<-as.data.frame(SNPs_raros_cr13)
    tabela_full_cr13<- tab_nomes %>% left_join(snp_info,by='snp_name')
    tabela_gene_cr13<- tabela_full_cr13[,-c(2,3,5,6,7)]
    
    #### FENÓTIPO
    ytr_1 = NULL #cria uma lista vazia
    for (i in 1:nrow(unr_phen)){
      ytr_1[i] = unr_phen$Q1[i]
    }
    
    
    #### TENTANDO AGRUPAR OS SNPS Q ESTÃO NO MESMO GENE
    genes<-as.character(data.frame(table(tabela_gene_cr13[,2]))[,1])
    
    
    ### TENTANDO FAZER A SELEÇÃO DOS SNPS RAROS DE TODO O CR13
    corr_ML_grm <- NULL
    corr_ML_gpcm<- NULL
    valor_p_grm<- NULL
    valor_p_gpcm<- NULL
    matriz_grm<- matrix(0,nrow=697,ncol=length(genes))
    matriz_gpcm<- matrix(0,nrow=697,ncol=length(genes))
    ponto_corte_grm<- NULL
    ponto_corte_gpcm<- NULL
    valor_p_arvore_grm<- NULL
    valor_p_arvore_gpcm<- NULL
    a_grm<- NULL
    b1_grm<- NULL
    a_gpcm<- NULL
    b1_gpcm<- NULL
    
    
    library(mirt)
    set.seed(2511)
    enableJIT(3)
    for (i in 1:length(genes)){
      snps<-which(tabela_gene_cr13[,2]==genes[i]) #aqui sai a posição na lista genes (tá em ordem alfabética) dos snps q estão no gene pedido da tabela_gene_cr13
      snps_val <- SNPs_raros_cr13[,snps]
      nomes_col_gene <- colnames(SNPs_raros_cr13[which(tabela_gene_cr13[,2]==genes[i])])
      resultado<- try({
        mod_grm <- mirt(snps_val, 1, itemtype = "graded", TOL = 0.0001,
                        technical = list(NCYCLES = 100000))
        mod_gpcm <- mirt(snps_val, 1, itemtype = "gpcm", TOL = 0.0001,
                         technical = list(NCYCLES = 100000))
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
      
      if ((cr==22) & (i==50)) {
        valor_p_grm[i] <- 'problema'
        valor_p_gpcm[i] <- 'problema'
        ponto_corte_grm[i]<- 'problema'
        ponto_corte_gpcm[i]<- 'problema'
        valor_p_arvore_grm[i]<- 'problema'
        valor_p_arvore_gpcm[i]<- 'problema'
        next
      } 
      else {
        a_grm<- c(a_grm,coef(mod_grm, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_grm<-c(b1_grm,coef(mod_grm, IRTpars = TRUE, simplify = TRUE)$items[,2])
        fs_grm<-fscores(mod_grm,full.scores=FALSE, method = 'ML')
        resp_mirt_df_grm<- as.data.frame(fs_grm)
        resp_mirt_df_grm$ID <- do.call(paste0, resp_mirt_df_grm[ ,nomes_col_gene])
        snps_val$ID <- do.call(paste0, snps_val[ ,nomes_col_gene])
        snps_val$Q1 <- unr_phen$Q1
        left_cr13_grm <- dplyr::left_join(snps_val, dplyr::select(resp_mirt_df_grm, c(ID, F1)), by = "ID")
        
        
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
        
        ajuste_grm<- lm(left_cr13_grm$Q1~(left_cr13_grm$F1))
        
        matriz_grm[,i]<-left_cr13_grm$F1
        valor_p_grm[i]<- anova(ajuste_grm)[[5]][1]
        
        ### categórico grm
        
        arvore_reg_grm <- rpart(left_cr13_grm$Q1~(left_cr13_grm$F1),  method = "anova", minbucket=10,maxdepth=1)
        
        if (!is.null(arvore_reg_grm$splits)) {
          ponto_corte_grm[i] <- arvore_reg_grm$splits[,"index"][1]
          left_cr13_grm$F1_cat <- ifelse(left_cr13_grm$F1 < as.numeric(ponto_corte_grm[i]), "negativo","positivo")
          left_cr13_grm$F1_cat <- factor(left_cr13_grm$F1_cat, levels = c("negativo", "positivo"))
          ajuste_arvore_grm <- lm(left_cr13_grm$Q1~(left_cr13_grm$F1_cat))
          valor_p_arvore_grm[i]<- anova(ajuste_arvore_grm)[[5]][1]
          
        } else {
          ponto_corte_grm[i] <- 'sem corte'  # ou qualquer outro valor padrão
          valor_p_arvore_grm[i] <- 'sem corte'
        }
        
        #### AGORA É GPCM
        a_gpcm<- c(a_gpcm,coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_gpcm<-c(b1_gpcm,coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)$items[,2])
        fs_gpcm<-fscores(mod_gpcm,full.scores=FALSE, method = 'ML')
        resp_mirt_df_gpcm<- as.data.frame(fs_gpcm)
        resp_mirt_df_gpcm$ID <- do.call(paste0, resp_mirt_df_gpcm[ ,nomes_col_gene])
        snps_val$ID <- do.call(paste0, snps_val[ ,nomes_col_gene])
        snps_val$Q1 <- unr_phen$Q1
        left_cr13_gpcm <- dplyr::left_join(snps_val, dplyr::select(resp_mirt_df_gpcm, c(ID, F1)), by = "ID")
        
        
        ## substitui -inf no F1 pelo [menor valor (sem ser o -inf) -1] e
        # +inf no F1 por [maior valor (sem ser +inf) +1]
        min_finito_gpcm <- min(left_cr13_gpcm$F1[is.finite(left_cr13_gpcm$F1)])-1
        max_finito_gpcm <- max(left_cr13_gpcm$F1[is.finite(left_cr13_gpcm$F1)])+1
        left_cr13_gpcm$F1[left_cr13_gpcm$F1 == -Inf] <- min_finito_gpcm
        left_cr13_gpcm$F1[left_cr13_gpcm$F1 == Inf] <- max_finito_gpcm
        
        
        corr_ML_gpcm[i]<- cor(left_cr13_gpcm$Q1,left_cr13_gpcm$F1)
        ajuste_gpcm<- lm(left_cr13_gpcm$Q1~(left_cr13_gpcm$F1))
        
        matriz_gpcm[,i]<-left_cr13_gpcm$F1
        valor_p_gpcm[i]<- anova(ajuste_gpcm)[[5]][1]
        
        ### categórico grm
        
        arvore_reg_gpcm <- rpart(left_cr13_gpcm$Q1~(left_cr13_gpcm$F1),  method = "anova", minbucket=10,maxdepth=1)
        
        
        if (!is.null(arvore_reg_gpcm$splits)) {
          ponto_corte_gpcm[i] <- arvore_reg_gpcm$splits[,"index"][1]
          left_cr13_gpcm$F1_cat <- ifelse(left_cr13_gpcm$F1 < as.numeric(ponto_corte_gpcm[i]), "negativo","positivo")
          left_cr13_gpcm$F1_cat <- factor(left_cr13_gpcm$F1_cat, levels = c("negativo", "positivo"))
          ajuste_arvore_gpcm <- lm(left_cr13_gpcm$Q1~as.factor(left_cr13_gpcm$F1_cat))
          valor_p_arvore_gpcm[i]<- anova(ajuste_arvore_gpcm)[[5]][1]
          
        } else {
          ponto_corte_gpcm[i] <- 'sem corte'  # ou qualquer outro valor padrão
          valor_p_arvore_gpcm[i] <- 'sem corte'
        }
        
        
      }
    }
    valor_p_grm
    valor_p_gpcm
    valor_p_arvore_grm
    valor_p_arvore_gpcm
    
    
    
    qtd_sucesso_grm<- length(which(is.na(valor_p_grm)==F))
    qtd_sucesso_gpcm<- length(which(is.na(valor_p_gpcm)==F))
    
    ########## PEGANDO OS GENES QUE TEM POUCA VARIABILIDADE E QUE NÃO RODOU O GRM E O GPCM, AGRUPANDO ELES
    ########## E RODANDO O GRM E O GPCM PRA ELES AGRUPADOS
    sign_gene_cr13<-cbind(genes, valor_p_grm,valor_p_gpcm, valor_p_arvore_grm,valor_p_arvore_gpcm, cr)
    colnames(sign_gene_cr13)<- c("Genes", "Valor-p GRM ML","Valor-p GPCM ML", "Valor-p Arvore GRM", "Valor-p Arvore GPCM", "Cromossomo")
    
    
    sign_gene_cr13_df <- data.frame(sign_gene_cr13)
    pos_snps_pouco_var <- which((is.na(sign_gene_cr13_df$Valor.p.GRM.ML)==TRUE) | (is.na(sign_gene_cr13_df$Valor.p.GPCM.ML)==TRUE) | (sign_gene_cr13_df$Valor.p.GRM.ML=='problema') | (sign_gene_cr13_df$Valor.p.GPCM.ML=='problema'))
    genes_na <- as.data.frame(sign_gene_cr13_df[pos_snps_pouco_var,1])
    colnames(genes_na)<-c("gene")
    
    snps_filtrados <- tabela_gene_cr13[tabela_gene_cr13$gene %in% genes_na$gene, ]
    ### PEGANDO SÓ OS NOMES DOS SNPS
    snps_dos_genes_na <- snps_filtrados$snp
    
    ### PEGANDO O GENÓTIPO DESSES SNPS (PEGA O GENÓTIPO DE TODOS OS Q NÃO RODARAM)
    snps_agrup <- SNPs_raros_cr13[, snps_dos_genes_na, drop = FALSE]
    nomes_agrup<- colnames(snps_agrup)
    #View(snps_agrup)
    valor_p_agrup_grm<- NULL
    valor_p_arvore_grm_agrup<- NULL
    ponto_corte_grm_agrup<- NULL
    valor_p_agrup_gpcm<- NULL
    valor_p_arvore_gpcm_agrup<- NULL
    ponto_corte_gpcm_agrup<- NULL
    a_grm_agrup<- NULL
    b1_grm_agrup<- NULL
    a_gpcm_agrup<- NULL
    b1_gpcm_agrup<- NULL
    agrup<- NULL
    enableJIT(3)
    if (ncol(snps_agrup)>20){
      seque<- seq(1,ncol(snps_agrup), by = 10)
      for(i in 1:(length(seque)-1)){
        snps_agrupados<- snps_agrup[,seque[i]:(seque[i+1]-1)]
        mod_grm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "graded", TOL = 0.0001,
                                  technical = list(NCYCLES = 100000)),silent=TRUE)
        mod_gpcm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "gpcm", TOL = 0.0001,
                                   technical = list(NCYCLES = 100000)),silent=TRUE)
        nomes_agrup<- colnames(snps_agrupados)
        agrup<- c(agrup,snps_filtrados[seque[i]:(seque[i+1]-1),2])
        
        if ((cr==11) & (i==17)) {
          valor_p_agrup_grm[i]<- 'problema'
          valor_p_arvore_grm_agrup[i]<- 'problema'
          ponto_corte_grm_agrup[i]<- 'problema'
          a_grm_agrup[i]<- 'problema'
          b1_grm_agrup[i]<- 'problema'
          next
        }
        
        if (inherits(mod_grm_agrup, "try-error")) {
          valor_p_agrup_grm[i]<- 'funcao falhou'
          valor_p_arvore_grm_agrup[i]<- 'funcao falhou'
          ponto_corte_grm_agrup[i]<- 'funcao falhou'
          a_grm_agrup[i]<- 'funcao falhou'
          b1_grm_agrup[i]<- 'funcao falhou'
          
        }
        else{
          a_grm_agrup<- c(a_grm_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
          b1_grm_agrup<-c(b1_grm_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
          fs_grm_agrup<-fscores(mod_grm_agrup,full.scores=FALSE, method = 'ML')
          resp_mirt_df_grm_agrup<- as.data.frame(fs_grm_agrup)
          resp_mirt_df_grm_agrup$ID <- do.call(paste0, resp_mirt_df_grm_agrup[ ,nomes_agrup])
          snps_agrup$ID <- do.call(paste0, snps_agrup[ ,nomes_agrup])
          snps_agrup$Q1 <- unr_phen$Q1
          
          left_cr13_grm_agrup <- dplyr::left_join(snps_agrup, dplyr::select(resp_mirt_df_grm_agrup, c(ID, F1)), by = "ID")
          
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
            ajuste_grm_agrup<- lm(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1))
            
            valor_p_agrup_grm[i]<- anova(ajuste_grm_agrup)[[5]][1]
            
            
            arvore_reg_grm_agrup<- rpart(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
            
            if (!is.null(arvore_reg_grm_agrup$splits)) {
              ponto_corte_grm_agrup[i] <- arvore_reg_grm_agrup$splits[,"index"][1]
              left_cr13_grm_agrup$F1_cat <- ifelse(left_cr13_grm_agrup$F1 < as.numeric(ponto_corte_grm_agrup[i]), "negativo","positivo")
              left_cr13_grm_agrup$F1_cat <- factor(left_cr13_grm_agrup$F1_cat, levels = c("negativo", "positivo"))
              ajuste_arvore_grm_agrup <- lm(left_cr13_grm_agrup$Q1~as.factor(left_cr13_grm_agrup$F1_cat))
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
          a_gpcm_agrup[i]<- 'funcao falhou'
          b1_gpcm_agrup[i]<- 'funcao falhou'
          
        } else {
          a_gpcm_agrup<- c(a_gpcm_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
          b1_gpcm_agrup<-c(b1_gpcm_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
          
          
          #### AGORA É GPCM
          fs_gpcm_agrup<-fscores(mod_gpcm_agrup,full.scores=FALSE, method = 'ML')
          resp_mirt_df_gpcm_agrup<- as.data.frame(fs_gpcm_agrup)
          resp_mirt_df_gpcm_agrup$ID <- do.call(paste0, resp_mirt_df_gpcm_agrup[ ,nomes_agrup])
          snps_agrup$ID <- do.call(paste0, snps_agrup[ ,nomes_agrup])
          snps_agrup$Q1 <- unr_phen$Q1
          left_cr13_gpcm_agrup <- dplyr::left_join(snps_agrup, dplyr::select(resp_mirt_df_gpcm_agrup, c(ID, F1)), by = "ID")
          
          
          
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
            ajuste_gpcm_agrup<- lm(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1))
            valor_p_agrup_gpcm[i]<- anova(ajuste_gpcm_agrup)[[5]][1]
            
            
            arvore_reg_gpcm_agrup <- rpart(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
            
            if (!is.null(arvore_reg_gpcm_agrup$splits)) {
              ponto_corte_gpcm_agrup[i] <- arvore_reg_gpcm_agrup$splits[,"index"][1]
              left_cr13_gpcm_agrup$F1_cat <- ifelse(left_cr13_gpcm_agrup$F1 < as.numeric(ponto_corte_gpcm_agrup[i]), "negativo","positivo")
              left_cr13_gpcm_agrup$F1_cat <- factor(left_cr13_gpcm_agrup$F1_cat, levels = c("negativo", "positivo"))
              ajuste_arvore_gpcm_agrup <- lm(left_cr13_gpcm_agrup$Q1~as.factor(left_cr13_gpcm_agrup$F1_cat))
              valor_p_arvore_gpcm_agrup[i]<- anova(ajuste_arvore_gpcm_agrup)[[5]][1]
              
            } else {
              ponto_corte_gpcm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
              valor_p_arvore_gpcm_agrup[i] <- 'sem corte'
            }
          } 
        }  
      }
      i=i+1
      snps_agrupados<- snps_agrup[,seque[i]:(ncol(snps_agrup)-2)]
      mod_grm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "graded", TOL = 0.0001,
                                technical = list(NCYCLES = 100000)),silent=TRUE)
      mod_gpcm_agrup <- try(mirt(snps_agrupados, 1, itemtype = "gpcm", TOL = 0.0001,
                                 technical = list(NCYCLES = 100000)),silent=TRUE)
      nomes_agrup<- colnames(snps_agrupados)
      #agrup<- snps_filtrados[seque[i]:seque[i+1]-1,2]
      agrup<- c(agrup,snps_filtrados[seque[i]:(ncol(snps_agrup)-2),2])
      
      
      if (inherits(mod_grm_agrup, "try-error")) {
        valor_p_agrup_grm[i]<- 'funcao falhou'
        valor_p_arvore_grm_agrup[i]<- 'funcao falhou'
        ponto_corte_grm_agrup[i]<- 'funcao falhou'
        a_grm_agrup<- 'funcao falhou'
        b1_grm_agrup<- 'funcao falhou'
      } else{
        a_grm_agrup<- c(a_grm_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_grm_agrup<-c(b1_grm_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
        fs_grm_agrup<-fscores(mod_grm_agrup,full.scores=FALSE, method = 'ML')
        resp_mirt_df_grm_agrup<- as.data.frame(fs_grm_agrup)
        resp_mirt_df_grm_agrup$ID <- do.call(paste0, resp_mirt_df_grm_agrup[ ,nomes_agrup])
        snps_agrup$ID <- do.call(paste0, snps_agrup[ ,nomes_agrup])
        snps_agrup$Q1 <- unr_phen$Q1
        
        left_cr13_grm_agrup <- dplyr::left_join(snps_agrup, dplyr::select(resp_mirt_df_grm_agrup, c(ID, F1)), by = "ID")
        
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
          ajuste_grm_agrup<- lm(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1))
          
          valor_p_agrup_grm[i]<- anova(ajuste_grm_agrup)[[5]][1]
          
          
          arvore_reg_grm_agrup<- rpart(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
          
          if (!is.null(arvore_reg_grm_agrup$splits)) {
            ponto_corte_grm_agrup[i] <- arvore_reg_grm_agrup$splits[,"index"][1]
            left_cr13_grm_agrup$F1_cat <- ifelse(left_cr13_grm_agrup$F1 < as.numeric(ponto_corte_grm_agrup[i]), "negativo","positivo")
            left_cr13_grm_agrup$F1_cat <- factor(left_cr13_grm_agrup$F1_cat, levels = c("negativo", "positivo"))
            ajuste_arvore_grm_agrup <- lm(left_cr13_grm_agrup$Q1~as.factor(left_cr13_grm_agrup$F1_cat))
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
        a_gpcm_agrup<- 'funcao falhou'
        b1_gpcm_agrup<- 'funcao falhou'
        
      } else {
        a_gpcm_agrup<- c(a_gpcm_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_gpcm_agrup<-c(b1_gpcm_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
        
        
        #### AGORA É GPCM
        fs_gpcm_agrup<-fscores(mod_gpcm_agrup,full.scores=FALSE, method = 'ML')
        resp_mirt_df_gpcm_agrup<- as.data.frame(fs_gpcm_agrup)
        resp_mirt_df_gpcm_agrup$ID <- do.call(paste0, resp_mirt_df_gpcm_agrup[ ,nomes_agrup])
        snps_agrup$ID <- do.call(paste0, snps_agrup[ ,nomes_agrup])
        snps_agrup$Q1 <- unr_phen$Q1
        left_cr13_gpcm_agrup <- dplyr::left_join(snps_agrup, dplyr::select(resp_mirt_df_gpcm_agrup, c(ID, F1)), by = "ID")
        
        
        
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
          ajuste_gpcm_agrup<- lm(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1))
          valor_p_agrup_gpcm[i]<- anova(ajuste_gpcm_agrup)[[5]][1]
          
          
          arvore_reg_gpcm_agrup <- rpart(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
          
          if (!is.null(arvore_reg_gpcm_agrup$splits)) {
            ponto_corte_gpcm_agrup[i] <- arvore_reg_gpcm_agrup$splits[,"index"][1]
            left_cr13_gpcm_agrup$F1_cat <- ifelse(left_cr13_gpcm_agrup$F1 < as.numeric(ponto_corte_gpcm_agrup[i]), "negativo","positivo")
            left_cr13_gpcm_agrup$F1_cat <- factor(left_cr13_gpcm_agrup$F1_cat, levels = c("negativo", "positivo"))
            ajuste_arvore_gpcm_agrup <- lm(left_cr13_gpcm_agrup$Q1~as.factor(left_cr13_gpcm_agrup$F1_cat))
            valor_p_arvore_gpcm_agrup[i]<- anova(ajuste_arvore_gpcm_agrup)[[5]][1]
            
          } else {
            ponto_corte_gpcm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
            valor_p_arvore_gpcm_agrup[i] <- 'sem corte'
          }
        } 
      }
    }
    
    valor_p_agrup_grm
    valor_p_agrup_gpcm
    valor_p_arvore_gpcm_agrup
    valor_p_arvore_grm_agrup
    
    ###############
    if (ncol(snps_agrup)<=20){
      i=1
      seque=1
      set.seed(2511)
      mod_grm_agrup <- try(mirt(snps_agrup, 1, itemtype = "graded", TOL = 0.0001,
                                technical = list(NCYCLES = 100000)),silent=TRUE)
      mod_gpcm_agrup <- try(mirt(snps_agrup, 1, itemtype = "gpcm", TOL = 0.0001,
                                 technical = list(NCYCLES = 100000)),silent=TRUE)
      nomes_agrup<- colnames(snps_agrup)
      agrup<- c(agrup,snps_filtrados[,2])
      
      if (inherits(mod_grm_agrup, "try-error")) {
        valor_p_agrup_grm[i]<- 'funcao falhou'
        valor_p_arvore_grm_agrup[i]<- 'funcao falhou'
        ponto_corte_grm_agrup[i]<- 'funcao falhou'
        a_grm_agrup<- 'funcao falhou'
        b1_grm_agrup<- 'funcao falhou'
        
      }
      else{
        a_grm_agrup<- c(a_grm_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_grm_agrup<-c(b1_grm_agrup,coef(mod_grm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
        fs_grm_agrup<-fscores(mod_grm_agrup,full.scores=FALSE, method = 'ML')
        resp_mirt_df_grm_agrup<- as.data.frame(fs_grm_agrup)
        resp_mirt_df_grm_agrup$ID <- do.call(paste0, resp_mirt_df_grm_agrup[ ,nomes_agrup])
        snps_agrup$ID <- do.call(paste0, snps_agrup[ ,nomes_agrup])
        snps_agrup$Q1 <- unr_phen$Q1
        
        left_cr13_grm_agrup <- dplyr::left_join(snps_agrup, dplyr::select(resp_mirt_df_grm_agrup, c(ID, F1)), by = "ID")
        
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
          ajuste_grm_agrup<- lm(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1))
          
          valor_p_agrup_grm[i]<- anova(ajuste_grm_agrup)[[5]][1]
          
          
          arvore_reg_grm_agrup<- rpart(left_cr13_grm_agrup$Q1~(left_cr13_grm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
          
          if (!is.null(arvore_reg_grm_agrup$splits)) {
            ponto_corte_grm_agrup[i] <- arvore_reg_grm_agrup$splits[,"index"][1]
            left_cr13_grm_agrup$F1_cat <- ifelse(left_cr13_grm_agrup$F1 < as.numeric(ponto_corte_grm_agrup[i]), "negativo","positivo")
            left_cr13_grm_agrup$F1_cat <- factor(left_cr13_grm_agrup$F1_cat, levels = c("negativo", "positivo"))
            ajuste_arvore_grm_agrup <- lm(left_cr13_grm_agrup$Q1~as.factor(left_cr13_grm_agrup$F1_cat))
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
        a_gpcm_agrup<- 'funcao falhou'
        b1_gpcm_agrup<- 'funcao falhou'
        
      } else {
        a_gpcm_agrup<- c(a_gpcm_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_gpcm_agrup<-c(b1_gpcm_agrup,coef(mod_gpcm_agrup, IRTpars = TRUE, simplify = TRUE)$items[,2])
        
        
        #### AGORA É GPCM
        fs_gpcm_agrup<-fscores(mod_gpcm_agrup,full.scores=FALSE, method = 'ML')
        resp_mirt_df_gpcm_agrup<- as.data.frame(fs_gpcm_agrup)
        resp_mirt_df_gpcm_agrup$ID <- do.call(paste0, resp_mirt_df_gpcm_agrup[ ,nomes_agrup])
        snps_agrup$ID <- do.call(paste0, snps_agrup[ ,nomes_agrup])
        snps_agrup$Q1 <- unr_phen$Q1
        left_cr13_gpcm_agrup <- dplyr::left_join(snps_agrup, dplyr::select(resp_mirt_df_gpcm_agrup, c(ID, F1)), by = "ID")
        
        
        
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
          ajuste_gpcm_agrup<- lm(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1))
          valor_p_agrup_gpcm[i]<- anova(ajuste_gpcm_agrup)[[5]][1]
          
          
          arvore_reg_gpcm_agrup <- rpart(left_cr13_gpcm_agrup$Q1~(left_cr13_gpcm_agrup$F1),  method = "anova", minbucket=10,maxdepth=1)
          
          if (!is.null(arvore_reg_gpcm_agrup$splits)) {
            ponto_corte_gpcm_agrup[i] <- arvore_reg_gpcm_agrup$splits[,"index"][1]
            left_cr13_gpcm_agrup$F1_cat <- ifelse(left_cr13_gpcm_agrup$F1 < as.numeric(ponto_corte_gpcm_agrup[i]), "negativo","positivo")
            left_cr13_gpcm_agrup$F1_cat <- factor(left_cr13_gpcm_agrup$F1_cat, levels = c("negativo", "positivo"))
            ajuste_arvore_gpcm_agrup <- lm(left_cr13_gpcm_agrup$Q1~as.factor(left_cr13_gpcm_agrup$F1_cat))
            valor_p_arvore_gpcm_agrup[i]<- anova(ajuste_arvore_gpcm_agrup)[[5]][1]
            
          } else {
            ponto_corte_gpcm_agrup[i] <- 'sem corte'  # ou qualquer outro valor padrão
            valor_p_arvore_gpcm_agrup[i] <- 'sem corte'
          }
        }
      }  
    }
    valor_p_arvore_grm_agrup
    valor_p_arvore_gpcm_agrup
    
    
    
    alpha_grm<- rep((0.05/(length(seque)+qtd_sucesso_grm)), length(genes))
    alpha_gpcm<- rep((0.05/(length(seque)+qtd_sucesso_gpcm)), length(genes))
    
    sign_gene_cr13<-cbind(genes, valor_p_grm,valor_p_gpcm, valor_p_arvore_grm,valor_p_arvore_gpcm,alpha_grm, alpha_gpcm, cr)
    colnames(sign_gene_cr13)<- c("Genes", "Valor-p GRM ML","Valor-p GPCM ML", "Valor-p Arvore GRM", "Valor-p Arvore GPCM", "Nível sign grm", "Nível sign gpcm", "Cromossomo")
    #View(sign_gene_cr13)
    
    ##### UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM GENE RARO É SIGNIFICATIVO
    ind_gene_raro_sign_bf_sem_cov<- which(as.numeric(sign_gene_cr13[,2])< (0.05/(length(seque) + qtd_sucesso_grm))|
                                            as.numeric(sign_gene_cr13[,3])<(0.05/(length(seque) + qtd_sucesso_gpcm))|
                                            as.numeric(sign_gene_cr13[,4])<(0.05/(length(seque) + qtd_sucesso_grm))|
                                            as.numeric(sign_gene_cr13[,5])<(0.05/(length(seque) + qtd_sucesso_gpcm)))
    
    gene_raro_sign_bf_sem_cov<- as.matrix(sign_gene_cr13[ind_gene_raro_sign_bf_sem_cov,])
    #View(gene_raro_sign_bf_sem_cov)
    #View(sign_gene_cr13)
    
    
    
    valor_p_agrup_grm_exten<- NULL
    valor_p_agrup_gpcm_exten<- NULL
    valor_p_arvore_grm_agrup_exten<- NULL
    valor_p_arvore_gpcm_agrup_exten<- NULL
    
    if (ncol(snps_agrup)>20){
      alpha_grm<- rep(0.05/(length(seque) + qtd_sucesso_grm), (ncol(snps_agrup)-2))
      alpha_gpcm<- rep(0.05/(length(seque) + qtd_sucesso_grm), (ncol(snps_agrup)-2))
      for (i in 1:(length(seque)-1)){
        valor_p_agrup_grm_exten<- c(valor_p_agrup_grm_exten,rep(valor_p_agrup_grm[i],(seque[i+1]-seque[i])))
        valor_p_agrup_gpcm_exten<- c(valor_p_agrup_gpcm_exten,rep(valor_p_agrup_gpcm[i],(seque[i+1]-seque[i])))
        valor_p_arvore_grm_agrup_exten<- c(valor_p_arvore_grm_agrup_exten,rep(valor_p_arvore_grm_agrup[i],(seque[i+1]-seque[i])))
        valor_p_arvore_gpcm_agrup_exten<- c(valor_p_arvore_gpcm_agrup_exten,rep(valor_p_arvore_gpcm_agrup[i],(seque[i+1]-seque[i])))
      }
      i=i+1
      valor_p_agrup_grm_exten<- c(valor_p_agrup_grm_exten,rep(valor_p_agrup_grm[i],((ncol(snps_agrup)-2)-(seque[i]-1))))
      valor_p_agrup_gpcm_exten<- c(valor_p_agrup_gpcm_exten,rep(valor_p_agrup_gpcm[i],((ncol(snps_agrup)-2)-(seque[i]-1))))
      valor_p_arvore_grm_agrup_exten<- c(valor_p_arvore_grm_agrup_exten,rep(valor_p_arvore_grm_agrup[i],((ncol(snps_agrup)-2)-(seque[i]-1))))
      valor_p_arvore_gpcm_agrup_exten<- c(valor_p_arvore_gpcm_agrup_exten,rep(valor_p_arvore_gpcm_agrup[i],((ncol(snps_agrup)-2)-(seque[i]-1))))
    } else{
      alpha_grm<- rep(0.05/(length(seque) + qtd_sucesso_grm), length(seque))
      alpha_gpcm<- rep(0.05/(length(seque) + qtd_sucesso_grm), length(seque))
      valor_p_agrup_grm_exten<- c(valor_p_agrup_grm_exten,rep(valor_p_agrup_grm,length(agrup)))
      valor_p_agrup_gpcm_exten<- c(valor_p_agrup_gpcm_exten,rep(valor_p_agrup_gpcm,length(agrup)))
      valor_p_arvore_grm_agrup_exten<- c(valor_p_arvore_grm_agrup_exten,rep(valor_p_arvore_grm_agrup,length(agrup)))
      valor_p_arvore_gpcm_agrup_exten<- c(valor_p_arvore_gpcm_agrup_exten,rep(valor_p_arvore_gpcm_agrup,length(agrup)))
    }
    
    sign_gene_cr13_agrup<-cbind(agrup, valor_p_agrup_grm_exten,valor_p_agrup_gpcm_exten,valor_p_arvore_grm_agrup_exten, valor_p_arvore_gpcm_agrup_exten,alpha_grm, alpha_gpcm, cr)
    colnames(sign_gene_cr13_agrup)<- c("Genes AGRUP", "Valor-p GRM ML AGRUP","Valor-p GPCM ML AGRUP", "Valor-p GRM AGRUP Arvore", "Valor-p GPCM AGRUP Arvore","Nível sign GRM","Nível sign GPCM", "Cromossomo")
    
    
    #View(sign_gene_cr13_agrup)
    
    
    
    ##### UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM GENE RARO É SIGNIFICATIVO
    ind_gene_raro_sign_bf_sem_cov_agrup<- which(as.numeric(sign_gene_cr13_agrup[,2])<(0.05/(length(seque) + qtd_sucesso_grm))|
                                                  as.numeric(sign_gene_cr13_agrup[,3])<(0.05/(length(seque) + qtd_sucesso_grm))|
                                                  as.numeric(sign_gene_cr13_agrup[,4])<(0.05/(length(seque) + qtd_sucesso_grm))|
                                                  as.numeric(sign_gene_cr13_agrup[,5])<(0.05/(length(seque) + qtd_sucesso_grm)))
    
    gene_raro_sign_bf_sem_cov_agrup<- as.matrix(sign_gene_cr13_agrup[ind_gene_raro_sign_bf_sem_cov_agrup,])
    #View(gene_raro_sign_bf_sem_cov_agrup)
    
    #### REMOVER DUPLICAÇÃO POR GENE E VALOR-P IGUAL
    gene_raro_sign_bf_sem_cov_agrup_final <- gene_raro_sign_bf_sem_cov_agrup[!duplicated(gene_raro_sign_bf_sem_cov_agrup), ]
    
    #View(gene_raro_sign_bf_sem_cov_agrup_final)
    
    
    
    ######## AGORA É O MODELO COM COVARIÁVEIS SELECIONADAS PELO LASSO
    
    
    #### SELECIONANDO SNPS COMUNS PELO LASSO
    set.seed(2511)
    x<- as.matrix(SNPs_comuns_cr13)
    X<- cbind(x,fuma,idade,sexo)
    cv.out_1 <- cv.glmnet(as.matrix(X), ytr_1, alpha = 1, family="gaussian")
    plot(cv.out_1, ylab = "EQP", xlab = "log(lambda)")
    bestlam_1 <- cv.out_1$lambda.min
    bestlam_1 # o melhor valor de lambda é 0.04905032
    ### ARMAZENA O MELHOR LÂMBIDA E O GRAU DE LIBERDADE
    lasso.mod_1 <- glmnet(as.matrix(X), ytr_1, alpha = 1, lambda = bestlam_1,family="gaussian")
    lasso.coef_1  <- predict(lasso.mod_1, type = 'coefficients', s = bestlam_1)
    
    coef_lasso_cromo_1 = as.matrix(lasso.coef_1)
    ###REMOVE A PRIMEIRA LINHA DE "COEF_LASSO_CROMO_1" QUE ERA O VALOR DE "ALPHA"
    betas<- matrix(0,nrow=ncol(X),ncol=2)
    betas_cromo_1 = coef_lasso_cromo_1[c(-1)]
    betas[,1]<-colnames(X)
    betas[,2]<- betas_cromo_1
    snps_comuns_sign_cr13 = betas[which(betas[,2]!=0),]
    ### CONTÉM TODOS OS BETAS SIGN DE TODOS OS SNPs DO CROMOSSOMO 13
    snps_comuns_sign_cr13<- as.data.frame(snps_comuns_sign_cr13)
    colnames(snps_comuns_sign_cr13)<- c('Variável comum','Coeficiente beta')
    
    
    
    
    
    
    ##### JUNTAR SNPS RAROS E COMUNS SELECIONADOS PELAS METODOLOGIAS
    indice_var_comum_sign = which(betas[,2]!=0)
    var_comum_sign<- X[,indice_var_comum_sign]
    
    #### TESTANDO SIGNIFICANCIA DOS GENES DE SNPS RAROS USANDO COMO COVARIÁVEIS OS SNPS COMUNS E VARIÁVEIS SELECIONADAS PELO LASSO
    
    corr_ML_grm <- NULL
    corr_ML_gpcm<- NULL
    valor_p_grm<- NULL
    valor_p_gpcm<- NULL
    matriz_grm<- matrix(0,nrow=697,ncol=length(genes))
    matriz_gpcm<- matrix(0,nrow=697,ncol=length(genes))
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
    for (i in 1:length(genes)){
      snps<-which(tabela_gene_cr13[,2]==genes[i]) #aqui sai a posição na lista genes (tá em ordem alfabética) dos snps q estão no gene pedido da tabela_gene_cr13
      snps_val <- SNPs_raros_cr13[,snps]
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
      if ((cr==22) & (i==50)) {
        valor_p_grm[i] <- 'problema'
        valor_p_gpcm[i] <- 'problema'
        ponto_corte_grm[i]<- 'problema'
        ponto_corte_gpcm[i]<- 'problema'
        valor_p_arvore_grm[i]<- 'problema'
        valor_p_arvore_gpcm[i]<- 'problema'
        next
      } 
      else {
        a_grm_cov<-c(a_grm_cov,coef(mod_grm, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_grm_cov<-c(b1_grm_cov,coef(mod_grm, IRTpars = TRUE, simplify = TRUE)$items[,2])
        a_gpcm_cov<-c(a_gpcm_cov,coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)$items[,1])
        b1_gpcm_cov<-c(b1_gpcm_cov,coef(mod_gpcm, IRTpars = TRUE, simplify = TRUE)$items[,2])
        
        
        fs_grm<-fscores(mod_grm,full.scores=FALSE, method = 'ML')
        resp_mirt_df_grm<- as.data.frame(fs_grm)
        resp_mirt_df_grm$ID <- do.call(paste0, resp_mirt_df_grm[ ,nomes_col_gene])
        snps_val$ID <- do.call(paste0, snps_val[ ,nomes_col_gene])
        snps_val$Q1 <- unr_phen$Q1
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
        snps_val$Q1 <- unr_phen$Q1
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
    valor_p_arvore_grm
    valor_p_arvore_gpcm
    valor_p_grm
    valor_p_gpcm
    
    
    
    qtd_sucessos_grm_cov<-length(which(is.na(valor_p_grm)==FALSE))
    qtd_sucessos_gpcm_cov<- length(which(is.na(valor_p_gpcm)==FALSE))
    
    
    ###### PEGANDO OS GENES Q NÃO RODARAM ANTES, AGRUPANDO E RODANDO COM COVARIAVEIS
    
    
    sign_gene_cr13<-cbind(genes, valor_p_grm,valor_p_gpcm, valor_p_arvore_grm, valor_p_arvore_gpcm, cr)
    colnames(sign_gene_cr13)<- c("Genes", "Valor-p GRM C COV","Valor-p GPCM C COV", "Valor-p Arvore GRM C COV", "Valor-p Arvore GPCM C COV", "Cromossomo")
    
    sign_gene_cr13_cov_df <- data.frame(sign_gene_cr13)
    #View(sign_gene_cr13)
    pos_snps_pouco_var <- which((is.na(sign_gene_cr13_cov_df$Valor.p.GRM.C.COV)==TRUE) | (is.na(sign_gene_cr13_cov_df$Valor.p.GPCM.C.COV)==TRUE)| (sign_gene_cr13_cov_df$Valor.p.GRM.C.COV=='problema') | (sign_gene_cr13_cov_df$Valor.p.GPCM.C.COV=='problema'))
    #pos_snps_pouco_var ### genes com pouca variabilidade
    genes_na <- as.data.frame(sign_gene_cr13_df[pos_snps_pouco_var,1])
    colnames(genes_na)<-c("gene")
    
    snps_filtrados <- tabela_gene_cr13[tabela_gene_cr13$gene %in% genes_na$gene, ]
    ### PEGANDO SÓ OS NOMES DOS SNPS
    snps_dos_genes_na <- snps_filtrados$snp
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
        
        if ((cr==11) & (i==17)) {
          valor_p_agrup_grm[i]<- 'problema'
          valor_p_arvore_grm_agrup[i]<- 'problema'
          ponto_corte_grm_agrup[i]<- 'problema'
          a_grm_cov_agrup[i]<- 'problema'
          b1_grm_cov_agrup[i]<- 'problema'
          next
        }
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
          snps_agrup_cov$Q1 <- unr_phen$Q1
          
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
          snps_agrup_cov$Q1 <- unr_phen$Q1
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
        snps_agrup_cov$Q1 <- unr_phen$Q1
        
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
        snps_agrup_cov$Q1 <- unr_phen$Q1
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
        snps_agrup_cov$Q1 <- unr_phen$Q1
        
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
        snps_agrup_cov$Q1 <- unr_phen$Q1
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
    
    alpha_grm<- rep((0.05/(length(seque)+qtd_sucessos_grm_cov)),length(genes))
    alpha_gpcm<- rep((0.05/(length(seque)+qtd_sucessos_gpcm_cov)),length(genes))
    
    
    sign_gene_cr13<-cbind(genes, valor_p_grm,valor_p_gpcm, valor_p_arvore_grm, valor_p_arvore_gpcm, alpha_grm, alpha_gpcm, cr)
    colnames(sign_gene_cr13)<- c("Genes", "Valor-p GRM C COV","Valor-p GPCM C COV", "Valor-p Arvore GRM C COV", "Valor-p Arvore GPCM C COV", "Nivel de sign grm", "Nível de sign gpcm", "Cromossomo")
    #View(sign_gene_cr13)
    ##### UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM GENE RARO É SIGNIFICATIVO
    ind_gene_raro_sign_bf_com_cov<- which(as.numeric(sign_gene_cr13[,2])<(0.05/(length(seque) + qtd_sucessos_grm_cov))|
                                            as.numeric(sign_gene_cr13[,3])<(0.05/(length(seque) + qtd_sucessos_gpcm_cov))|
                                            as.numeric(sign_gene_cr13[,4])<(0.05/(length(seque) + qtd_sucessos_grm_cov))|
                                            as.numeric(sign_gene_cr13[,5])<(0.05/(length(seque) + qtd_sucessos_gpcm_cov)))
    
    
    gene_raro_sign_bf_com_cov<- as.matrix(sign_gene_cr13[ind_gene_raro_sign_bf_com_cov,])
    #View(gene_raro_sign_bf_com_cov)
    
    
    
    
    
    
    
    
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
    
    #View(gene_raro_sign_bf_com_cov) #resultado gene a gene
    #View(gene_raro_sign_bf_com_cov_agrup_final) #resultado com genes agrupados para SNPs q não rodavam
    if (dim(as.matrix(gene_raro_sign_bf_sem_cov))[2]==1){
      gene_raro_sign_bf_sem_cov<- matrix(gene_raro_sign_bf_sem_cov,nrow=1)
    }
    #cat('',gene_raro_sign_bf_sem_cov,file=paste("gene_raro_sem_cov_repli",repli,".txt",sep=""),append=T)
    write.table(
      gene_raro_sign_bf_sem_cov,
      file = paste0("gene_raro_sem_cov_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    
    if (dim(as.matrix(gene_raro_sign_bf_sem_cov_agrup_final))[2]==1){
      gene_raro_sign_bf_sem_cov_agrup_final<- matrix(gene_raro_sign_bf_sem_cov_agrup_final,nrow=1)
    }
    #cat('',gene_raro_sign_bf_sem_cov_agrup_final,file=paste("gene_raro_sem_cov_agrup_repli",repli,".txt",sep=""),append=T)
    write.table(
      gene_raro_sign_bf_sem_cov_agrup_final,
      file = paste0("gene_raro_sem_cov_agrup_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    
    
    if (dim(as.matrix(gene_raro_sign_bf_com_cov))[2]==1){
      gene_raro_sign_bf_com_cov<- matrix(gene_raro_sign_bf_com_cov,nrow=1)
    }
    #cat('',gene_raro_sign_bf_com_cov,file=paste("gene_raro_com_cov_repli",repli,".txt",sep=""),append=T)
    write.table(
      gene_raro_sign_bf_com_cov,
      file = paste0("gene_raro_com_cov_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    
    
    if (dim(as.matrix(gene_raro_sign_bf_com_cov_agrup_final))[2]==1){
      gene_raro_sign_bf_com_cov_agrup_final<- matrix(gene_raro_sign_bf_com_cov_agrup_final,nrow=1)
    }
    #cat('',gene_raro_sign_bf_com_cov_agrup_final,file=paste("gene_raro_com_cov_agrup_repli",repli,".txt",sep=""),append=T)
    write.table(
      gene_raro_sign_bf_com_cov_agrup_final,
      file = paste0("gene_raro_com_cov_agrup_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    
    
    #### exportando os coef a e b 
    tabelas_a = list(summary(a_grm),
                     summary(a_gpcm),
                     summary(a_grm_cov),
                     summary(a_gpcm_cov))
    
    nomes_tabelas_a = c("Summary a grm sem cov",
                        "Summary a gpcm sem cov",
                        "Summary a grm com cov",
                        "Summary a gpcm com cov")
    
    cat('Cromossomo',cr,"\n",file=paste("summary_coef_a_repli",repli,".txt",sep=""),append=T)
    
    for (i in seq_along(tabelas_a)) {
      cat(nomes_tabelas_a[i], "\n", file=paste("summary_coef_a_repli",repli,".txt",sep=""),append=T)
      capture.output(print(tabelas_a[[i]], row.names = FALSE), file=paste("summary_coef_a_repli",repli,".txt",sep=""),append=T)
      cat("\n", file=paste("summary_coef_a_repli",repli,".txt",sep=""),append=T)
    }
    # Linha em branco extra entre cromossomos
    cat("\n", file = paste("summary_coef_a_repli",repli,".txt",sep=""),append=T)
    
    #### AGORA AGRUPADO
    tabelas_a_agrup = list(summary(a_grm_agrup),
                           summary(a_gpcm_agrup),
                           summary(a_grm_cov_agrup),
                           summary(a_gpcm_cov_agrup))
    
    
    
    nome_tabelas_a_agrup = c("Summary a grm sem cov agrup",
                             "Summary a gpcm sem cov agrup",
                             "Summary a grm com cov agrup",
                             "Summary a gpcm com cov agrup")
    
    cat('Cromossomo',cr,"\n",file=paste("summary_coef_a_agrup_repli",repli,".txt",sep=""),append=T)
    for (i in seq_along(tabelas_a_agrup)) {
      cat(nome_tabelas_a_agrup[i], "\n", file=paste("summary_coef_a_agrup_repli",repli,".txt",sep=""),append=T)
      capture.output(print(tabelas_a_agrup[[i]], row.names = FALSE), file=paste("summary_coef_a_agrup_repli",repli,".txt",sep=""),append=T)
      cat("\n", file=paste("summary_coef_a_agrup_repli",repli,".txt",sep=""),append=T)
    }
    
    # Linha em branco extra entre cromossomos
    cat("\n", file=paste("summary_coef_a_agrup_repli",repli,".txt",sep=""),append=T)
    
    
    #### exportando os coef a e b 
    tabelas_b = list(summary(b1_grm),
                     summary(b1_gpcm),
                     summary(b1_grm_cov),
                     summary(b1_gpcm_cov))
    
    nomes_tabelas_b = c("Summary b1 grm sem cov",
                        "Summary b1 gpcm sem cov",
                        "Summary b1 grm com cov",
                        "Summary b1 gpcm com cov")
    
    cat('Cromossomo',cr,"\n",file=paste("summary_coef_b1_repli",repli,".txt",sep=""),append=T)
    for (i in seq_along(tabelas_b)) {
      cat(nomes_tabelas_b[i], "\n", file=paste("summary_coef_b1_repli",repli,".txt",sep=""),append=T)
      capture.output(print(tabelas_b[[i]], row.names = FALSE), file=paste("summary_coef_b1_repli",repli,".txt",sep=""),append=T)
      cat("\n", file=paste("summary_coef_b1_repli",repli,".txt",sep=""),append=T)
    }
    # Linha em branco extra entre cromossomos
    cat("\n", file=paste("summary_coef_b1_repli",repli,".txt",sep=""),append=T)
    
    
    
    
    tabelas_b_agrup <- list(summary(b1_grm_agrup),
                            summary(b1_gpcm_agrup),
                            summary(b1_grm_cov_agrup),
                            summary(b1_gpcm_cov_agrup))
    nomes_tabelas_b_agrup <- c("Summary b1 grm sem cov agrup",
                               "Summary b1 gpcm sem cov agrup",
                               "Summary b1 grm com cov agrup",
                               "Summary b1 gpcm com cov agrup")
    
    cat('Cromossomo',cr,"\n",file=paste("summary_coef_b1_agrup_repli",repli,".txt",sep=""),append=T)
    for (i in seq_along(tabelas_b_agrup)) {
      cat(nomes_tabelas_b_agrup[i], "\n", file=paste("summary_coef_b1_agrup_repli",repli,".txt",sep=""),append=T)
      capture.output(print(tabelas_b_agrup[[i]], row.names = FALSE), file=paste("summary_coef_b1_agrup_repli",repli,".txt",sep=""),append=T)
      cat("\n", file=paste("summary_coef_b1_agrup_repli",repli,".txt",sep=""),append=T)
    }
    # Linha em branco extra entre cromossomos
    
    cat("\n", file=paste("summary_coef_b1_agrup_repli",repli,".txt",sep=""),append=T)
    
  }
}