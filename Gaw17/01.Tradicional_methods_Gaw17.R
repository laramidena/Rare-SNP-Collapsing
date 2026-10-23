library(dplyr)
library(psych)
library(readr)
#library(ggplot2)
library(grDevices)
library(SKAT)
library(glmnet)
library(compiler)

######## LEITURA DA BASE DE DADOS DO CROMOSSOMO 1 E DO FENÓTIPO 


snp_info <- read_csv("~/Dados/snp_info")

enableJIT(3)
for (repli in c(1:10)){
  enableJIT(3)
  for (cr in c(1:22)){
    unr_phen <- read_csv(paste0("~/Dados/unr_phen.",repli))
    c1_snps <- read_csv(paste0("~/Dados/c",cr, "_snps.unr"))
    fuma<- unr_phen$SMOKE
    idade<- unr_phen$AGE
    sexo<- unr_phen$SEX-1
    
    #boxplot(unr_phen$Q1, col = "lightblue", xlab = "Q1")
    #boxplot.stats(unr_phen$Q1)$out
    
    ############# RETIRANDO SNPS SEM VARIAÇÃO NO CROMOSSOMO 1
    mono_1<-numeric(ncol(c1_snps)-1)
    for (i in 1:(ncol(c1_snps)-1)) mono_1[i]<-length(table(c1_snps[,i+1]))
    checa<-table(mono_1)
    j<-1
    excl<-numeric(checa[1])
    for (i in 1:(ncol(c1_snps)-1)){
      if (mono_1[i]==1) {
        excl[j]<-i
        j<-j+1}}
    excl<-excl+1
    c1_snps<-c1_snps[,-excl] 
    #ncol(c1_snps) 
    
    
    ######## CLASSIFICACAO DAS BASE DE DADOS EM 0, 1 E 2 
    dados_1<-matrix(0,nrow=nrow(c1_snps),ncol = ncol(c1_snps))
    enableJIT(3)
    for (i in 1:nrow(dados_1)){
      for (j in 1: ncol(dados_1)){
        if (c1_snps[i,j]=="A/A"|c1_snps[i,j]=="T/T") dados_1[i,j]=2
        if (c1_snps[i,j]=="G/G"|c1_snps[i,j]=="C/C") dados_1[i,j]=0
        if (c1_snps[i,j]!="A/A" & c1_snps[i,j]!="T/T" &
            c1_snps[i,j]!="G/G" & c1_snps[i,j]!="C/C") dados_1[i,j]=1
      }
    }
    #View(dados_1)
    #ncol(dados_1)
    #table(dados_1[,7])
    nomes_1 = colnames(c1_snps)
    colnames(dados_1)<- nomes_1
    
    
    ######## VENDO A FREQ ALÉLICA DE TODOS OS SNPS DO CR1
    freq_alelica_SNPs_cr1 <-matrix(0,nrow=2,ncol = ncol(c1_snps))
    row.names(freq_alelica_SNPs_cr1)<- c('freq alelo A','freq alelo a')
    colnames(freq_alelica_SNPs_cr1)<- colnames(c1_snps)
    #View(freq_alelica_SNPs_cr1)
    enableJIT(3)
    for (i in 1:ncol(dados_1)){
      freq_alelica_SNPs_cr1[1,i]<- 
        (0*length(dados_1[which(dados_1[,i]=='2'),i])/(2*nrow(dados_1))) + 
        (1*length(dados_1[which(dados_1[,i]=='1'),i])/(2*nrow(dados_1))) +
        (2*length(dados_1[which(dados_1[,i]=='0'),i])/(2*nrow(dados_1)))
      
      freq_alelica_SNPs_cr1[2,i]<- 
        (2*length(dados_1[which(dados_1[,i]=='2'),i])/(2*nrow(dados_1))) + 
        (1*length(dados_1[which(dados_1[,i]=='1'),i])/(2*nrow(dados_1))) +
        (0*length(dados_1[which(dados_1[,i]=='0'),i])/(2*nrow(dados_1)))
      
    }
    #View(freq_alelica_SNPs_cr1)
    
    ######## SEPARANDO OS SNPS EM RAROS E COMUNS E IDENTIFICANDO-OS
    indice_SNPs_raros_cr_1<-which(freq_alelica_SNPs_cr1[,1:ncol(dados_1)]<0.05)
    for (i in 1:length(indice_SNPs_raros_cr_1)){
      if (indice_SNPs_raros_cr_1[i]%%2==0) indice_SNPs_raros_cr_1[i] = indice_SNPs_raros_cr_1[i]/2
      else indice_SNPs_raros_cr_1[i] = (indice_SNPs_raros_cr_1[i]+1)/2
    }
    
    SNPs_raros_cr1<- dados_1[,indice_SNPs_raros_cr_1]
    #View(SNPs_raros_cr1)
    colnames(SNPs_raros_cr1)<- nomes_1[indice_SNPs_raros_cr_1]
    SNPs_comuns_cr1<- dados_1[,-c(indice_SNPs_raros_cr_1)]
    #View(SNPs_comuns_cr1)
    
    
    ######## CALCULANDO O MAF DE CADA SNP DO CR 1 PRA FAZER BOXPLOT
    #DO MAF POR CROMOSSOMO
    maf_1<- NULL
    for (i in 1:ncol(freq_alelica_SNPs_cr1)){
      maf_1[i]<- min(freq_alelica_SNPs_cr1[1,i],freq_alelica_SNPs_cr1[2,i])
    }
    
    ######## FENÓTIPO
    ytr_1 = NULL #cria uma lista vazia
    
    for (i in 1:nrow(unr_phen)){
      ytr_1[i] = unr_phen$Q1[i]
    }
    #ytr_1
    
    
    ######## SEPARANDO OS SNPS RAROS DE ACORDO COM O GENE
    gene_1<- matrix(0,2,ncol(dados_1))
    
    ## ARRUMAR AQUI - DEFESA
    tab_nomes_1<-as.matrix(colnames(SNPs_raros_cr1))
    colnames(tab_nomes_1)<-c("snp_name")
    tab_nomes_1<-as.data.frame(tab_nomes_1)
    SNPs_raros_cr1<-as.data.frame(SNPs_raros_cr1)
    tabela_full_cr1<- tab_nomes_1 %>% left_join(snp_info,by='snp_name')
    tabela_gene_cr1<- tabela_full_cr1[,-c(2,3,5,6,7)]
    
    
    
    ######## AGRUPANDO OS SNPS Q ESTÃO NO MESMO GENE
    genes_cr1<-as.character(data.frame(table(tabela_gene_cr1[,2]))[,1])
    
    
    ######### PRIMEIRA ABORDAGEM: APLICA AS METODOLOGIAS NO MODELO
    #AJUSTADO APENAS COM O INTERCEPTO
    
    
    ######## TESTANDO SIGNIFICANCIA DOS GENES DE SNPS RAROS PARA CADA
    #METODOLOGIA EM UM MODELO SÓ COM O INTERCEPTO
    valorp_Burden<-NULL
    valorp_SKAT<-NULL
    valorp_SKAT_O<-NULL
    SKAT_O<-NULL
    rho_SKAT_O<-NULL
    #valorp_gene<- matrix(0,nrow=length(genes),ncol=2)
    for (i in 1:length(genes_cr1)){
      snps<-which(tabela_gene_cr1[,2]==genes_cr1[i]) # posição na lista genes 
      #(tá em ordem alfabética) dos snps que estão no gene pedido da tabela_gene_cr1
      obj<-SKAT_Null_Model(ytr_1 ~ 1, out_type="C")
      valorp_Burden[i]<-SKAT(as.matrix(SNPs_raros_cr1[,snps]), obj, r.corr=1)$p.value 
      #ordem dos genes em ordem alfabética
      valorp_SKAT[i]<-SKAT(as.matrix(SNPs_raros_cr1[,snps]), obj, r.corr=0)$p.value
      SKAT_O<-SKAT(as.matrix(SNPs_raros_cr1[,snps]), obj, method="optimal.adj")
      if (is.null(SKAT_O$param$minp)==FALSE){
        valorp_SKAT_O[i]<- SKAT_O$param$minp
        rho_SKAT_O[i]<-SKAT_O$param$rho_est
      }
    }
    sign_gene_cr1<-cbind(genes_cr1, valorp_Burden, valorp_SKAT, valorp_SKAT_O,
                         rho_SKAT_O,repli,cr)
    colnames(sign_gene_cr1)<- c("Genes", "Valor-p Burden", "Valor-p SKAT",
                                "Valor-p SKAT-O", "Rho SKAT-O","Replica","Cromossomo")
    
    
    
    #View(sign_gene_cr1)
    
    ######### SEGUNDA ABORDAGEM: SELECIONA AS COVARIÁVEIS PELO LASSO E APLICA
    #AS METODOLOGIAS NO MODELO AJUSTADO COM AS COVARIÁVEIS SELECIONADAS PELO LASSO
    
    
    
    ######## SELECIONANDO SNPS COMUNS PELO LASSO
    set.seed(202310)
    x_cr1<- as.matrix(SNPs_comuns_cr1)
    X_cr1<- cbind(x_cr1,fuma,idade,sexo)
    cv.out_cr1 <- cv.glmnet(as.matrix(X_cr1), ytr_1, alpha = 1, 
                            family="gaussian")
    plot(cv.out_cr1, ylab = "EQP", xlab = "log(lambda)")
    bestlam_1_cr1 <- cv.out_cr1$lambda.min
    
    ### ARMAZENA O MELHOR LÂMBIDA E O GRAU DE LIBERDADE
    lasso.mod_cr1 <- glmnet(as.matrix(X_cr1), ytr_1, alpha = 1, lambda = bestlam_1_cr1,
                            family="gaussian")
    lasso.coef_cr1  <- predict(lasso.mod_cr1, type = 'coefficients', s = bestlam_1_cr1)
    coef_lasso_cromo_1 = as.matrix(lasso.coef_cr1)
    
    ######## REMOVE O VALOR DO INTERCEPTO NO LASSO E DETERMINA QUAIS VARIÁVEIS SÃO
    #SIGNIFICATIVAS
    betas_cr1<- matrix(0,nrow=ncol(X_cr1),ncol=2)
    betas_cromo_1 = coef_lasso_cromo_1[c(-1)]
    betas_cr1[,1]<-colnames(X_cr1)
    betas_cr1[,2]<- betas_cromo_1
    snps_comuns_sign_cr1 = betas_cr1[which(betas_cr1[,2]!=0),]
    snps_comuns_sign_cr1<- as.data.frame(snps_comuns_sign_cr1)
    snps_comuns_sign_cr1 <- cbind(snps_comuns_sign_cr1,repli,cr)
    
    colnames(snps_comuns_sign_cr1)<- c('Variável comum','Coeficiente beta','Replica','Cromossomo')
    indice_var_comum_sign_cr1 = which(betas_cr1[,2]!=0)
    var_comum_sign_cr1<- X_cr1[,indice_var_comum_sign_cr1]
    
    #View(var_comum_sign_cr1)
    
    ######## TESTANDO SIGNIFICANCIA DOS GENES DE SNPS RAROS USANDO COMO COVARIÁVEIS
    #OS SNPS COMUNS E VARIÁVEIS SELECIONADAS PELO LASSO
    valorp_Burden_cov<-NULL
    valorp_SKAT_cov<-NULL
    valorp_SKAT_O_cov<-NULL
    SKAT_O_cov<-NULL
    rho_SKAT_O_cov<-NULL
    #valorp_gene<- matrix(0,nrow=length(genes),ncol=2)
    for (i in 1:length(genes_cr1)){
      snps<-which(tabela_gene_cr1[,2]==genes_cr1[i]) #posição na lista genes
      #(está em ordem alfabética) dos snps que estão no gene pedido da tabela_gene_cr1
      obj<-SKAT_Null_Model(ytr_1 ~ var_comum_sign_cr1, out_type="C")
      valorp_Burden_cov[i]<-SKAT(as.matrix(SNPs_raros_cr1[,snps]),
                                 obj, r.corr=1)$p.value #ordem dos genes em ordem alfabética
      valorp_SKAT_cov[i]<-SKAT(as.matrix(SNPs_raros_cr1[,snps]), obj,
                               r.corr=0)$p.value
      SKAT_O_cov<-SKAT(as.matrix(SNPs_raros_cr1[,snps]), obj,   method="optimal.adj")
      if (is.null(SKAT_O_cov$param$minp)==FALSE){
        valorp_SKAT_O_cov[i]<- SKAT_O_cov$param$minp
        rho_SKAT_O_cov[i]<-SKAT_O_cov$param$rho_est
      }
    }
    sign_gene_cr1_cov<-cbind(genes_cr1, valorp_Burden_cov, valorp_SKAT_cov,
                             valorp_SKAT_O_cov, rho_SKAT_O_cov,repli,cr)
    colnames(sign_gene_cr1_cov)<- c("Genes", "Valor-p Burden com Cov",
                                    "Valor-p SKAT com Cov", "Valor-p SKAT-O com Cov", "Rho SKAT-O com Cov",
                                    "Replica","Cromossomo")
    
    
    #View(sign_gene_cr1_cov)
    
    
    
    ######## UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM GENE RARO
    #É SIGNIFICATIVO NO TESTE DE BURDEN
    ind_gene_raro_sign_bf_sem_cov_cr1_burden<- 
      which(as.numeric(sign_gene_cr1[,2])<0.05/length(genes_cr1))
    
    gene_raro_sign_bf_sem_cov_cr1_burden<- 
      as.matrix(sign_gene_cr1[ind_gene_raro_sign_bf_sem_cov_cr1_burden,c(1,2,6,7)])
    
    #View(gene_raro_sign_bf_sem_cov_cr1_burden)
    
    if (dim(as.matrix(gene_raro_sign_bf_sem_cov_cr1_burden))[2]==1){
      nomes<-rownames(gene_raro_sign_bf_sem_cov_cr1_burden)
      gene_raro_sign_bf_sem_cov_cr1_burden<- matrix(gene_raro_sign_bf_sem_cov_cr1_burden,nrow=1)
      colnames(gene_raro_sign_bf_sem_cov_cr1_burden)<-nomes
    }
    
    ind_gene_raro_sign_bf_com_cov_cr1_burden<- 
      which(as.numeric(sign_gene_cr1_cov[,2])<0.05/length(genes_cr1))
    gene_raro_sign_bf_com_cov_cr1_burden<- 
      as.matrix(sign_gene_cr1_cov[ind_gene_raro_sign_bf_com_cov_cr1_burden,c(1,2,6,7)])
    
    
    #View(gene_raro_sign_bf_com_cov_cr1_burden)
    
    
    if (dim(as.matrix(gene_raro_sign_bf_com_cov_cr1_burden))[2]==1){
      nomes<-rownames(gene_raro_sign_bf_com_cov_cr1_burden)
      gene_raro_sign_bf_com_cov_cr1_burden<- matrix(gene_raro_sign_bf_com_cov_cr1_burden,nrow=1)
      colnames(gene_raro_sign_bf_com_cov_cr1_burden)<-nomes
    }
    
    ######## UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM
    #GENE RARO É SIGNIFICATIVO NO SKAT
    ind_gene_raro_sign_bf_sem_cov_cr1_skat<-
      which(as.numeric(sign_gene_cr1[,3])<0.05/length(genes_cr1))
    gene_raro_sign_bf_sem_cov_cr1_skat<-
      as.matrix(sign_gene_cr1[ind_gene_raro_sign_bf_sem_cov_cr1_skat,c(1,3,6,7)])
    
    #View(gene_raro_sign_bf_sem_cov_cr1_skat)
    
    if (dim(as.matrix(gene_raro_sign_bf_sem_cov_cr1_skat))[2]==1){
      nomes<-rownames(gene_raro_sign_bf_sem_cov_cr1_skat)
      gene_raro_sign_bf_sem_cov_cr1_skat<- matrix(gene_raro_sign_bf_sem_cov_cr1_skat,nrow=1)
      colnames(gene_raro_sign_bf_sem_cov_cr1_skat)<-nomes
    }
    
    
    ind_gene_raro_sign_bf_com_cov_cr1_skat<- 
      which(as.numeric(sign_gene_cr1_cov[,3])<0.05/length(genes_cr1))
    gene_raro_sign_bf_com_cov_cr1_skat<- 
      as.matrix(sign_gene_cr1_cov[ind_gene_raro_sign_bf_com_cov_cr1_skat,c(1,3,6,7)])
    
    
    #View(gene_raro_sign_bf_com_cov_cr1_skat)
    
    if (dim(as.matrix(gene_raro_sign_bf_com_cov_cr1_skat))[2]==1){
      nomes<-rownames(gene_raro_sign_bf_com_cov_cr1_skat)
      gene_raro_sign_bf_com_cov_cr1_skat<- matrix(gene_raro_sign_bf_com_cov_cr1_skat,nrow=1)
      colnames(gene_raro_sign_bf_com_cov_cr1_skat)<-nomes
    }
    
    ######## UTILIZANDO A CORREÇÃO DE BONFERRONI E VENDO SE ALGUM 
    #GENE RARO É SIGNIFICATIVO NO SKAT-O
    ind_gene_raro_sign_bf_sem_cov_cr1_skat_o<- 
      which(as.numeric(sign_gene_cr1[,4])<0.05/length(genes_cr1))
    gene_raro_sign_bf_sem_cov_cr1_skat_o<- 
      as.matrix(sign_gene_cr1[ind_gene_raro_sign_bf_sem_cov_cr1_skat_o,c(1,4,5,6,7)])
    
    #View(gene_raro_sign_bf_sem_cov_cr1_skat_o)
    
    if (dim(as.matrix(gene_raro_sign_bf_sem_cov_cr1_skat_o))[2]==1){
      nomes<-rownames(gene_raro_sign_bf_sem_cov_cr1_skat_o)
      gene_raro_sign_bf_sem_cov_cr1_skat_o<- matrix(gene_raro_sign_bf_sem_cov_cr1_skat_o,nrow=1)
      colnames(gene_raro_sign_bf_sem_cov_cr1_skat_o)<-nomes
    }
    
    ind_gene_raro_sign_bf_com_cov_cr1_skat_o<- 
      which(as.numeric(sign_gene_cr1_cov[,4])<0.05/length(genes_cr1))
    gene_raro_sign_bf_com_cov_cr1_skat_o<- 
      as.matrix(sign_gene_cr1_cov[ind_gene_raro_sign_bf_com_cov_cr1_skat_o,c(1,4,5,6,7)])
    
    #View(gene_raro_sign_bf_com_cov_cr1_skat_o)
    
    if (dim(as.matrix(gene_raro_sign_bf_com_cov_cr1_skat_o))[2]==1){
      nomes<-rownames(gene_raro_sign_bf_com_cov_cr1_skat_o)
      gene_raro_sign_bf_com_cov_cr1_skat_o<- matrix(gene_raro_sign_bf_com_cov_cr1_skat_o,nrow=1)
      colnames(gene_raro_sign_bf_com_cov_cr1_skat_o)<-nomes
    }
    
    
    ####### RESULTADOS
    
    ##### PRIMEIRA ABORDAGEM: MODELO AJUSTADO SÓ COM O INTERCEPTO
    
    ##### EXPORTANDO
    write.table(
      gene_raro_sign_bf_sem_cov_cr1_burden,
      file = paste0("resul_modelo_s_cov_burden_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    write.table(
      gene_raro_sign_bf_sem_cov_cr1_skat,
      file = paste0("resul_modelo_s_cov_skat_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    write.table(
      gene_raro_sign_bf_sem_cov_cr1_skat_o,
      file = paste0("resul_modelo_s_cov_skat_o_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    
    ##### SEGUNDA ABORDAGEM: MODELO AJUSTADO COM AS COVARIÁVEIS 
    #SELECIONADAS PELO LASSO
    
    write.table(
      gene_raro_sign_bf_com_cov_cr1_burden,
      file = paste0("resul_modelo_c_cov_burden_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    write.table(
      gene_raro_sign_bf_com_cov_cr1_skat,
      file = paste0("resul_modelo_c_cov_skat_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    write.table(
      gene_raro_sign_bf_com_cov_cr1_skat_o,
      file = paste0("resul_modelo_c_cov_skat_o_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
    
    
    
    #### RESULTADOS DAS VARIÁVEIS E SNPS COMUNS SELECIONADOS COMO 
    #SIGNIFICATIVASPELO LASSO
    write.table(
      snps_comuns_sign_cr1,
      file = paste0("result_var_comuns_sign_repli", repli, ".txt"),
      sep = ",",
      row.names = FALSE,
      col.names = (cr == 1),  # só escreve cabeçalho na 1ª vez
      append = (cr != 1)
    )
  }
}

1+1

