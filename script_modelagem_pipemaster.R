#Análises de modelagem demográfica
#carregando pacotes necessários
library(PipeMaster)
library(ggplot2)
library(matrixStats)

#1.0-diretório de trabalho
work_dir=paste("/media/lucatristaom/30FEAC51FEAC1160/Backup_linux/data/models/newmodels/")
#1.1 OS ARQUIVOS FASTA ESTÃO NO LINUX
fasta_files=paste("/home/lucatristaom/Desktop/Mestrado/Backup/backup_servidor/phasing/alignments/new_1_alignments/")

#Windows
work_dir=paste ("D:/Backup_linux/data/models/newmodels/")
#1.2-designando populações a priori e recuperando estatísticas sumárias
setwd(work_dir)
pops <- read.delim(paste("pop_assign_k3.txt", sep = ""), header = FALSE, sep = "\t")

obs_stats <- obs.sumstat.ngs(model = M1, path.to.fasta = "./", 
                        pop.assign = pops, moments = F)

#2-criando modelos ou recuperando, caso já estejam criados, pular para seção 3. 
#topologia seguida é a da árvore filogenética ((1,2),3)

#Isolamento completo sem fluxo gênico
M1 <- main.menu(M1)
M1 <- dget("M1.txt", file="M1.txt")
#Isolamento com migração entre 2 e 3 (leste e oeste), contato secundário.
M2 <- main.menu(M2)
M2 <- dget ("M2.txt", file="M2.txt")
#Isolamento com migração em todas as populações.
M3 <- main.menu(M3)
M3 <- dget("M3.txt", file="M3.txt")
#Isolamento sem migração e com expansão populacional.
M4 <- main.menu(M4)
M4 <- dget("M4.txt", file="M4.txt")
#Isolamento com expansão populacional e migração entre todas as populações.
M5 <- main.menu(M5)
M5 <- dget("M5.txt", file="M5.txt")
#Isolamento sem migração e com gargalo populacional
M6 <- main.menu()
M6 <- dget("M6.txt", file="M6.txt")
#Isolamento com migração entre todas as populações e com gargalo populacional
M7 <- main.menu()
M7 <-dget("M7.txt",file="M7.txt")

#3-inserindo os loci nos modelos.
  setwd(work_dir)
  M1 <- get.data.structure(M1, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)
  M2 <- get.data.structure(M2, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)
  M3 <- get.data.structure(M3, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)
  M4 <- get.data.structure(M4, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)
  M5 <- get.data.structure(M5, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)
  M6 <- get.data.structure(M6, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)
  M7 <- get.data.structure(M7, path.to.fasta = "./", pop.assign = pops, sanger=FALSE)

#3.2 salvando os modelos com dados genômicos carregados
if (file.exits("M1.txt")==F){
  dput(M1, file="M1.txt")
  dput(M2, file="M2.txt")
  dput(M3, file="M3.txt")
  dput(M4, file="M4.txt")
  dput(M5, file="M5.txt")
  dput(M6, file="M6.txt")
  dput(M7, file="M7.txt")
  
  print("Modelos com dados genômicos salvos")
}  
 
#3.3 recuperando todos os modelos com a estrutura
if (file.exists("M1.txt")==T){
  M1 <- dget("M1.txt")
  M2 <- dget("M2.txt")
  M3 <- dget("M3.txt")
  M4 <- dget("M4.txt")
  M5 <- dget("M5.txt")
  M6 <- dget("M6.txt")
  M7 <- dget("M7.txt")
  
  print("Modelos carregados com dados genômicos")
}

#4 Obtendo estatísticas descritivas dos dados.
setwd(work_dir)
if (file.exists("stats_obs2.csv")==F) {
  stat <- obs.sumstat.ngs(model = M1, path.to.fasta = "./", 
                          pop.assign = pop_assign_k3, moments = F)
  write.csv(stat, "stats_obs2.csv")
} else {
  obs_stat <- read.csv("stats_obs2.csv")
}

#4.1 Selecionando apenas as estatísticas desejáveis.
cols <- c(grep("average_segs", colnames(obs_stat)),
          grep("average_pairwise_fst_", colnames(obs_stat)),
          grep("average_pi", colnames(obs_stat)),
          grep("average_w", colnames(obs_stat)),
          grep("average_tajd", colnames(obs_stat)),
          grep("average_shared_", colnames(obs_stat)),
          grep("average_private_", colnames(obs_stat)),
          grep("average_FayWuH", colnames(obs_stat)))

obs_stat_corrigido <- t(data.frame(obs_stat[,-cols]))

#4.2 salvando as estatísticas selecionadas
write.table(obs_stat_corrigido, "Observed_data_now.txt", quote=F, row.names=F, col.names=F)

#4.3 carregando as estatísticas de interesse
tabela_nova_stats <- read.table("Observed_data_now.txt", header=T)

#5. rodando a modelagem em loop
setwd(work_dir)

if (file.exists("SIMS_M1.txt")==F){
  #M1 Isolamento completo.
  sim.msABC.sumstat(model=M1, path="./",
                    nsim.blocks = 25, use.alpha = F,
                    output.name = "M1",
                    append.sims = F, ncores=7, block.size = 300)
  
  #M2 Isolamento com migração (contato secundário) entre leste e oeste.
  sim.msABC.sumstat(M2, path="./",
                    nsim.blocks = 25, use.alpha = F,
                    output.name = "M2",
                    append.sims = F, ncores=7, block.size = 300)
  #M3 Isolamento com migração igual para todas as populações.
  sim.msABC.sumstat(M3, path="./",
                    nsim.blocks = 25, use.alpha=F,
                    output.name="M3",
                    append.sims=F, ncores=7, block.size=300)
  #M4. Isolamento com migração igual para todas as pops.
  sim.msABC.sumstat(M4, path="./",
                    nsim.blocks=25, use.alpha=F,
                    output.name="M4",
                    append.sims=F, ncores=7, block.size = 300)
  
  #M5. Isolamento com migração e expansão populacional.
  sim.msABC.sumstat(M5, path="./",
                    nsim.blocks=25, use.alpha=F,
                    output.name="M5", 
                    append.sims=F, ncores=7, block.size=300)
  
  #M6. Isolamento com migração e gargalo populacional.
  sim.msABC.sumstat(M6, path="./",
                    nsim.blocks=25, use.alpha=F,
                    output.name="M6",
                    append.sims=F, ncores=7, block.size=300)
  #M7. Isolamento com gargalo populacional e migração entre todas as populações.
  sim.msABC.sumstat(M7, path="./",
                    nsim.blocks=25, use.alpha=F,
                    output.name="M7",
                    append.sims=F, ncores=7, block.size=300)
}

#5.1 LENDO AS SIMULAÇÕES
setwd(work_dir)
M1_sim_est <- read.delim("SIMS_M1.txt")
M2_sim_est <- read.delim("SIMS_M2.txt")
M3_sim_est <- read.delim("SIMS_M3.txt")
M4_sim_est <- read.delim("SIMS_M4.txt")
M5_sim_est <- read.delim("SIMS_M5.txt")
M6_sim_est <- read.delim("SIMS_M6.txt")
M7_sim_est <- read.delim("SIMS_M7.txt")

#5.2 Selecionando apenas as variáveis de interesse nas simulações
M1_sim <- M1_sim_est[,colnames(M1_sim_est) %in% colnames(tabela_nova_stats)]
M2_sim <- M2_sim_est[,colnames(M2_sim_est) %in% colnames(tabela_nova_stats)]
M3_sim <- M3_sim_est[,colnames(M3_sim_est) %in% colnames(tabela_nova_stats)]
M4_sim <- M4_sim_est[,colnames(M4_sim_est) %in% colnames(tabela_nova_stats)]
M5_sim <- M5_sim_est[,colnames(M5_sim_est) %in% colnames(tabela_nova_stats)]
M6_sim <- M6_sim_est[,colnames(M6_sim_est) %in% colnames(tabela_nova_stats)]
M7_sim <- M7_sim_est[,colnames(M7_sim_est) %in% colnames(tabela_nova_stats)]

#5.3 concatenando simulações com variáveis de interesse
setwd(work_dir)
modelos_concatenados_reduzidos <- rbind(M1_sim, M2_sim, M3_sim, M4_sim, M5_sim, M6_sim, M7_sim)

#5.4 salvando simulações concatenados
write.csv(modelos_concatenados_reduzidos, file="modelos_concatenados_variaveis_reduzidas.csv")

#5.5 lendo as simulações
modelos_concatenados_reduzidos <- read.csv("modelos_concatenados_variaveis_reduzidas.csv", header=TRUE, sep=",")

#5.6 calculando um índice com as simulações reduzidas
index <- c(rep("M1", nrow(M1_sim)), rep("M2", nrow(M2_sim)), rep("M3", nrow(M3_sim)), rep("M4", nrow(M4_sim)), rep("M5", nrow(M5_sim)), rep("M6", nrow(M6_sim)), rep("M7", nrow(M7_sim)))

#6 PCA para visualização do fit das simulações e dados observados
PCA_models <- plotPCs(model=modelos_reduzidos, index=index, observed=tabela_nova_stats, subsample=1)
PCA_models

#ABC#
#7 Inferência de um melhor modelo pelo ABC.
prob <- postpr(target=tabela_nova_stats, sumstat=modelos_reduzidos, index=index, method="neuralnet", tol=0.1)
summary(prob)

#7.1 validações cruzadas para determinar confiança nas inferências
CV <- cv4postpr(sumstat=modelos_reduzidos, index=index, method="neuralnet", tol=0.1, nval=100)
summary(CV)
plot(CV)

#7.2 Precisão das validações
acc <- summary(CV)
sum(diag(acc$conf.matrix$tol0.1))/60

#8 ESTIMATIVA DE PARÂMETROS
#Selecionando somente os parâmetros estimados
param_M3 <- M3_sim_est[,1:9]

#8.1 Estimando posterioris dos parâmetros para o melhor modelo (M3)
posterior_M3 <- abc(target=tabela_nova_stats,
                    param=param_M3,
                    sumstat=M3_sim,
                    method="neuralnet",
                    tol=0.1)
summary(posterior_M3)

#8.2 Plotando distribuição das probabilidades a posteriori dos parâmetros em função dos priors
for(i in 1:ncol(param_M3)){
  plot(density(posterior_M3$unadj.values[,i]), col=2, main = colnames(param_M3)[i])
  lines(density(param_M3[,i]))
}
#8.3 validação cruzada das estimativas dos parâmetros
cv_abc_M3 <- cv4abc(param = param_M3,
                    sumstat = M3_sim,
                    nval = 100,
                    method = "neuralnet",
                    tol = 0.1)

plot(cv_abc_M3)
summary(cv_abc_M3)