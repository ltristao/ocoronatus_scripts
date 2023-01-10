#SCRIPT PARA ANÁLISES DE ESTRUTURA POPULACIONAL E GENÉTICA.

#BIBLIOTECAS UTILIZADAS
library("rnaturalearth")
library("tidyverse")
library("LEA")
library("adegenet")
library("vcfR")
library("ggplot2")
library("M3C")
library("genepop")

#INSERINDO DADOS EM FORMATO VCF PARA PRESERVAR INFORMAÇÃO DE INDIVÍDUOS
vcf_ocoronatus2 <- read.vcfR("D:/Backup_linux/data/snmf/results/snmf_rebuild_sem_swainsoni/arquivo_final_sem_swainsoni_para_snmf.vcf")

#INSERINDOS DADOS DE LOCALIDADE E INFORMATIVOS SOBRE CADA INDIVÍDUO
amostras_sem_swainsoni <- read.csv("C:/Users/lucaf/Dropbox/PC/Desktop/Faculdade/Mestrado/Dissertação/data/snmf_results/planilha localizacao_AEndemism_csv_atualizado_com_localidades_modificado_sem_contaminados_sem_swainsoni.csv", header=TRUE, sep=",")
as_tibble(amostras_sem_swainsoni)

#CONVERSÃO PARA FORMATO GENLIGHT 
ocoronatus_genlight2 <- vcfR2genlight(vcf_ocoronatus2)

#INFORMAÇÃO SUMÁRIA DO ARQUIVO VCF DISPONÍVEL EM AMBIENTE R
ocoronatus_genlight2

#GENPLOT
plot_alelos <- glPlot(ocoronatus_genlight2)
plot_alelos

#ANÁLISE DE COMPONENTES PRINCIPAIS
PCA_ocoronatus <- glPca(ocoronatus_genlight2)
scatter(PCA_ocoronatus) #3 populações visualmente encontradas no PCA.

#REFORMATAR PARA O GGPLOT
PCA_ocoronatus_tibble <- as_tibble(PCA_ocoronatus$scores, rownames="individuos") %>% 
  mutate(populacao=amostras_sem_swainsoni$Area_endemismo)
PCA_ocoronatus_tibble

#PLOTANDO O PCA
PCA_plot <- 
  ggplot(PCA_ocoronatus_tibble, aes(
    x=PC1,
    y=PC2,
    fill=populacao
  )) + 
  geom_point(shape=21, size=3, label=amostras_sem_swainsoni$N)+
  theme_bw(base_size=16)+
  scale_fill_manual(values=my.colorsdapc, labels=c("Belem","Guiana","Inambari","Napo","Nenhuma","Rondonia","Tapajós","Xingu"))

PCA_plot



#ANÁLISE DISCRIMINANTE

#ENCONTRAR CLUSTERS A PARTIR DE K-MEANS. INCLUIR O MÁXIMO DE NÚMERO DE PCS QUE CONTRIBUEM PARA A VARIÂNCIA.
num_clust <- find.clusters(ocoronatus_genlight2)
num_clust

#RODAR DAPC COM MELHOR CLUSTER, INCLUIR ATÉ 80%, MAIS QUE ISSO É EXCEÇÃO.
ocoronatus_dapc2 <- dapc(ocoronatus_genlight2, num_clust$grp)
ocoronatus_dapc2

#SCATTERPLOT DOS DADOS, RÁPIDA VISUALIZAÇÃO.
scatter(ocoronatus_dapc2, posi.da="bottomleft")

#TIDYVERSE, INCLUI DADOS ADICIONAIS NA MATRIZ
dapc_dados_df2 <- as_tibble(ocoronatus_dapc2$ind.coord, rownames="individuos") %>% 
  mutate(populacao=amostras_sem_swainsoni$Area_endemismo,
         group=ocoronatus_dapc2$grp)
dapc_dados_df2
dapc_dados_df2$populacao

#CORES PRO DAPC, CORRESPONDEM ÀS USADAS NO MAPA
my.colorsdapc = c("purple","red","blue","yellow","gray58","forestgreen","cyan","orange")

#PLOT DO DAPC
dapc_plot <- 
  ggplot(dapc_dados_df2, aes(
    x=LD1,
    y=LD2,
    fill=populacao
  )) + 
  geom_point(shape=21, size=3)+
  theme_bw(base_size=16)+
  scale_fill_manual(values=my.colorsdapc, labels=c("Belem","Guiana","Inambari","Napo","Nenhuma","Rondonia","Tapajós","Xingu"))

dapc_plot

#ESTATÍSTICA SUMÁRIA

vcfdnabin <- vcfR2DNAbin(vcf_ocoronatus2)
vcfdnabin

chromR_ocoronatus <- create.chromR(vcf_ocoronatus2,name="UCEs")
genepop_stats <- gt2popsum (chromR_ocoronatus)


#ANÁLISE DE ESTRUTURA GENÉTICA.
#CONVERTENDO VCF PARA FORMATO GENO.
vcf2geno(input.file='vcf_ocoronatus2.vcf', output.file='vcf_ocoronatus2.geno')

#SNMF
#CARREGANDO ARQUIVO LOG
setwd("/home/luca/Ocoronatus/snpcalling/alignments/vcftools/")
sink("sNMF_corridas_alpha_1_log.log")

#RODANDO A ANÁLISE. REPETIR PARA VALORES DE ALPHA = 10, 100, 200, 300, 400, 500 E 1000.
resultado_dissertacao = snmf("vcf_ocoronatus2.geno", K=1:10, entropy=TRUE, CPU=30, repetitions=100, alpha=1, project = "new")

#OBTENDO MELHOR CORRIDA PARA CADA "K", REPETIR PARA CADA VALOR DE ALPHA.
ce1 = cross.entropy(resultado_dissertacao, K = 1)
ce2 = cross.entropy(resultado_dissertacao, K = 2)
ce3 = cross.entropy(resultado_dissertacao, K = 3)
ce4 = cross.entropy(resultado_dissertacao, K = 4)
ce5 = cross.entropy(resultado_dissertacao, K = 5)
ce6 = cross.entropy(resultado_dissertacao, K = 6)
ce7 = cross.entropy(resultado_dissertacao, K = 7)
ce8 = cross.entropy(resultado_dissertacao, K = 8)
ce9 = cross.entropy(resultado_dissertacao, K = 9)
ce10 = cross.entropy(resultado_dissertacao, K = 10)
best1 = which.min(ce1)
best2 = which.min(ce2)
best3 = which.min(ce3)
best4 = which.min(ce4)
best5 = which.min(ce5)
best6 = which.min(ce6)
best7 = which.min(ce7)
best8 = which.min(ce8)
best9 = which.min(ce9)
best10 = which.min(ce10)

#PLOTANDO MELHOR "K" (NÚMERO DE POPULAÇÕES ANCESTRAIS), REPETIR PARA CADA VALOR DE ALPHA.
plot(resultado_dissertacao, lwd=5, col="red",pch=1, type="b")

#FORMATAÇÃO DAS MELHORES CORRIDAS DOS MELHORES VALORES DE "K" PARA ELABORAÇÃO DOS GRÁFICOS.
#RESSALTA-SE QUE O SCRIPT DEVE SER CORRIDO PARA OS ARQUIVOS QUE FOREM OBTIDOS PELO SNMF.
#SALVE A MELHOR CORRIDA EM FORMATO CSV.

#LEITURA DO RESULTADO DO SNMF.
a_1_q_matrix_nova_k4 <- read.csv(file="D:/Backup_linux/data/snmf/results/snmf_rebuild_sem_swainsoni/a1/arquivo_vcf_a1_para_snmf_r71.4.csv")
head(a_1_q_matrix_nova_k4)

#INCLUSÃO DE INFORMAÇÃO NECESSÁRIA.
a_1_q_matrix_k4_nova_rotulada <- a_1_q_matrix_nova_k4 %>%
  as_tibble()%>%
  mutate(
    Nº=amostras_sem_swainsoni$N,
    AE=amostras_sem_swainsoni$Area_endemismo,
    LAT=amostras_sem_swainsoni$LAT,
    LON=amostras_sem_swainsoni$LON)
a_1_q_matrix_k4_nova_rotulada

#FORMATAÇÃO PARA O GGPLOT.
a_1_q_matrix_k4_nova_long <- a_1_q_matrix_k4_nova_rotulada %>% 
  pivot_longer(cols=starts_with("K"), names_to="pop", values_to="q")
a_1_q_matrix_k4_nova_long

#ORDENAÇÃO DOS QUOFICIENTES DE ANCESTRALIDADE DA MATRIZ DE DADOS, AGRUPANDO OS DADOS POR PROXIMIDADE DOS VALORES.
a_1_q_matrix_k4_nova_ordenado <- a_1_q_matrix_k4_nova_long %>%
  group_by(Nº) %>%
  mutate(likely_assignment=pop[which.max(q)],
         assignment_prob=max(q)) %>%
  arrange(likely_assignment, assignment_prob)%>%
  ungroup()%>%
  mutate(Nº=forcats::fct_inorder(factor(Nº)))
a_1_q_matrix_k4_nova_ordenado

#CORES PARA O GRÁFICO DO SNMF.
my.colors_a1_new_k4 <- c("gray80","gray58","gray29","gray7")

#GRÁFICOS DO SNMF.
plot_snmf_a1_nova_k4 <- a_1_q_matrix_k4_nova_ordenado %>% 
  ggplot()+
  ggtitle("sNMF K=4, alpha=1")+
  geom_col(aes(x=Nº,y=q,fill=pop))+
  scale_fill_manual(values=my.colors_a1_new_k4, labels=c("K1","K2","K3","K4"))+
  labs(fill="K")+
  xlab("Indivíduos")+
  ylab("Coeficiente de ancestralidade")+
  theme_minimal()+
  theme(panel.spacing.x = unit(.7, "lines"),
        strip.background=element_rect(fill="transparent",color="black"),
        panel.background = element_blank(),
        panel.grid=element_blank()
  )
plot_snmf_a1_nova_k4
