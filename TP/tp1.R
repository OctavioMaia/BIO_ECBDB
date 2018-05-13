#packages
source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("ALL")
biocLite("GEOquery")
biocLite("limma")
biocLite("genefilter")
install.packages("pillar")
library(Biobase)
library(GEOquery)
library(genefilter)
library(limma)

############################## 1- LEITURA E PROCESSAMENTO DE DADOS ###################################
# 1.1- Leitura de dados
gds4794 <- getGEO('GDS4794',destdir = ".")

# 1.2- Metadata
"O objeto GDS é composto por metadata (header do ficheiro)
e a tabela com os dados de expressão."

meta = Meta(gds4794)
channel_count = meta$channel_count
description = meta$description
feature_count = meta$feature_count
platform = meta$platform
sample_count = meta$sample_count
sample_organism = meta$sample_organism
sample_type = meta$sample_type
title = meta$title
type = meta$type

col = colnames(Table(gds4794))
col

sample = Table(gds4794)[1:10,1:6]
sample

# 1.3- Dados de Expressão & Pré-processamento
"Dados de expressão: matriz numérica que guarda os valores numéricos de expressão
(genes nas linhas, amostras nas colunas). Consultado com a função exprs

-Nomes dos genes (features): featureNames
-Nomes das amostras: sampleNames
-Lista de atributos: varMetaData (retorna um dataframe com dois campos: nome e descrição do atributo)
-Valores de atributos: phenoData
"

#Transformação Logarítmica
eset <- GDS2eSet(gds4794, do.log2 = T)
eset

feature_names = featureNames(eset)
feature_names[1:10]

sample_names = sampleNames(eset)
sample_names[1:10]

varMetadata(eset)

pheno_data = pData(eset)
pheno_data$sample
pheno_data$disease.state
pheno_data$tissue
pheno_data$description

abstract(eset)

eset["121_at","GSM1060771"]
exprs(eset["121_at","GSM1060771"])

#Filtro Flat Pattern
"Filtra genes cujo desvio padrão dos valores de expressão é maior do
que duas vezes a mediana dos desvios padrão de todos os genes"

#remover dados omissos
eset_complete = eset[complete.cases((exprs(eset))),]
exp = exprs(eset_complete)

#desvios padrão de todos os genes
"Como o desvio padrão se trata de uma medida de dispersão em torno da média populacional,
valores muito baixos não retratam bem o domínio de valores estudados. Assim sendo podem ser
excluidos dos dados de expressão."

sds = rowSds(exp)

#mediana dos desvios padrão
m = median(sds)

hist(sds, breaks=50, col = "blue")
abline(v=m, col="red", lwd=4, lty=2)
abline(v=m*2, col="green",lwd=4,lty=2)

"
dados de expressão filtrados
  - remoção de probesets (dados omissos, sem informação em todas as colunas)
  - remoção de genes com domínio de valoers muito baixos
"
eset_filtered = eset_complete[sds <= 2*m,]
eset_filtered

#normalização de dados
exp = exprs(eset_filtered)
exprs(eset_filtered) = scale(exp)
eset_filtered

################################## 2- EXPRESSÃO DIFERENCIAL ##########################################
"Identificar o conjunto de gentes que têm níveis diferentes de expressão comparando 
células normais vs cancerígena através de testes estatísticos de hipóteses.
  - hipóteste nula: as médias dos níveis de transcrição das duas situações são idênticas

Para a análise de expressão diferencial vamos comparar tecido pulmonar canceroso (SCLC) com tecido
normal (incluindo tecido de pulmão)."


eset_filtered$disease.state = factor(eset_filtered$disease.state)
table(eset_filtered$disease.state)


design = model.matrix(~eset_filtered$disease.state)
fit = lmFit(eset_filtered, design)
fit2 = eBayes(fit)
summary(decideTests(fit))
diff = topTable(fit2, coef=2, 10)
diff

library("ALL")
data(ALL)
exp = exprs(ALL)
maximos = apply(exp,1,max)
minimos = apply(exp,1,min)
vl = maximos/minimos > 2
ALLm2 = ALL[vl,]
exp_m2 = exprs(ALLm2)
exprs(ALLm2) = scale(exp_m2)
varMetadata(ALL)
