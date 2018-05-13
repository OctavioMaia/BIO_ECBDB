#packages
source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("limma")
biocLite("genefilter")
install.packages("pillar")
library(Biobase)
library(GEOquery)
library(genefilter)

############################## 1- LEITURA E PROCESSAMENTO DE DADOS ###################################
# 1.1- Leitura de dados
gds4794 <- getGEO('GDS4794',destdir = ".")

# 1.2- Metadata
"O objeto GDS � composto por metadata (header do ficheiro)
e a tabela com os dados de express�o."

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

# 1.3- Dados de Express�o & Pr�-processamento
"Dados de express�o: matriz num�rica que guarda os valores num�ricos de express�o
(genes nas linhas, amostras nas colunas). Consultado com a fun��o exprs

-Nomes dos genes (features): featureNames
-Nomes das amostras: sampleNames
-Lista de atributos: varMetaData (retorna um dataframe com dois campos: nome e descri��o do atributo)
-Valores de atributos: phenoData
"

#Transforma��o Logar�tmica
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
"Filtra genes cujo desvio padr�o dos valores de express�o � maior do
que duas vezes a mediana dos desvios padr�o de todos os genes"

#remover dados omissos
eset_complete = eset[complete.cases((exprs(eset))),]
exp = exprs(eset_complete)

#desvios padr�o de todos os genes
"Como o desvio padr�o se trata de uma medida de dispers�o em torno da m�dia populacional,
valores muito baixos n�o retratam bem o dom�nio de valores estudados. Assim sendo podem ser
excluidos dos dados de express�o."

sds = rowSds(exp)

#mediana dos desvios padr�o
m = median(sds)

hist(sds, breaks=50, col = "blue")
abline(v=m, col="red", lwd=4, lty=2)
abline(v=m*2, col="green",lwd=4,lty=2)

"
dados de express�o filtrados
  - remo��o de probesets (dados omissos, sem informa��o em todas as colunas)
  - remo��o de genes com dom�nio de valoers muito baixos
"
eset_filtered = eset_complete[sds <= 2*m,]
eset_filtered

#normaliza��o de dados
exp = exprs(eset_filtered)
exprs(eset_filtered) = scale(exp)
eset_filtered

################################## 2- EXPRESS�O DIFERENCIAL ##########################################