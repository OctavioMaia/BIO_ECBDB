#packages
source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
biocLite("limma")
biocLite("genefilter")
biocLite("DESeq2")
install.packages("pillar")
install.packages("survival")
install.packages("stringi")
biocLite("GSEABase")
biocLite("goseq")
biocLite("GOstats")
biocLite("Category")
biocLite("GEOmetadb")
biocLite("RSQLite")
biocLite("geneLenDataBase")
biocLite("hgu133plus2.db")
biocLite("rafalib")
biocLite("gplots")
library(Biobase)
library(GEOquery)
library(genefilter)
library(limma)
library(DESeq2)
library(GSEABase)
library(annotate)
library(goseq)
library(GOstats)
library(GEOmetadb)
library(hgu133plus2.db)
library(rafalib)
library(RColorBrewer)
library(gplots)

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
eset_filtered = eset_complete[sds >= 2*m,]
eset_filtered

#normalização de dados
exp = exprs(eset_filtered)
exprs(eset_filtered) = scale(exp)
eset_filtered

################################## 2- EXPRESSÃO DIFERENCIAL ##########################################
"Identificar o conjunto de genes que têm níveis diferentes de expressão comparando 
células normais vs cancerígena através de testes estatísticos de hipóteses.

Para a análise de expressão diferencial vamos comparar tecido pulmonar canceroso (SCLC) com tecido
normal (incluindo tecido de pulmão) através de um modelo de variância linear dos genes.
Uma vez que os valores de expressão genética seguem uma distribuição normal, basta realizar múltiplos testes estatísticos
que permitem determinar a correlação entre genes.
No nosso caso de teste os genes estão divididos em dois grandes grupos (tecido pulmonar canceroso e tecido normal)
  - Criar uma matriz que inclui separadamente os coeficientes dos grupos 'SCLC' e 'normal' e extrair a sua diferença
  através de um contraste.
"
# 1 - Expressão diferencial
tt = rowttests(eset_filtered,"disease.state")
tt$p.value
rank = order(tt$p.value)
p20 = rank[1:20]
tt$p.value[p20]
featureNames(eset_filtered[p20])

# 2 - Análise estatistica - package limma
"Ao testar N genes em simultâneo, o número esperado de falsos positivos é muito alto.
Este package utiliza uma abordagem de modelos lineares para analisar os níveis de expressão linear
através de uma matriz design, que indica quais tipo de tecido correspondem a cada amostra. A matriz
de contraste indica quais as comparações que vão ser feitas entre cada amostra.
Para os testes estatisticos é utilizado o modelo de Bayes de forma a moderar o erro estimado de cada
mudança, o que permite obter melhores resultados, especialmente para um número pequeno de amostras."

eset_filtered$disease.state <- factor(eset_filtered$disease.state,levels = unique(eset_filtered$disease.state))
eset_filtered$disease.state
design <- model.matrix(~eset_filtered$disease.state)
colnames(design) = c('SCLC','normal')
design


fit = lmFit(eset_filtered,design)
fit
cont.matrix = makeContrasts('SCLC','normal',levels=design)
cont.matrix
fit2 = contrasts.fit(fit,cont.matrix)
fit2 = eBayes(fit2)
help("topTable")
diff = topTable(fit2,coef = 2,10, adjust.method = 'BH')

"Como podemos verificar, existem 584 genes que foram considerados diferencialmente expressos 
com as condições determinadas."

test = topTable(fit2,coef = 2,1000, adjust.method = 'BH')
dim(diff[which(diff$adj.P.Val<0.05),])[1]


# 3- Análise de enriquecimento
"Genes identificados na análise anterior são comparados com outros conjuntos de genes,
onde cada um destes contém genes com funções biológicas semelhantes, com objetivo de verificar
se nos genes identificados existe enriquecimento estatisticamente significativo de algum/vários conjuntos.

De forma a obter a package correta foi feita uma procura pela plataforma GPL570
(informação disponível no DataSet Browser) no website bioconductor.
  - http://bioconductor.org/packages/release/data/annotation/html/hgu133plus2.db.html
"

annotation(eset_filtered) <- "hgu133plus2"
eset_filtered

filt = nsFilter(eset_filtered, require.entrez = T, remove.dupEntrez = T, var.func = IQR,
                var.cutoff = 0.5, feature.exclude = "^AFFX")

eset_enrich = filt$eset
eset_enrich

affyUniverse = featureNames(eset_enrich)
affyUniverse
entrezUniverse = unlist(mget(affyUniverse,hgu133plus2ENTREZID))
entrezUniverse

fit.enr = lmFit(eset_enrich,design)
cont.matrixenr = makeContrasts('SCLC','normal',levels=design)
fit2.enr = contrasts.fit(fit.enr,cont.matrixenr)
fit2.enr = eBayes(fit2.enr)
diff.enr = topTable(fit2.enr,coef = 2,1000, adjust.method = 'BH')

smPV = which(diff.enr$adj.P.Val < 0.05)
pvalFiltered = eset_enrich[smPV,]


selectedEntrezIds = unlist((mget(featureNames(pvalFiltered),hgu133plus2ENTREZID)))
selectedEntrezIds

pvalFiltered
eset_enrich
"Parâmetros:
  - testDirection: over, uma vez que o total de genes diferencialmente expressos comparado com
  o total de genes é elevado.
  - ontology: BP (Biological process), ontologia ao qual os termos GO pertencem."

params <- new("GOHyperGParams", geneIds=selectedEntrezIds, universeGeneIds=entrezUniverse,
             annotation="hgu133plus2.db",ontology="BP", pvalueCutoff=0.05, testDirection="over")

"Criar parâmetros e correr os testes estatísticos hipergeométricos (grupos de genes do Gene ontology),
genes sobre expressos."

hgOver <- hyperGTest(params)
summary(hgOver)

"
Através da análise de enriquecimento é possivel ver os processos biológicos em que os genes estão envolvidos
(através da coluna term)."

######################################### 3- CLUSTERING ##############################################
"Se um gene desconhecido i é similar em termos de expressão a um gene conhecido j, é muito provável
que estejam envolvidos na mesma via metabólica.
  - Criar matriz de distâncias entre cada par de genes
  - Pares de genes com distâncias pequenas partilham os mesmos padrões de expressão
O Clustering revela genes com padrões de expressão semelhantes, logo potencialmente
relacionados funcionalmente"
diff = topTable(fit2,coef = 2,50, adjust.method = 'BH')

#ordenar genes por menor p-value
#adj.P.Val corresponde ao p.value ajustado para testes múltiplos
rank = which(diff$adj.P.Val < 0.05)
p10 = eset_filtered[rank]

#Criar matriz de distâncias Euclidianas
dEuc <- dist(exprs(p10))
dEuc

# 1 - Hierarchical clustering
"Uma vez computada a matriz de distâncias, o algoritmo de clustering vai servir para agrupar genes.
Neste tipo de clustering, cada amostra é atrubuida ao seu próprio grupo e o algoritmo itera
juntando, em cada iteração, os dois clusters mais semelhantes até haver apenas um grupo grande.

Na variável d, temos a distância entre genes mas não temos a distância entre grupos. A distância entre
grupos vai ser calculada pair-wise com a ajuda da primitiva hclust, que retorna um objeto descrevendo
o algoritmo falado anteriormente. Para visualizar este algoritmo, recorremos à função plot.

Os algoritmos de clustering definidos a seguir permitem, tendo em conta a condição experimental dos
dados de expressão (tecido normal ou canceroso), quais dos genes tem condições semelhantes, quais os
tecidos com dados de expressão semelhantes."

#método complete: distância entre dois clusters é a distância máxima entre dois elementos
cl.hier <- hclust(dEuc)
plot(cl.hier)

#método simple: distância entre dois clusters é a distância mínima entre dois elementos
cl.hier_single <- hclust(dEuc, method='single')
plot(cl.hier)

#método average
cl.hier_average <- hclust(dEuc, method='average')
plot(cl.hier)

"Os plots anteriores não permitem saber se os genes foram ou não separados por estado de doença. No método
seguinte é possível verificar, visualmente, que à primeira vista, parece que os genes foram separados por
estado de doença, com algumas excepções. Com este dendograma é possível verificar a distância entre grupos
a partir da altura definida do lado esquerdo. De forma a definir os clusters, é necessário cortar a àrvore a
uma dada altura de forma a agrupar todas as amostras à distância definida pelo corte."

myplclust(cl.hier, labels=as.character(p10$disease.state), lab.col=as.numeric(p10$disease.state)+2, cex=0.4)

#linha de corte
abline(h=10)
hclusters <- cutree(cl.hier, h=10)
table(true=p10$disease.state[1:length(hclusters)],cluster=hclusters)

#podemos também determinar o número de clusters que queremos
hclusters <- cutree(cl.hier, k= 4)
table(true=p10$disease.state[1:length(hclusters)],cluster=hclusters)


#k-means
"Escolher o número de k clusters e atribuir um gene a cada cluster aleatoriamente. De seguida
calcula os novos centros e a distância de cada gene a esse novo centro. Por fim cada gene é atribuido
ao centro mais próximo. Este processo é repetido até o algoritmo estar estável.
Ao analisar os resultados para 9 clusters, notamos que obtemos resultados semelhantes."
set.seed(1)
km <- kmeans(exprs(p10),centers = 9)

plot(km$cluster, col = as.numeric(p10$disease.state)+2, pch=16)
plot(km$cluster, col = km$cluster, pch=16)
table(true=p10$disease.state[1:length(km$cluster)],cluster=km$cluster)



#heatmap
"A árvore da linha representa os genes enquanto a da coluna representa o tecido. As cores na 
heatmap representam as intensidades (relações) entre genes."

#definir as cores
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))
cols <- palette(brewer.pal(8, "Dark2"))[as.numeric(p10$disease.state)]
head(cbind(colnames(exprs(p10)),cols))

rv <- rowVars(exprs(p10))
idx <- order(-rv)

heatmap.2(exprs(p10)[idx,], labCol=p10$disease.state,
          trace="none", ColSideColors=cols, col=hmcol)

"Analisando o heatmap é visível uma boa separação entre estados de doença (normal e SCLC) tal como nos
métodos anteriores."

####################################### 4- CLASSIFICAÇÂO #############################################

model_knn = train(t(exprs(eset_filtered)), eset_filtered$disease.state, method = "knn", trControl=trainControl("cv", number = 5))
pred_knn = predict(model_knn, t(exprs(eset_filtered)))
mk1=confusionMatrix(pred_knn, eset_filtered$disease.state)
mk1$table; mk1$overal[1]




model_tree = train(t(exprs(eset_filtered)), eset_filtered$disease.state, method = "rpart", trControl=trainControl("cv", number = 5))
pred_tree = predict(model_tree, t(exprs(eset_filtered)))
mt1 = confusionMatrix(pred_tree, eset_filtered$disease.state)
mt1$table; mt1$overal[1]



model_knn2 = train(t(exprs(eset_filtered)), eset_filtered$tissue, method = "knn", trControl=trainControl("cv", number = 5))
pred_knn2 = predict(model_knn2, t(exprs(eset_filtered)))
mk2=confusionMatrix(pred_knn2, eset_filtered$tissue)
mk2$table; mk2$overal[1]



model_tree2 = train(t(exprs(eset_filtered)), eset_filtered$tissue, method = "rpart", trControl=trainControl("cv", number = 5))
pred_tree2 = predict(model_tree2, t(exprs(eset_filtered)))
mt2=confusionMatrix(pred_tree2, eset_filtered$tissue)
mt2$table; mt2$overal[1]


model_svm2 = train(t(exprs(eset_filtered)), eset_filtered$tissue, method = "svmLinear", trControl=trainControl("cv", number = 5))
pred_svm2 = predict(model_svm2, t(exprs(eset_filtered)))
ms2=confusionMatrix(pred_svm2, eset_filtered$tissue)
ms2$table; ms2$overal[1]