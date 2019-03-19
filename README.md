# WGCNA piepine 

## 基本概念
加权基因共表达网络分析 (WGCNA, Weighted correlation network analysis)是用来描述不同样品之间
基因关联模式的系统生物学方法，可以用来鉴定高度协同变化的基因集, 并根据基因集的内连性和基因集
与表型之间的关联鉴定候补生物标记基因或治疗靶点。
相比于只关注差异表达的基因，WGCNA利用数千或近万个变化最大的基因或全部基因的信息识别感兴趣的
基因集，并与表型进行显著性关联分析。一是充分利用了信息，二是把数千个基因与表型的关联转换为数
个基因集与表型的关联，免去了多重假设检验校正的问题。

## 基本原理
从方法上来讲，WGCNA分为表达量聚类分析和表型关联两部分，主要包括基因之间相关系数计算、基因模块
的确定、共表达网络、模块与性状关联四个步骤。

第一步计算任意两个基因之间的相关系数（Person Coefficient）。为了衡量两个基因是否具有相似表达模式，
一般需要设置阈值来筛选，高于阈值的则认为是相似的。但是这样如果将阈值设为0.8，那么很难说明0.8和0.79
两个是有显著差别的。因此，WGCNA分析时采用相关系数加权值，即对基因相关系数取N次幂，使得网络中的基因
之间的连接服从无尺度网络分布(scale-freenetworks)，这种算法更具生物学意义。

第二步通过基因之间的相关系数构建分层聚类树，聚类树的不同分支代表不同的基因模块，不同颜色代表不同的
模块。基于基因的加权相关系数，将基因按照表达模式进行分类，将模式相似的基因归为一个模块。这样就可以
将几万个基因通过基因表达模式被分成了几十个模块，是一个提取归纳信息的过程。

理解WGCNA，需要先理解下面几个术语和它们在WGCNA中的定义。

- 共表达网络：定义为加权基因网络。点代表基因，边代表基因表达相关性。加权是指对相关性值进行冥次
运算 (冥次的值也就是软阈值 (power, pickSoftThreshold这个函数所做的就是确定合适的power))。无向网
络的边属性计算方式为 abs(cor(genex, geney)) ^ power；有向网络的边属性计算方式为 
(1+cor(genex, geney)/2) ^ power; sign hybrid的边属性计算方式为cor(genex, geney)^power if cor>0 else 0。
这种处理方式强化了强相关，弱化了弱相关或负相关，使得相关性数值更符合无标度网络特征，更具有生物意义。
如果没有合适的power，一般是由于部分样品与其它样品因为某种原因差别太大导致的，可根据具体问题移除部分
样品或查看后面的经验值。

- Module(模块)：高度內连的基因集。在无向网络中，模块内是高度相关的基因。在有向网络中，模块内是
高度正相关的基因。把基因聚类成模块后，可以对每个模块进行三个层次的分析：1. 功能富集分析查看其
功能特征是否与研究目的相符；2. 模块与性状进行关联分析，找出与关注性状相关度最高的模块；3. 模块与
样本进行关联分析，找到样品特异高表达的模块。

- Connectivity (连接度)：类似于网络中 “度” (degree)的概念。每个基因的连接度是与其相连的基因的边属性之和。

- Module eigengene E: 给定模型的第一主成分，代表整个模型的基因表达谱。这个是个很巧妙的梳理，我们之前
讲过PCA分析的降维作用，之前主要是拿来做可视化，现在用到这个地方，很好的用一个向量代替了一个矩阵，方便
后期计算。(降维除了PCA，还可以看看tSNE)

- Intramodular connectivity: 给定基因与给定模型内其他基因的关联度，判断基因所属关系。

- Module membership: 给定基因表达谱与给定模型的eigengene的相关性。

- Hub gene: 关键基因 (连接度最多或连接多个模块的基因)。

- Adjacency matrix (邻接矩阵)：基因和基因之间的加权相关性值构成的矩阵。

- TOM (Topological overlap matrix)：把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得的新距离矩阵，
这个信息可拿来构建网络或绘制TOM图。

## 分析步骤
*单组数据的关联分析*
### step1 数据格式
1、表达数据文件的格式
即表达量，行为样本，列放对应的表达量
2、性状文件
行为样本，列为性状值
```shell
>head(datExpr)
              MMT00082832 MMT00082847 MMT00082850 MMT00082869 MMT00082877 MMT00082899
      F2_2    0.0310000   0.1290000   0.0467000  0.00991000   0.0291000    -0.00927
      F2_3    0.0106000   0.1130000  -0.0252000  0.03190000   0.0408000    -0.12100
      F2_14  -0.1310000   0.2550000  -0.1230000  0.08800000   0.0892000    -0.11400
      F2_15   0.0882000   0.0790000   0.0002760 -0.04820000   0.0493000    -0.05010
      F2_19   0.2950000   0.1270000  -0.0560000 -0.02890000  -0.0389000     0.00718
      F2_20   0.1161963   0.1180381  -0.1171272 -0.09774204  -0.0745188     0.31857
>head(trait)
       Aneurysm Aortic_cal_M Aortic_cal_L CoronaryArtery_Cal Myocardial_cal
       F2_290       16            0           17                  0              0
       F2_291       16            4            0                  2              4
       F2_292        0            0           11                  0              0
       F2_293        0            0            0                  0            236
       F2_294       12           10            0                  0              0
       F2_295       17            2            0                  0              0
```
### step2 样本聚类
如果有离群样本，那么就需要去掉
```shell
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK 
if (!gsg$allOK)
    {
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0) 
           printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
           datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
    }
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
png("step1.sampleClustering.png",height=800*3,width=800*3,res=72*3)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
pdf(file = "step1.sampleClustering.pdf", width = 12, height = 9);
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
```
![样本聚类图](https://github.com/nanshanjin/WGCNA/blob/master/step1.sampleClustering.png)
### step3 确定power值
soft thresholding power值即软阈值的筛选原则是使构建的网络更符合无标度网络特征
横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
网络越符合无标度特征 (non-scale)

```shell
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
关键就是理解pickSoftThreshold函数及其返回的对象，最佳的beta值就是sft$powerEstimate
![power图](https://github.com/nanshanjin/WGCNA/blob/master/step2.Model-Fit.png)

### step4 构建共表达矩阵
构建共表达矩阵是WGCNA分析的核心，其核心思想是把表达矩阵的基因分类成十几个模块
大体思路：计算基因间的邻接性，根据邻接性计算基因间的相似性，然后推出基因间的相异性系数，
并据此得到基因间的系统聚类树。然后按照混合动态剪切树的标准，设置每个基因模块最少的基因数目为30。

```shell
net = blockwiseModules(datExpr, power = 6,#6
                     corType="pearson",#or bicor
                     TOMType = "unsigned", minModuleSize = 30,#30
                     reassignThreshold = 0, mergeCutHeight = 0.25,
                     numericLabels = TRUE, pamRespectsDendro = FALSE,
                     saveTOMs = TRUE,
                     saveTOMFileBase = "TOM", 
                     verbose = 3)
```
然后就是分类模块的可视化
这里用不同的颜色来代表那些所有的模块，其中灰色默认是无法归类于任何模块的那些基因，
如果灰色模块里面的基因太多，那么前期对表达矩阵挑选基因的步骤可能就不太合适。

```shell
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
          "Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE
          , guideHang = 0.05)

```

![模块图](https://github.com/nanshanjin/WGCNA/blob/master/step2.plotDendroAndColors.png)
### step5 模块和性状的关联分析
通过模块与各种表型的相关系数，可以很清楚的挑选自己感兴趣的模块进行下游分析了。
这个图就是把moduleTraitCor这个矩阵给用热图可视化一下。
```shell
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
       xLabels = names(datTraits),
       yLabels = names(MEs),
       ySymbols = names(MEs),
       colorLabels = FALSE,
       colors = greenWhiteRed(50),
       textMatrix = textMatrix,
       setStdMargins = FALSE,
       cex.text = 0.5,
       zlim = c(-1,1),
       main = paste("Module-trait relationships"))
```
![热图](https://github.com/nanshanjin/WGCNA/blob/master/step3.labeledHeatmap.png)
### step6 挑选感兴趣的的模块内的基因进行分析
```shell
############################################################
################## select the trait that your interests 
############################################################
# Define variable weight containing the weight column of datTrait
#这里选取了weight_g性状
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

#####module选择对应显著的module
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
pdf(file = "step3.verboseScatterplot.pdf", width = 12, height = 9);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
```
![散点图](https://github.com/nanshanjin/WGCNA/blob/master/step3.verboseScatterplot.png)
### step7 选取指定模块基因名

```shell
# Select module
module = "brown";
# Select module probes 选择模块对应的基因
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
```
### step8 模块导出
```shell
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Select module
module = "brown";
# Select module probes 选择模块对应的基因
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read导出到 VisANT
vis = exportNetworkToVisANT(modTOM,
file = paste("VisANTInput-", module, ".txt", sep=""),
                weighted = TRUE,
                 threshold = 0,) 
#####如果模块包含的基因太多，网络太复杂，还可以进行筛选
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
      file = paste("VisANTInput-", module, "-top30.txt", sep=""),
      weighted = TRUE,
      threshold = 0,
)
```
导出到cytoscape
```shell
cyt = exportNetworkToCytoscape(
        modTOM,
        edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
        nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
        weighted = TRUE,
        threshold = 0.02,
        nodeNames = modProbes, 
        nodeAttr = moduleColors[inModule]
        );
```
