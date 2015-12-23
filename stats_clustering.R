tropical = c('darkorange','dodgerblue','hotpink','limegreen','yellow')
palette(tropical)
par(pch=19)

library(devtools)
library(Biobase)
library(dendextend)

con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file = con)
close(con)

bm = bodymap.eset
pdata = pData(bm)
edata = as.data.frame(exprs(bm))
fdata = fData(bm)

edata = edata[rowMeans(edata) > 5000,]
dim(edata)
edata = log2(edata+1)


# H cluster
# distance between every pair samples
dist1 = dist(t(edata))

colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(as.matrix(dist1),col=colramp,Colv = NA,Rowv = NA)

hclust1 = hclust(dist1)
plot(hclust1,hang=-1)

dend = as.dendrogram(hclust1)
dend = color_labels(hclust1,5,1:5)
plot(dend)

labels_colors(dend) = c(rep(1,10),rep(2,9))
plot(dend)

# Kmeans cluster
kmeans1 = kmeans(edata,centers = 3)
names(kmeans1)
matplot(t(kmeans1$centers),col=1:3,type="l",lwd=3)
table(kmeans1$cluster)
kmeans1$cluster[1:10]

newdata = as.matrix(edata)[order(kmeans1$cluster),]
heatmap(newdata,col=colramp,Colv = NA,Rowv = NA)

kmeans2 = kmeans(edata,centers = 3)
table(kmeans1$cluster,kmeans2$cluster)
