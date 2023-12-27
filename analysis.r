# DIVYESH THIRUKONDA

load("gasch_expression_data.RData")
source("kmeans_cluster.R")

set.seed(1)

myVars <- sort(apply(expr_data, 1, var), TRUE)[1:2000]
cluster_data <- expr_data[names(myVars), ]


cluster_results <- kmeans_cluster(cluster_data, 50)
# Part 5 a. It took 60 iterations


classified_clusters <- (table(cluster_results))
print(classified_clusters)
for (i in seq_along(classified_clusters)) {
  cat("Cluster", i, ":", classified_clusters[[i]], "genes\n")
}

# Part 5 b.
# Cluster 1 : 16 genes
# Cluster 2 : 20 genes
# Cluster 3 : 131 genes
# Cluster 4 : 4 genes
# Cluster 5 : 127 genes
# Cluster 6 : 29 genes
# Cluster 7 : 125 genes
# Cluster 8 : 17 genes
# Cluster 9 : 15 genes
# Cluster 10 : 8 genes
# Cluster 11 : 13 genes
# Cluster 12 : 29 genes
# Cluster 13 : 2 genes
# Cluster 14 : 52 genes
# Cluster 15 : 1 genes
# Cluster 16 : 69 genes
# Cluster 17 : 112 genes
# Cluster 18 : 9 genes
# Cluster 19 : 27 genes
# Cluster 20 : 14 genes
# Cluster 21 : 62 genes
# Cluster 22 : 19 genes
# Cluster 23 : 1 genes
# Cluster 24 : 1 genes
# Cluster 25 : 71 genes
# Cluster 26 : 108 genes
# Cluster 27 : 12 genes
# Cluster 28 : 2 genes
# Cluster 29 : 57 genes
# Cluster 30 : 43 genes
# Cluster 31 : 87 genes
# Cluster 32 : 14 genes
# Cluster 33 : 37 genes
# Cluster 34 : 33 genes
# Cluster 35 : 15 genes
# Cluster 36 : 110 genes
# Cluster 37 : 14 genes
# Cluster 38 : 43 genes
# Cluster 39 : 2 genes
# Cluster 40 : 28 genes
# Cluster 41 : 24 genes
# Cluster 42 : 98 genes
# Cluster 43 : 18 genes
# Cluster 44 : 7 genes
# Cluster 45 : 1 genes
# Cluster 46 : 96 genes
# Cluster 47 : 14 genes
# Cluster 48 : 159 genes
# Cluster 49 : 2 genes
# Cluster 50 : 2 genes

# Part 5 c, d, and e.
clusterGenes <- c()
for (cluster in seq_along(classified_clusters[1:10])) {
    clusterMatrix <- cluster_data[cluster_results == cluster, ]
    clusterAvg <- colMeans(clusterMatrix)
    myfilename <- paste0("Cluster", cluster, ".png")
    png(filename=myfilename) #EC
    matplot(t(clusterMatrix),  xlab = "Stress Condition", ylab = "Gene Expression", title = cat("Cluster", cluster, ",", nrow(clusterMatrix), "genes"), type="l")
    lines(clusterAvg, lwd="2")
    dev.off()

    # for part e, choose cluster 2 instead of 1:10 in the for loop declaration
    clusterGenes <- (gene_info[rownames(clusterMatrix), c("description", "gene_function")])
}
print(colnames(expr_data)[100:106])