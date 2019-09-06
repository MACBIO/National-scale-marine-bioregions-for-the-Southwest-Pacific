#load data
Grid.4km = read.csv(“C:/data.csv”, header = T)

# convert NA to 0
Grid.20km[is.na(Grid.20km)] <- 0

#scale the grid [0,1]
Grid.20km.scaled = apply(Grid.20km, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

# double the importance of depth
Grid.20km.scaled[,’depth’] = Grid.20km.scaled[,’depth’] * 2

# create data frame to write to
ID = data.frame(attr(Grid.20km, “row.names”)-1)
colnames(ID) = c(“cell_id”)
Grid.means = ID
Grid.20km$ID = NULL

# kmeans clustering

# set number of clusters to use with kmeans
k = 5000

# do kmeans clustering
means = kmeans(Grid.20km.scaled, k, nstart = 20, iter.max = 1000, algorithm = “MacQueen”)
Grid.means = cbind(Grid.means, data.frame(means$cluster))
colnames(Grid.means) = c(“CellID”,”kmeans”)

# do heirarchical clustering using kmeans cluster centers

# create of distance matrix
similarity = dist(means$centers)

# hclust analysis
compare = hclust(similarity)
groups = data.frame(cutree(compare, h = .4))
colnames(groups) = c(“h_4”)

# create lookup table to assign hclust clusters to kmeans clusters
Grid.hclust = cbind(groups)
hclustOutput = as.data.frame((seq(1,nrow(groups))))
hclustOutput = cbind(hclustOutput, groups)
colnames(hclustOutput) = c(“kmeans”,”h_4”)

# assign hclust clusters to each CellID
Grid.means = merge(Grid.means, hclustOutput, by = “kmeans”)

# save output
write.csv(Grid.means, “C:/Grid20km_kmeans_5000_h4.csv”) 
