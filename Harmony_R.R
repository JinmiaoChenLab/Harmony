library(ggplot2)
library(tools)
library(flowCore)
library(FlowSOM)
library(dplyr)
library(smallvis)
library(Rtsne)
library('Biobase')
library("gplots")
library(Seurat)
install.packages("reticulate")
library(reticulate)

raw_data = "03_PB22903.txt"
r_data <- read.table(raw_data, sep="\t", header=T,row.names = NULL)
r_data <- as.matrix(r_data)
# Check data and data column names -- for this script to work, the first row must be the column names
fcsfilename <- paste(raw_data,".fcs")

# Create FCS file metadata - column names with descriptions
metadata <- data.frame(name=dimnames(r_data)[[2]],
                       desc=paste('this is column',dimnames(r_data)[[2]],'from your CSV'))
# Create flowframe with tSNE data
data.ff <- new("flowFrame",
               exprs=r_data,
               parameters=AnnotatedDataFrame(metadata))
write.FCS(data.ff, fcsfilename)

file <- paste(raw_data,".fcs")
data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))


# select protein marker columns to use for clustering
#marker_cols <- c(4, 80, 92, 96, 100, 104, 112)
marker_cols <- 6:12
# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)

asinh_scale <- 5
data <- data[, marker_cols]

summary(data)

# create flowFrame object (required input format for FlowSOM)

data_FlowSOM <- flowCore::flowFrame(data)

###################
### RUN FLOWSOM ###
###################

# set seed for reproducibility

set.seed(1234)

# run FlowSOM (initial steps prior to meta-clustering)

out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out)
out <- FlowSOM::BuildMST(out)

# optional visualization

FlowSOM::PlotStars(out)

# extract cluster labels (pre meta-clustering) from output object

labels_pre <- out$map$mapping[, 1]

# specify final number of clusters for meta-clustering (can also be selected 
# automatically, but this often does not perform well)

k <- 8

# run meta-clustering

# note: In the current version of FlowSOM, the meta-clustering function 
# FlowSOM::metaClustering_consensus() does not pass along the seed argument 
# correctly, so results are not reproducible. We use the internal function 
# ConsensusClusterPlus::ConsensusClusterPlus() to get around this. However, this
# will be fixed in the next update of FlowSOM (version 1.5); then the following 
# (simpler) code can be used instead:
#seed <- 1234
#out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)

seed <- 1234
out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, seed = seed)
out <- out[[k]]$consensusClass


# extract cluster labels from output object

labels <- out[labels_pre]

# summary of cluster sizes and number of clusters

table(labels)
length(table(labels))

# save cluster labels

res <- data.frame(cluster = labels)
r_data$cluster <- res$cluster

write.table(r_data, file = paste0(raw_data,"_cluster.txt"), 
            row.names = FALSE, quote = FALSE, sep = "\t")

###########################
### RUN VR Spatial Plot ###
###########################

#cluster_data <- paste0(raw_data,"_cluster.txt")
#df = read.table(cluster_data, header = T, sep = "\t")
df <- r_data
p <- ggplot(df, aes(x = Cell.X.Position, y = Cell.Y.Position, color = as.factor(cluster))) +
  geom_point(aes(colour = as.factor(cluster)),
             size = 1) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#F8766D",
                                "2" = "#CD9600",
                                "3" = "#7CAE00",
                                "4" = "#00BE67",
                                "5" = "#00BFC4",
                                "6" = "#00A9FF",
                                "7" = "#C77CFF",
                                "8" = "#FF61CC")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  ggtitle(paste0(raw_data)) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))
p

################
### HEAT MAP ###
################

raw_data = read.table(paste0(raw_data,"_cluster.txt"), header = T, sep = "\t")
#Select all markers and cluster columns
df = raw_data[,c(6,7,8,9,10,11,12,17)]

p = df[df$cluster == "1", ]
w = apply(p,2,mean)

for (i in 2:8){
  p = df[df$cluster == i, ]
  x = apply(p,2,mean)
  w = rbind(w,x)
}

rownames(w) <- w[,"cluster"]
w = w[,c(1:7)]
w = w[,c(1,3,4,5,6,7)]
#w = scale(w)
w = w[,c(6,4,5,2,1,3)]

colnames(w) <- c("CD3", "CD8", "PD1", "CD103", "Ecad", "Glypican")

heatmap.2(w, scale = "col", col = bluered(100), 
          trace = "none", density.info = "none", Rowv=FALSE, Colv=FALSE)


#################
### RUN t-SNE ###
#################

set.seed(1234)
# prepare data for Rtsne (matrix format required)
data_Rtsne <- data
data_Rtsne <- as.matrix(data_Rtsne)

head(data_Rtsne)
dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm; runtime: 2-3 min)

# note initial PCA is not required, since we do not have too many dimensions
# (i.e. not thousands, which may be the case in other domains)

set.seed(1234)
out_Rtsne <- Rtsne(data_Rtsne, pca = FALSE, verbose = TRUE)




###################
### CREATE PLOT ###
###################

# load cluster labels (if not still loaded)
labels <- paste(raw_data,"_cluster.txt")
#data_labels <- read.table(file_labels, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

labels <- r_data[, "cluster"]

# select points used by Rtsne

labels_plot <- labels[ix][!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne
length(data_Rtsne)
# prepare Rtsne output data for plot

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot[, Cluster_type] <- as.factor(labels_plot)

head(data_plot)


# plot 2-dimensional t-SNE projection

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = data_plot$cluster)) + 
  geom_point(aes(colour = as.factor(cluster)),
             size = 1) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#F8766D",
                                "2" = "#CD9600",
                                "3" = "#7CAE00",
                                "4" = "#00BE67",
                                "5" = "#00BFC4",
                                "6" = "#00A9FF",
                                "7" = "#C77CFF",
                                "8" = "#FF61CC")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  ggtitle(paste0("t-SNE: ",raw_data)) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))
p


##############
#### UMAP ####
##############

write.table(data,file = 'raw_umap', sep = "\t")

#Source Python UMAP Script
#Python UMAP saves UMAP output file as "umap_data.txt"
# UMAP parameters can be changed in the python script "UMAP Script.py"
source_python("UMAP Script.py")
umap <- run_umap('raw_umap')
umap_data <- read.table("umap_data.txt", sep="\t", header=T,row.names = NULL)
umap_data$cluster <- res$cluster
umap_data$cluster <- as.factor(umap_data$cluster)

ggplot(umap_data, aes(x = umap.dim.1, y = umap.dim.2, color = umap_data$cluster)) + 
  geom_point(aes(colour = as.factor(cluster)),
             size = 1) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#F8766D",
                                "2" = "#CD9600",
                                "3" = "#7CAE00",
                                "4" = "#00BE67",
                                "5" = "#00BFC4",
                                "6" = "#00A9FF",
                                "7" = "#C77CFF",
                                "8" = "#FF61CC")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  ggtitle(paste0("UMAP: ",raw_data)) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))
p
