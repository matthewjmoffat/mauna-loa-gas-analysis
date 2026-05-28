## <<===============>> EXPLORATORY DATA ANALYSIS (EDA) <<===============>> ##

# Import required libraries
library(corrplot)
library(ggplot2)
library(mice)
library(VIM)
library(naniar)
library(GGally)
library(psych)

# ====> IMPORTING <====
dataset <- read.csv("/Users/matthewmoffat/Documents/Durham/Epiphany term/DEVUL/Assignments/Assignment 2/MaunaLoa_miss.csv")
names(dataset)

# ============> POINT 1: DISTRIBUTIONS <============

# Histograms of all variables
par(mfrow = c(2, 3))
for (colname in names(dataset)) {
  if (is.numeric(dataset[[colname]])) {
    hist(dataset[[colname]], main = paste("Histogram of", colname),
         col = "skyblue",
         ylab = "Frequency",
         xlab = paste(colname, "value"))
  } else {
    print(paste("Skipping non-numeric column:", colname))
  }
}
par(mfrow = c(1,1))

# Boxplots of all variables
par(mfrow = c(2, 3))
for (colname in names(dataset)) {
  if (is.numeric(dataset[[colname]])) {
    boxplot(dataset[[colname]], main = paste("Boxplot of", colname),
            col = "orange",
            ylab = paste(colname, "count"),
            xlab = paste(colname))
  } else {
    print(paste("Skipping non-numeric column:", colname))
  }
}
par(mfrow = c(1,1))

# Correlation heatmap
cor_matrix <- cor(dataset[, sapply(dataset, is.numeric)], use = "complete.obs")
par(mfrow = c(1,1), mar = c(5,4,6,2) + 0.2)
corrplot(cor_matrix, method = "color",
         diag = FALSE,
         type = "upper",
         tl.col = "black",
         tl.cex = 1.5,
         number.cex = 1.5,
         tl.srt = 45,
         addCoef.col = "black",
         col = colorRampPalette(c("orange", "white", "royalblue"))(150),
         title = "Correlation Heatmap")

# Line plots
par(mfrow = c(2, 3))
for (colname in names(dataset)) {
  if (is.numeric(dataset[[colname]])) {
    plot(dataset[[colname]], type = "l",
         main = paste("Line Plot of", colname),
         col = "navy",
         ylab = colname,
         xlab = "Index")
  } else {
    print(paste("Skipping non-numeric column:", colname))
  }
}
par(mfrow = c(1,1))

# Density + scatter matrix
pairs.panels(dataset[2:6],
             method = "pearson",
             hist.col = "#00AFBB",
             density = TRUE,
             ellipses = TRUE)

# ========> POINT 2: EXTREME VALUES AND OUTLIERS <========

# Function to detect outliers using IQR
detect_outliers <- function(column) {
  Q1 <- quantile(column, 0.25, na.rm = TRUE)
  Q3 <- quantile(column, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  lower_bound <- Q1 - 1.5 * IQR_val
  upper_bound <- Q3 + 1.5 * IQR_val
  return(column[column < lower_bound | column > upper_bound])
}

# Apply function to each numeric variable
outliers_list <- lapply(dataset[, c("CO", "CO2", "Methane", "NitrousOx", "CFC11")],
                        detect_outliers)
outliers_list

# Remove the anomalous 2013 data (rows 131-135)
data_to_remove <- c(131:135)
cleaned_data <- dataset[-data_to_remove, ]
summary(cleaned_data)

# Check if missingness is associated with extreme values
par(mfrow = c(1,1))
boxplot(CO2 ~ is.na(CO), data = cleaned_data)
par(mfrow = c(1,1))
boxplot(CO ~ is.na(CO2), data = cleaned_data)

# ========> POINT 3: DEALING WITH MISSING VALUES <========

# Visualise missing values
par(mfrow = c(1,1))
miss_plot <- gg_miss_var(cleaned_data)
print(miss_plot +
  ggtitle("Missing Value Count by Variable") +
  xlab("Number of Missing Values") +
  ylab("Variables"))

# Visualise missing values as proportions
par(mfrow = c(1,1))
miss_perc_plot <- vis_miss(cleaned_data)
print(miss_perc_plot +
  ggtitle("Missing Values by Proportion") +
  xlab("Proportion of Missing Values") +
  ylab("Observation"))

# Missing data table
pattern_plot <- md.pattern(cleaned_data)
md.pairs(cleaned_data)

# Proportion of missing values using VIM
par(mfrow = c(1,1))
aggr(cleaned_data, col=mdc(1:2), numbers=TRUE, sortVars=TRUE,
     labels=names(cleaned_data),
     cex.axis=.7, gap=3, ylab=c("Proportion of missingness", "Missingness Pattern"),
     prop = TRUE)

# Margin plot
par(mfrow = c(1,1))
marginplot(cleaned_data[, c("CO", "CO2")],
           col = mdc(1:2),
           cex.numbers = 1.2, pch = 19,
           main = "Margin Plot of Variables Displaying Missingness: CO and CO2")

# Contingency table
table(is.na(cleaned_data$CO), is.na(cleaned_data$CO2))

# MCAR Test (Little's)
mcar_test <- mcar_test(cleaned_data[, c("CO", "CO2")])
print("Little's MCAR test:")
print(mcar_test)

# Fisher's Test
missing_CO  <- as.numeric(is.na(cleaned_data$CO))
missing_CO2 <- as.numeric(is.na(cleaned_data$CO2))
fisher.test(missing_CO, missing_CO2)

# ----------> Multiple Imputation via PMM <----------
mice_dataset <- mice(cleaned_data, m=10, seed=123, method = "pmm")
summary(mice_dataset)

# Check and adjust predictor matrix
pred <- mice_dataset$predictorMatrix
print(pred)
pred["CO",  "Date"] <- 0
pred["CO2", "Date"] <- 0

# Impute with adjusted predictor matrix (6 columns: Date, CO, CO2, Methane, NitrousOx, CFC11)
mice_dataset_new <- mice(cleaned_data, m=10, seed=123,
                         method = c("", "pmm", "pmm", "", "", ""),
                         predictorMatrix = pred)
summary(mice_dataset_new)

# Visualise the 10 iterations
par(mfrow = c(1,1))
xyplot(mice_dataset_new, CO2 ~ CO | .imp, pch = 20, cex = 1.4,
       main = "Scatterplot showing for 10 Iterations of Imputed Values")
par(mfrow = c(1,1))
densityplot(mice_dataset_new,
            main = "Density Plots for 10 Iterations of Imputation via PMM")

# Fit linear regression models to imputed datasets
model_final <- with(mice_dataset_new, lm(Methane ~ CO + CO2 + CFC11 + NitrousOx))
summary(pool(model_final))


## <<===============>> DIMENSION REDUCTION (PCA) <<===============>> ##

library(factoextra)
library(ggfortify)
library(dplyr)
library(gridExtra)

# Define completed imputed dataset
imputed_data <- complete(mice_dataset_new, action = 1)
str(imputed_data)

# Ensure Date is in Date format
imputed_data$Date <- as.Date(imputed_data$Date)

# Engineer YearRangethird column
imputed_data <- imputed_data %>%
  mutate(Year = format(Date, "%Y"),
         YearRangethird = case_when(
           Year >= 2000 & Year <= 2005 ~ "2000-2005",
           Year >= 2006 & Year <= 2012 ~ "2006-2012",
           Year >= 2013 & Year <= 2019 ~ "2013-2019",
           TRUE ~ "Other"
         ))

imputed_data$YearRangethird <- factor(imputed_data$YearRangethird,
                                       levels = c("2000-2005", "2006-2012", "2013-2019"))

# ============> PCA <============

dat    <- imputed_data[, 2:6]
pr.out <- prcomp(dat, scale = TRUE)
names(pr.out)
summary(pr.out)
pr.out$rotation

# Variance explained
pr.var <- pr.out$sdev^2
pr.var
sum(pr.var)

# Scree plot
par(mfrow = c(1,1))
fviz_screeplot(pr.out, addlabels = TRUE, main = "Scree Plot of Principal Components")

# Biplot
par(mfrow = c(1,1))
fviz_pca_biplot(pr.out, axes = c(1,2), repel = TRUE, col.ind = 'navy')

# Biplot grouped by year range
par(mfrow = c(1,1))
fviz_pca_biplot(pr.out,
                col.ind = imputed_data$YearRangethird, palette = "jco",
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "YearRange",
                title = "PCA Biplot Grouped by Year Range")

# Autoplot
par(mfrow = c(1,1))
autoplot(pr.out, data = imputed_data, colour = "YearRangethird",
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)

# Variable contributions to PC1 and PC2
par(mfrow = c(1,1))
plot1 <- fviz_contrib(pr.out, choice = 'var', axes = 1, top = 5,
                       title = "Contribution of Original Variables to Dim-1")
plot2 <- fviz_contrib(pr.out, choice = 'var', axes = 2, top = 5,
                       title = "Contribution of Original Variables to Dim-2")
grid.arrange(plot1, plot2, nrow = 1)

# Kaiser's Rule
pr.var <- pr.out$sdev^2
print(pr.var)
pr.out$rotation


## <<===============>> CLUSTERING <<===============>> ##

library(cluster)
library(mclust)
library(factoextra)
library(fpc)
library(GGally)

# Retain first 2 PCA dimensions
pca_data <- pr.out$x[, 1:2]

# <<============ K-MEANS ============>>
set.seed(6543654)

# Silhouette and WSS to find optimal K
par(mfrow = c(1,1))
fviz_nbclust(pca_data, kmeans, method = "silhouette") + labs(title = "SIL/kmeans")
par(mfrow = c(1,1))
fviz_nbclust(pca_data, kmeans, method = "wss") + labs(title = "WSS/kmeans")

# Apply K-Means with k=3
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)

# Silhouette scores
par(mfrow = c(1,1))
sil_score_kmean <- silhouette(kmeans_result$cluster, dist(pca_data))
fviz_silhouette(sil_score_kmean) + labs(title = "Silhouette Score for K-means Clustering")
mean(sil_score_kmean[, 3])

# WCSS
wcss <- sum(kmeans_result$tot.withinss)
print(paste("WCSS for K-means clustering:", wcss))

# Silhouette width across K values
par(mfrow = c(1,1))
sil.width <- rep(NA, 6)
for (i in 2:6) {
  kmeans_temp <- kmeans(pca_data, centers = i, nstart = 25)
  sil <- silhouette(kmeans_temp$cluster, dist(pca_data))
  sil.width[i] <- mean(sil[, 3])
}
plot(2:6, sil.width[2:6], type = "b", pch = 19, col = "blue",
     xlab = "Number of Clusters (K)", ylab = "Average Silhouette Width",
     main = "Silhouette Analysis for Optimal K")
abline(v = which.max(sil.width), col = "red", lty = 2)

par(mfrow = c(1,1))
fviz_nbclust(pca_data, kmeans, method = "silhouette") +
  labs(title = "Silhouette Method for K-means Clustering")

# Visualise K-Means clusters
par(mfrow = c(1,1))
fviz_cluster(kmeans_result, data = pr.out$x[, 1:2],
             ellipse.type = "convex",
             palette = "jco",
             legend.title = "K-Means\nCluster",
             main = "K-Means Clustering in PCA Space")

# Overlay on PCA biplot
par(mfrow = c(1,1))
fviz_pca_biplot(pr.out,
                col.ind = factor(kmeans_result$cluster),
                palette = "jco",
                addEllipses = TRUE,
                label = "var",
                col.var = "black",
                repel = TRUE,
                legend.title = "K-Means\nCluster",
                title = "K-Means Clusters overlaid on PCA Biplot")

# <<============ AGGLOMERATIVE HIERARCHICAL CLUSTERING ============>>

diss      <- dist(pca_data)
agg_clust <- hclust(diss, method = "ward.D2")

# Dendrogram
par(mfrow = c(1,1))
plot(agg_clust, main = "Agglomerative Clustering Dendrogram",
     xlab = "", ylab = "Height (Dissimilarity)", hang = -1, label = FALSE)
K <- 3
rect.hclust(agg_clust, k = K, border = 2:6)

# Silhouette across K values for hierarchical
par(mfrow = c(1,1))
dist_mat <- dist(pca_data)
hc <- hclust(dist_mat, method = "ward.D2")

sil.width.hclust <- rep(NA, 6)
for (k in 2:6) {
  cluster_cut <- cutree(hc, k = k)
  sil <- silhouette(cluster_cut, dist_mat)
  sil.width.hclust[k] <- summary(sil)$avg.width
}
plot(2:6, sil.width.hclust[2:6], type = "b", pch = 19, col = "darkgreen",
     xlab = "Number of Clusters (K)", ylab = "Average Silhouette Width",
     main = "Silhouette Analysis for Hierarchical Clustering")
abline(v = which.max(sil.width.hclust), col = "red", lty = 2)

# Gap statistic
par(mfrow = c(1,1))
gap_stat_hc <- clusGap(pca_data, FUNcluster = hcut, K.max = 6, B = 100)
print(gap_stat_hc)
fviz_gap_stat(gap_stat_hc)

# Cut tree at optimal K
agg_clusters <- cutree(agg_clust, k = K)

# Silhouette analysis
par(mfrow = c(1,1))
sil <- silhouette(agg_clusters, diss)
fviz_silhouette(sil) + ggtitle("Silhouette Plot (Agglomerative)")
mean(sil[, 3])

# Cophenetic correlation
coph_diss <- cophenetic(agg_clust)
cor(diss, coph_diss)

# Agglomerative coefficient
agnes_result <- agnes(pca_data, method = "ward")
agnes_result$ac

# Overlay on PCA biplot
par(mfrow = c(1,1))
fviz_pca_biplot(pr.out,
                col.ind = factor(agg_clusters),
                palette = "jco",
                addEllipses = TRUE,
                label = "var",
                col.var = "black",
                repel = TRUE,
                legend.title = "Hierarchical\nCluster",
                title = "Hierarchical (Agglomerative) Clusters overlaid on PCA Biplot")

# Visualise hierarchical clusters in PCA space
par(mfrow = c(1,1))
fviz_cluster(list(data = pca_data, cluster = agg_clusters),
             geom = "point",
             ellipse.type = "convex",
             palette = "jco",
             legend.title = "Hierarchical\nCluster",
             main = "Hierarchical (Agglomerative) Clustering in PCA Space using Ward's Method")

# <<====> PAIRWISE SCATTERPLOT MATRICES <====>>

# K-Means clusters
par(mfrow = c(1,1))
ggpairs(data.frame(imputed_data[, 2:6]),
        aes(colour = factor(kmeans_result$cluster)),
        corr = FALSE,
        title = "Pairwise Scatterplot Matrix of Clusters (K-means and Hierarchical)")

# Hierarchical clusters
par(mfrow = c(1,1))
ggpairs(data.frame(imputed_data[, 2:6]),
        aes(colour = factor(agg_clusters))) +
  ggtitle("Pairwise Scatterplot Matrix (Dendrogram Clusters)")

