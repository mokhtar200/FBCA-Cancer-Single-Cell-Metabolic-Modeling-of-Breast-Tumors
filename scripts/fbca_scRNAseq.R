# ===============================
# 1. Load Libraries
# ===============================
library(Seurat)
library(dplyr)
library(ggplot2)
library(GSVA)
library(Matrix)
library(patchwork)

set.seed(123)

# ===============================
# 2. Load scRNA-seq Data
# ===============================
counts <- read.csv(
  "data/counts_matrix.csv",
  row.names = 1
)

seu <- CreateSeuratObject(
  counts = counts,
  project = "FBCA_Cancer",
  min.cells = 3,
  min.features = 200
)

# ===============================
# 3. Quality Control
# ===============================
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- subset(
  seu,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 6000 &
    percent.mt < 10
)

# ===============================
# 4. Normalization & Clustering
# ===============================
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(seu, dims = 1:20)

DimPlot(seu, label = TRUE)

# ===============================
# 5. Define Metabolic Gene Sets
# ===============================
metabolic_sets <- list(
  Glycolysis = c("HK2","PFKP","ALDOA","GAPDH","ENO1","PKM","LDHA"),
  OxPhos = c("NDUFA1","NDUFB8","COX5A","ATP5F1"),
  TCA = c("CS","ACO2","IDH3A","SDHA","MDH2"),
  Glutamine = c("GLS","GLUD1","SLC1A5"),
  FA_Synthesis = c("ACACA","FASN","SCD")
)

# ===============================
# 6. Single-Cell Metabolic Scoring
# ===============================
expr <- as.matrix(seu@assays$RNA@data)

metabolic_scores <- gsva(
  expr,
  metabolic_sets,
  method = "ssgsea"
)

seu <- AddMetaData(seu, t(metabolic_scores))

FeaturePlot(
  seu,
  features = c("Glycolysis","OxPhos","TCA"),
  ncol = 3
)

# ===============================
# 7. Define Metabolic States
# ===============================
seu$Metabolic_State <- ifelse(
  seu$Glycolysis >
    quantile(seu$Glycolysis, 0.75),
  "Glycolytic",
  "OxPhos-like"
)

DimPlot(seu, group.by = "Metabolic_State")

# ===============================
# 8. FBCA Grid Construction
# ===============================
n_cells <- ncol(seu)
grid_size <- ceiling(sqrt(n_cells))

fbca_grid <- data.frame(
  cell = colnames(seu),
  x = sample(1:grid_size, n_cells, replace = TRUE),
  y = sample(1:grid_size, n_cells, replace = TRUE),
  Glycolysis = seu$Glycolysis,
  OxPhos = seu$OxPhos
)

fbca_grid$ATP <- fbca_grid$Glycolysis - fbca_grid$OxPhos

# ===============================
# 9. FBCA Transition Rules
# ===============================
update_state <- function(atp){
  if(atp > 0.5){
    return("Proliferative")
  } else if(atp < -0.3){
    return("Dormant")
  } else {
    return("Survival")
  }
}

# ===============================
# 10. FBCA Simulation
# ===============================
time_steps <- 15
history <- list()

for(t in 1:time_steps){
  fbca_grid$State <- sapply(
    fbca_grid$ATP,
    update_state
  )
  fbca_grid$Time <- t
  history[[t]] <- fbca_grid
}

fbca_results <- bind_rows(history)

write.csv(
  fbca_results,
  "results/fbca_simulation.csv",
  row.names = FALSE
)

# ===============================
# 11. Visualization
# ===============================
p <- ggplot(
  fbca_results,
  aes(x = x, y = y, color = State)
) +
  geom_point(size = 2) +
  facet_wrap(~Time) +
  theme_minimal() +
  labs(
    title = "FBCA Simulation of Tumor Metabolic States",
    color = "Cell State"
  )

ggsave(
  "results/figures/fbca_simulation.png",
  p,
  width = 10,
  height = 8
)
