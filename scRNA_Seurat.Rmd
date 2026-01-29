IGC Basic scRNA Seurat Object Setup Guide
================
Kat Lande
2026-01-29

- [0 Load Packages and Data](#0-load-packages-and-data)
- [1 Setup Individual Objects](#1-setup-individual-objects)
- [2 Merging Samples](#1-merging-samples)
  - [2.1 Simple Merging](#21-simple-merging)
  - [2.2 Harmony Integration](#22-harmony-integration)


# 0 Load Packages
``` r
# Using R-4.4.3
library(Seurat) # v5.0.1.9001
# Matrix version: 1.7-3
# SeuratObject version: 5.1.0
library(ggplot2) # v3.5.2
library(dplyr) # v1.1.4
library(reshape2) # v1.4.4
```


# 1 Setup Individual Objects
Create a vector of input paths and read them all in iteratively to set up individual objects for each sample.
``` r
# To read in multiple samples iteratively, create a named vector of input paths
# e.g.,:
input_paths <- 
  c(Sample1="/path/to/cellranger/sample1/outs/filtered_feature_bc_matrix",
    Sample2="/path/to/cellranger/sample2/outs/filtered_feature_bc_matrix",
    Sample3="/path/to/cellranger/sample3/outs/filtered_feature_bc_matrix")

# create an empty list to store the unfiltered objects
unfilt_obj_list <- list()

# read in each sample and add it to the object list:
for(i in 1:length(input_paths)){
  
  message(names(input_paths)[i])
  cnt <- Read10X(data.dir = input_paths[[i]])
  obj <- CreateSeuratObject(counts = cnt, project = names(input_paths)[i])
  
  # logNormalize an SCtransform each sample individually before normalization and integration
  # this is the basic architecture we
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(object = obj) # log normalize - creates RNA data slot
  obj <- SCTransform(obj, verbose = FALSE) # SCTransform - creates SCT assay
  obj <- PrepSCTFindMarkers(obj) # PrepSCTFindMarkers ONCE - so you don't need to do this downstream
  
  # add a % mitochrondrial read column to the meta data, most commonly one of these two lines will work for mouse/human:
  # you may need to edit the syntax, "^mt-" searches for genes that start with "mt-" and defines them as mitochondrial reads.
  obj[["percent.mt"]]  <- PercentageFeatureSet(obj, pattern = "^mt-") # usually works for mouse
  # obj[["percent.mt"]]  <- PercentageFeatureSet(obj, pattern = "^MT-") # usually works for human
  
  # add each object to the unfiltered object list and name its entry:
  unfilt_obj_list <- append(unfilt_obj_list, obj)
  names(unfilt_obj_list)[length(unfilt_obj_list)] <- names(input_paths)[i]
  
  # free memory and release garbage after each sample:
  rm(obj, cnt)
  gc() 
}

# save the unfiltered object use for future use:
saveRDS(unfilt_obj_list, "/path/to/dir/Unfiltered_Objects_All_Samples.rds")
```


# 2 Merging Samples
Use a simple merge to combine samples for filtering.
```r
# combine samples with a simple merge to look at filtering statistics:
merged <- merge(unfilt_obj_list[[1]], unfilt_obj_list[2:length(unfilt_obj_list)])
merged <- PrepSCTFindMarkers(merged) # PrepSCTFindMarkers ONCE - so you don't need to do this downstream

Idents(merged) <- "orig.ident"
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
```

```r
# define QC cut-offs based on distributions, and filter low quality cells from the object e.g.,:
subset(merged, subset = nFeature_RNA < 5000 &
         nFeature_RNA > 500 &
         nCount_RNA < 25000 &
         nCount_RNA > 750 &
         percent.mt < 5) -> filtered
```

## 2.1 Simple Merging
Cluster the filtered object on simple merging only using the SCT assay. Integration can over-correct data and remove real signatures. If simple merging alone is sufficient, it can be better to avoid integration all together.
```r
DefaultAssay(filtered) <- "SCT"
filtered <- RunPCA(filtered, npcs = 30, features = rownames(filtered))
filtered <- FindNeighbors(filtered, reduction = 'pca', dims = 1:30)
filtered <- FindClusters(filtered, resolution = 1) # increase the resolution to increase ncluster formation
filtered <- RunUMAP(filtered, reduction = 'pca', dims = 1:30)

# Do the clusters look okay? Are the samples overlapping?
DimPlot(filtered, group.by = "orig.ident") # look at the global dimplot 
DimPlot(filtered, split.by = "orig.ident") # look at the dimplot split out

# optional: save the merged object
# saveRDS(filtered, "/path/to/dir/SimpleMerge_Filtered.rds")

```

## 2.2 Harmony Integration
If the samples don't overlap well using a simple merge, or the clusters look otherwise bad, you likely need to properly integrate your samples. Here we use a Harmony integration:
```r
# these two steps were already run in the section above, 
# but make sure to run them if you skipped the simple merge section:
DefaultAssay(Object) <- "SCT"
filtered <- RunPCA(filtered, npcs = 30, features = rownames(filtered))
# start here to integrate with harmony:
harmony <- IntegrateLayers(object = filtered, method = 'HarmonyIntegration',
                           orig.reduction = 'pca',
                           assay = 'SCT',
                           normalization.method = 'SCT',
                           new.reduction = 'harmony')
harmony <- JoinLayers(harmony, assay = "RNA")
harmony <- FindNeighbors(harmony, reduction = 'harmony', dims = 1:30)
harmony <- FindClusters(harmony, resolution = 1) # increase the resolution to increase ncluster formation
harmony <- RunUMAP(harmony, reduction = 'harmony', dims = 1:30)

DimPlot(harmony, group.by = "orig.ident") # look at the global dimplot 
DimPlot(harmony, split.by = "orig.ident") # look at the dimplot split out

# optional: save the integrated object
# saveRDS(harmony, "/path/to/dir/Harmony_Filtered.rds")
```

Pick either the harmony or simple merge object for downstream analysis.










