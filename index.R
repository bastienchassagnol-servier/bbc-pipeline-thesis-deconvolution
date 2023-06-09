## ----html-output-config, include = FALSE, eval=knitr::is_html_output()----------------------------------------------------------------
# bibliography: "`r rbbt::bbt_write_bib('DeCovarT_vignette_bbt.bib', overwrite = TRUE)`
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.fullwidth = TRUE,
  # We use the implementation available in the
  message = FALSE,
  warning = FALSE,
  cache = FALSE, lazy.cache = FALSE,
  echo=TRUE
)
# rmarkdown::render("index.Rmd", output_format = 'all')


## ----pdf-output-config, include = FALSE, eval=knitr::is_latex_output()----------------------------------------------------------------
## knitr::opts_chunk$set(
##   collapse = TRUE,
##   comment = "#>",
##   fig.align = "center",
##   out.width = "90%",
##   message = FALSE,
##   warning = FALSE,
##   echo = FALSE,
##   cache = FALSE, lazy.cache = FALSE,
##   fig.pos = "H"
## )


## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------
options(max.print = "75")
library(ggplot2)
library(dplyr)
library(pkgdown) # generate automated links to R packages and functions
library(kableExtra)
library(flextable)

# corporate packages
library(bbcVerse)
library(bbcSysBio)

# define function to automatically colour some sections, with respect to the output of the document
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else {
    sprintf("<span style='color: %s;'>%s</span>", color, x)
  }
}

# rmarkdown::pandoc_convert(
#   "test.md",
#   to = "latex",
#   output = "out.tex",
#   options = "--standalone"
# )
# render("test.Rmd", output_format = "latex_document")
# 


## ----active-DT-options, echo = FALSE,include = FALSE----------------------------------------------------------------------------------
DT::datatable(matrix())


## ----active-plotly-options, echo = FALSE,include = FALSE, eval=knitr::is_html_output()------------------------------------------------
htmltools::tagList(
  plotly::ggplotly(ggplot(cars) + geom_point(aes(speed, dist))))


## ----theme, eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------
theme_set(theme_bbc(center_titles = TRUE))


## ----download-geo, eval = FALSE-------------------------------------------------------------------------------------------------------
## illumina_online_eset <- import_normalised_data(GEO = "GSE27073")
## 
## # not run, how to import raw data
## raw_illumina <-
##   import_raw_files(microarray_technology = "illumina",
##                    eset = illumina_online_eset,
##                    destdir = destdir)


## ----download-e-mtab, eval = FALSE----------------------------------------------------------------------------------------------------
## illumina_online_eset <- import_normalised_data(GEO = "GSE27073")
## 
## # not run, how to import raw data
## raw_illumina <-
##   import_raw_files(microarray_technology = "illumina",
##                    eset = illumina_online_eset,
##                    destdir = destdir)


## ----import-local-data, eval = FALSE--------------------------------------------------------------------------------------------------
## data_path <- "path_to_your_data"
## 
## count <- read.delim(paste0(data_path, "counts.txt", sep = ""))
## fpkm <- read.delim(paste0(data_path, "fpkm.txt", sep = ""))
## tpm <- read.delim(paste0(data_path, "tpm.txt", sep = ""))
## microarray <- read.delim(paste0(data_path, "microarray.txt", sep = ""))
## 
## meta_data <- read.delim(paste0(data_path, "meta_data.txt", sep = ""))
## gene_annotation <- read.delim(paste0(data_path, "gene_annotation.txt", sep = ""))
## clinical_data <- read.delim(paste0(data_path, "clinical_data.txt", sep = ""))


## ----store-eset, eval = FALSE---------------------------------------------------------------------------------------------------------
## # Put result in an eset object
## eset_microarray <-
##   bmkanalysis::create_eset_object(
##     gene_expression_data = microarray_data$Expression_Intensity_Probes_ByGeneLevel %>%
##       column_to_rownames("GeneID") %>%
##       as.matrix(),
##     phenotype_data = clinicaldetails %>%
##       column_to_rownames("ID"),
##     gene_annotation = NULL
##   )
## 


## ----create-eset-object, eval = FALSE-------------------------------------------------------------------------------------------------
## eset_count <-
##   eset_object(
##     expression_data = count,
##     phenotype_data = clinical_data,
##     feature_data = gene_annotation
##   )
## 


## -------------------------------------------------------------------------------------------------------------------------------------
library(org.Hs.eg.db)
library(AnnotationDbi)


## ----gene annotation, eval = TRUE-----------------------------------------------------------------------------------------------------
feature_annot <- AnnotationDbi::select(org.Hs.eg.db, # database
                                     keys = eset_count %>% rownames(),  # data to use for retrieval
                                     columns = c("SYMBOL", "GENETYPE","GENENAME"), # information to retreive for given data
                                     keytype = "SYMBOL") %>% 
                                     mutate(across(GENETYPE, ~ relevel(factor(GENETYPE), ref = "protein-coding"))) %>% 
                                     arrange(GENETYPE) %>% 
                                     mutate(duplicate = duplicated(SYMBOL)) %>%
                                     filter(duplicate == FALSE) %>%
                                     dplyr::select(-duplicate) %>%
                                     arrange(SYMBOL)
                                      

new_gene_annotation <-  left_join(gene_annotation, feature_annot, by = c(  "GeneID" = "SYMBOL" )) %>% 
                        data.frame()

rownames(new_gene_annotation) <- new_gene_annotation$GeneID 

eset_count <- update_eset_object(eset = eset_count, feature_data = new_gene_annotation)
saveRDS(eset_count,file = "./data/eset_raw.Rds")


## ----estimate-cutoff, include = FALSE, message = FALSE, warning = FALSE---------------------------------------------------------------
# esti_count <- estimate_cutoff_lowcounts(eset_count %>% exprs)
# esti_tpm <- estimate_cutoff_lowcounts_norm(eset_tpm %>% exprs)
esti_count <- 8.1
esti_tpm <- 1


## ----melt-data, eval = TRUE, include = FALSE, message = FALSE, warning = FALSE--------------------------------------------------------
count_melt <- pivot_and_join(eset_count)
tpm_melt <- pivot_and_join(eset_tpm) 


## ----kernel-plot, include = TRUE, message=FALSE, fig.height = 4, fig.width = 8, fig.align = 'center'----------------------------------

kernel_count <- draw_kernel(
  melt_data = count_melt,
  x_label_kernel = "expression",
  id_colname = "samples_id",
  col_color  = "Treatment",
  palette = "dark2",
  threshold = log2(esti_count+1),
  title = "Density plot of raw counts",
  xlab = "log2(Count + 1)",
  log2 = TRUE,
  add_log_value = 1,
)

kernel_tpm <- draw_kernel(
  melt_data = tpm_melt,
  x_label_kernel = "expression",
  id_colname = "samples_id",
  col_color  = "Treatment",
  palette = "dark2",
  threshold = log2(esti_tpm+1),
  title = "Density plot of tpm",
  xlab = "log2(tpm + 1)",
  log2 = TRUE,
  add_log_value = 1
)

ggpubr::ggarrange(kernel_count, kernel_tpm, ncol = 2)



## ----gene-filtering-background, include = F-------------------------------------------------------------------------------------------
NSample <- 6

eset_count_filtered <-
  filter_background(eset_count, esti_count, NSample)
eset_tpm_filtered <-
  filter_background(eset_tpm, esti_tpm, NSample)

common_gene <-
  intersect(rownames(eset_count_filtered), rownames(eset_tpm_filtered))

eset_count_filtered <-
  filter_features_from_values(eset = eset_count_filtered,
                             values = common_gene,
                             column = "GeneID")
eset_tpm_filtered <-
  filter_features_from_values(eset = eset_tpm_filtered,
                             values = common_gene,
                             column = "GeneID")



## ---- eval = TRUE---------------------------------------------------------------------------------------------------------------------
# List of samples to keep
samples_list <- pData(eset_count_filtered) %>% 
  filter(RIN > 6) %>% 
  dplyr::pull(SampleID) 

# Creation of filtered eSet options
eset_count_filtered_RIN <-
  filter_samples_from_values(eset = eset_count_filtered,
                             column = "SampleID",
                             values =   samples_list)

eset_count_filtered_RIN <-
  filter_samples_from_expr(eset = eset_count_filtered, RIN > 6)


## -------------------------------------------------------------------------------------------------------------------------------------
display_table(pData(eset_count_filtered) %>% dplyr::select(., SampleID, Donor, Treatment, Dose, RIN),
                   interactive = TRUE,
                   bar_column = "RIN",
                   highlighted_row = RIN < 6,
                   highlighted_color = "lightblue",
                   nmax_display = 10
                 )


## ----log2transfo, include = F---------------------------------------------------------------------------------------------------------
# counts_log2 <- log2(exprs(eset_count_filtered) + 1)
# 
# # Creation of an eSet with log2(.+1) expressions
# eset_count_filtered_log2 <-
#   create_eset_object(counts_log2,
#                      pData(eset_count_filtered),
#                      fData(eset_count_filtered))

eset_count_filtered_log2 <- transfrom_exprs(eset_count_filtered, 
                                                ~ log2(. + 1))



## ----vsttransfo, include = F----------------------------------------------------------------------------------------------------------
# Creation of a DESeq dds object to compute vst transformation
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = exprs(eset_count_filtered),
  colData = pData(eset_count_filtered),
  design = ~ 1 + Treatment
) # This model is the model used to direct the vst normalization and can include the variable of interest

# Estimation of vst normalized expression matrix
counts_vst <- SummarizedExperiment::assay(dds %>% DESeq2::vst(blind = F)) 

# Creation of an eSet with vst normalized expression
eset_count_filtered_vst <-
  eset_object(expression_data = counts_vst,
              phenotype_data = pData(eset_count_filtered), 
              feature_data = fData(eset_count_filtered))


## ---- fig.width=11, fig.height=6------------------------------------------------------------------------------------------------------
draw_boxplot(
  eset_tpm_filtered,
  col_y = NULL,
  col_x = "SampleID",
  col_color = "SampleID",
  palette = "dark2",
  add_point = "point",
  point_size = 1,
  y_transformation = ~ log2(. + 1),
  facet_by = "Treatment",
  scales_facet = "free") +
  labs(title = "Gene expression per sample\nby Treatment") + 
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



## ----vst-quality, fig.width=10, fig.height=6------------------------------------------------------------------------------------------
# Boxplot pre vst
g1 <- draw_boxplot(
  eset_count_filtered,
  col_y = NULL,
  col_x = "SampleID",
  col_color = "SampleID",
  palette = "dark2",
  add_point = "point",
  point_size = 1,
  y_transformation = ~ log2(. + 1),
  #facet_by = "Treatment",
  hline_threshold = log2(esti_count + 1)) +
  labs(title = "Gene expression per sample before normalization") +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank())



# Boxplot post vst
g2 <- draw_boxplot(
  eset_count_filtered_vst,
  col_y = NULL,
  col_x = "SampleID",
  col_color = "SampleID",
  palette = "dark2",
  add_point = "point",
  point_size = 1,
  #facet_by = "Treatment",
  hline_threshold = log2(esti_count + 1)) +
  labs(title = "Gene expression per sample after normalization") +
  theme(legend.position = "none",
        axis.title.x = element_blank(), 
        axis.text.x = element_blank())

ggpubr::ggarrange(g1, g2, nrow = 2)


## ---- include = F---------------------------------------------------------------------------------------------------------------------
library(sva)


## ----pca------------------------------------------------------------------------------------------------------------------------------
pca_results <- compute_pca(eset_count_filtered_vst)


## ----pca-details, fig.width=5, fig.height=5-------------------------------------------------------------------------------------------
draw_pca_eigenvalues(pca_results)


## ---- fig.width=5, fig.height=5-------------------------------------------------------------------------------------------------------
draw_pca_loadings(pca_results, select.var = list(contrib = 2))


## ----pca-graph1, fig.width=12, fig.height=5-------------------------------------------------------------------------------------------

g1 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "Donor",
  palette = "nejm",
  add_ellipse = TRUE,
  title = "PCA of the samples by Donor"
  ) 

g2 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "PCA of the samples by Treatment",
  shape_manual = c(5, 6, 7, 8)
  ) 


g3 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "RIN",
  palette = "spectral",
  add_ellipse = FALSE,
  title = "PCA of the samples by RIN"
)

ggpubr::ggarrange(g1, g2, g3,
                  ncol = 3,
                  common.legend = FALSE,
                  legend = "bottom")


## ----pca-graph2, fig.width = 7, fig.height = 5----------------------------------------------------------------------------------------
draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "PCA of the samples by Treatment",
  interactive = TRUE, 
  col_hover_info= c("Treatment", "Donor")
)  


## ----pca-graph3, fig.width = 7, fig.height = 5----------------------------------------------------------------------------------------
draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "RIN",
  palette = "spectral",
  add_ellipse = FALSE,
  title = "PCA of the samples by RIN",
  interactive = TRUE
  )  


## ----plsda, include = FALSE-----------------------------------------------------------------------------------------------------------
plsda_results <- compute_plsda(eset_count_filtered_vst, y_outcome_column = "Treatment")

plsda_results_multi <- compute_plsda(eset_count_filtered_vst, y_outcome_column = "Treatment", 
                                     multilevel = "Donor")


## ----plsda-details, fig.width=8, fig.height=5-----------------------------------------------------------------------------------------
draw_plsda_loadings(
  plsda_results,
  eset_count_filtered_vst,
  title = "PLSDA variable representation",
  comp = 1:2,
  gene_label_column = "GeneID",
  rad.in = 0.9,
  cutoff = 0.8,
  ndisplay = 50,
  interactive = FALSE
) + coord_equal() 


## ----plsda-details-2, fig.width=8, fig.height=5---------------------------------------------------------------------------------------
draw_plsda_importance(
  plsda_results,
  eset_count_filtered_vst,
  title = "PLSDA scatter plot for individuals",
  comp = 1,
  ndisplay = 20,
  gene_label_column = "GeneID",
  contrib = 'max',
  method = 'mean'
)


## ----plsda-graph,  fig.width=9, fig.height=5------------------------------------------------------------------------------------------

g1 <- draw_plsda_individuals(
  eset_count_filtered_vst,
  plsda_results,
  add_ellipse = TRUE,
  comp = 1:2,
  palette = "dark2",
  title = "PLSDA of the samples by Treatment"
)  

g2 <- draw_plsda_individuals(
  eset_count_filtered_vst,
  plsda_results_multi,
  add_ellipse = TRUE,
  comp = 1:2,
  palette = "dark2",
  title = "PLSDA of the samples by Treatment \nleveled by Donor"
)  

ggpubr::ggarrange(g1,
                  g2,
                  common.legend = T,
                  legend = "bottom")



## ----plsda-graph-2, fig.width = 8, fig.height = 5-------------------------------------------------------------------------------------

draw_plsda_individuals(
  eset_count_filtered_vst,
  plsda_results,
  add_ellipse = TRUE,
  comp = 1:2,
  palette = "dark2",
  title = "PLSDA of the samples by Treatment", 
  interactive = TRUE
)  



## ----sample-clustering, fig.width=5, fig.height=5-------------------------------------------------------------------------------------
# Create a DESeq object dds to be used for clustering
dds <-
  DESeq2::DESeqDataSetFromMatrix(
    countData = exprs(eset_count_filtered),
    colData = pData(eset_count_filtered_vst),
    design = ~ Treatment + Donor
  )

# Caculation of vst normalized expression according to a provided model
vsd <- DESeq2::vst(dds, blind = FALSE)

# Calculation of an Euclidean distance matrix between samples
sampleDists <- SummarizedExperiment::assay(vsd) %>% t() %>% dist()

# Formatting distance matrix for clustering
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$Donor, sep = "_")
colnames(sampleDistMatrix) <- paste(vsd$Treatment, vsd$Donor, sep = "_")
colors <-
  colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)

# Sample by sample heatmap and clustering
pheatmap::pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  show_rownames = T,
  clustering_distance_cols = sampleDists,
  col = colors)



## ----ComBat, include = FALSE----------------------------------------------------------------------------------------------------------

# Reference model for ComBat including the variables to be preserved (e.g variable of interests)
mod_combat <-
  model.matrix( ~ 1 + Treatment, data = pData(eset_count_filtered_vst))

# Applying ComBat
count_filtered_vst_combat <-
  sva::ComBat(
    dat = eset_count_filtered_vst %>% exprs %>% as.matrix,
    batch = eset_count_filtered_vst %>% pData %>% dplyr::pull(Donor) %>% as.character,
    mod = mod_combat,
    par.prior = TRUE,
    prior.plots = FALSE
  )

# Create normalized eSet
eset_count_filtered_vst_combat <-
  eset_object(
    expression_data = count_filtered_vst_combat,
    phenotype_data = pData(eset_count_filtered),
    feature_data = fData(eset_count_filtered)
  )



## ----ComBat-PCA,  fig.width=11, fig.height=5------------------------------------------------------------------------------------------

# Calulates PCA before and after batch removal
pca_results <- compute_pca(eset_count_filtered_vst)
pca_results_combat <- compute_pca(eset_count_filtered_vst_combat)

# PCA plots
d1 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "PCA of the samples before batch correction",
)  

d2 <- draw_pca_individuals(
  eset_count_filtered_vst_combat,
  pca_results_combat,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "PCA of the samples after batch correction",
)  

ggpubr::ggarrange(d1, d2,
                  ncol = 2 ,
                  common.legend = T,
                  legend = "bottom")


## ----residuals------------------------------------------------------------------------------------------------------------------------

# This model is passed to the linear model and corresponds to the variables one wants to remove the effects
mod_residual <- model.matrix( ~ 1 + Donor, data = pData(eset_count_filtered))

# Vst normalization 
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = exprs(eset_count_filtered),
  colData = pData(eset_count_filtered),
  design = ~ 1 + Treatment # This model is the model used for the normalization vst and can include the variable of interest
) 

dds.vst  <- SummarizedExperiment::assay(dds %>% DESeq2::vst(blind = F))

# Run of Limma linear regression model
dds.limma.fit <- limma::lmFit(dds.vst, design = mod_residual)

# Residuals of previous model corresponds to the normalized matrix
resid <- residuals(dds.limma.fit, y = dds.vst)

# Creation a a normalized eSet
eset_count_filtered_vst_residuals <-
  eset_object(
    expression_data = resid,
    phenotype_data = pData(eset_count_filtered),
    feature_data = fData(eset_count_filtered)
  )



## ----residuals-PCA,  fig.width=11, fig.height=5---------------------------------------------------------------------------------------

pca_results <- compute_pca(eset_count_filtered_vst)
pca_results_residuals <- compute_pca(eset_count_filtered_vst_residuals)

d1 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_results,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "PCA of the samples before batch correction",
)  

d2 <- draw_pca_individuals(
  eset_count_filtered_vst_residuals,
  pca_results_residuals,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "PCA of the samples after batch correction",
)  

ggpubr::ggarrange(d1, d2,
                  ncol = 2 ,
                  common.legend = TRUE,
                  legend = "bottom")


## ----sva, message = FALSE, include = FALSE--------------------------------------------------------------------------------------------
# Full model including all variables of interest (that we know have an effect on the outcome that we want to estimate)
mod_sva <-
  model.matrix( ~ 1 + Treatment, data = pData(eset_count_filtered_vst))

# Null model with only an intercept
mod_null <-  model.matrix( ~ 1, data = pData(eset_count_filtered_vst))

# Estimation of the number of surogate variables
n_sv = num.sv(eset_count_filtered_vst %>% exprs , mod_sva, method = "be")

# Calculation of the surogate variables
svobj = sva(eset_count_filtered_vst %>% exprs, mod_sva, mod_null, n.sv = n_sv)

# svobj$sv contains the surrogate variables that can be used as adjustments in DEA models or as known batch variables in the 2nd approach using Residuals to obtain a batch free eSet


## ----sva-eval, fig.width=11, fig.height=5---------------------------------------------------------------------------------------------
# Creation of a bogus pca_object to plot the surrogate variables
pcs <- svobj$sv
colnames(pcs) <- paste0("Dim.", 1:ncol(pcs))
rownames(pcs) <- colnames(eset_count_filtered_vst)

eig <- data.frame(eigenvalue = rep(0, n_sv))
eig <- eig %>% mutate(var = 100*eigenvalue/sum(eigenvalue)) %>% mutate(cumvar = cumsum(var))
colnames(eig) <- c("eigenvalue", "percentage of variance", "cumulative percentage of variance")
rownames(eig) <- paste0("comp ", 1:nrow(eig))

# Creation of a bogus PCA object
pca_object <- list(ind = list(coord = pcs), eig = eig)

# PCA plot of the SV
d1 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_object,
  col_color = "Treatment",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "Sample plot on the Surrogate Variable spaces",
)

d2 <- draw_pca_individuals(
  eset_count_filtered_vst,
  pca_object,
  col_color = "Donor",
  palette = "dark2",
  add_ellipse = TRUE,
  title = "Sample plot on the Surrogate Variable spaces",
)

ggpubr::ggarrange(d1,
                  d2,
                  ncol = 2 ,
                  common.legend = F,
                  legend = "bottom")


## ----Global-variables-1---------------------------------------------------------------------------------------------------------------
t_pvalue <- 0.05
t_FC <- 1.3


## ----DESeq2-no-cov--------------------------------------------------------------------------------------------------------------------
# Model without covariate
model <- formula(~  Treatment)

# Model with covariates (Treatment preferably last)
model <- formula(~ Donor + Treatment )

# Contrasts
contr_list = list(c("Treatment", "AAA", "Placebo"),
                     c("Treatment", "BBB", "Placebo"),
                     c("Treatment", "CCC", "Placebo"))

# Perform differential expression analysis with DESeq
DEA_res <- dea_deseq(eset_object = eset_count_filtered,
           model = model,
           contr_list = contr_list, 
           feature_colname = "Genes")

# Add the tables of DEGs to the results
DEA_res <- DEA_res %>% mutate(DEgenes = purrr::map(data,
                       ~ .x %>% subset_deg(., "FC", t_FC, "FDR", t_pvalue, order = "FC")))



## ----fixed-effect-model---------------------------------------------------------------------------------------------------------------
# Without covariates
model <- model.matrix( ~ 0 + Treatment, data = eset_count_filtered)
colnames(model) <- c("Placebo", "AAA", "BBB", "CCC")

# With covariates (Treatment preferably first)
model <- model.matrix( ~ 0 + Treatment + Donor, data = eset_count_filtered)
colnames(model) <-
  c("Placebo",
    "AAA",
    "BBB",
    "CCC",
    "D2",
    "D3",
    "D4",
    "D5",
    "D6")

# Contrasts
contr_matrix <- makeContrasts(
  AAA_vs_Placebo = AAA - Placebo,
  BBB_vs_Placebo = BBB - Placebo,
  CCC_vs_Placebo = CCC - Placebo,
  levels = colnames(model)
)

# Perform count data normalization to apply Limma
data_normalized <- dea_limma_normalization(eset_object = eset_count_filtered,
                                           model = model)

# Skip the dea_limma_normalization step if you already have normalized count and want to use them in Limma
# and inout data_norm = NULL in dea_limma

# Perform differential expression analysis with Limma
DEA_res <-
  dea_limma(eset_object = eset_count_filtered,
            data_norm = data_normalized,
            model = model,
            contr_matrix = contr_matrix, 
           feature_colname = "Genes")

# Add the tables of DEGs to the results
DEA_res <- DEA_res %>% dplyr::mutate(DEgenes = purrr::map(data,
                       ~ .x %>% subset_deg(., "FC", t_FC, "FDR", t_pvalue, order = "FC")))



## ----random-effect-RNASeq-------------------------------------------------------------------------------------------------------------
# Model specification
# /!\ do not put random variable here in the model, it is done below

model <- model.matrix( ~ 0 + Treatment , data = eset_count_filtered)

colnames(model) <-
  c("Placebo",
    "AAA",
    "BBB",
    "CCC")

# Contrasts
contr_matrix <- makeContrasts(
  AAA_vs_Placebo = AAA - Placebo,
  BBB_vs_Placebo = BBB - Placebo,
  CCC_vs_Placebo = CCC - Placebo,
  levels = colnames(model)
)

# Perform count data normalization to apply Limma
data_normalized <- dea_limma_normalization(eset_object = eset_count_filtered,
                                           model = model,
                                           random = TRUE,
                                           var_random = "Donor")

# Perform differential expression analysis with Limma
DEA_res <-
  dea_limma(eset_object = eset_count_filtered,
            data_norm = data_normalized,
            model = model,
            contr_matrix = contr_matrix,
            random = TRUE,
            var_random = "Donor", 
            feature_colname = "Genes")

# Add the tables of DEGs to the results
DEA_res <- DEA_res %>% dplyr::mutate(DEgenes = purrr::map(data,
                       ~ .x %>% subset_deg(., "FC", t_FC, "FDR", t_pvalue, order = "FC")))



## ----continuous-vs-categorical, eval = F----------------------------------------------------------------------------------------------
## # Model specification
## pData(eset_count_filtered_vst)$Dose = as.numeric(pData(eset_count_filtered_vst)$Dose)
## 
## model <- model.matrix( ~ 0 + Treatment + Treatment:Dose , data = eset_count_filtered)
## 
## colnames(model) <-
##   c("Placebo",
##     "AAA",
##     "BBB",
##     "CCC",
##     "Placebo_Dose",
##     "AAA_Dose",
##     "BBB_Dose",
##     "CCC_Dose")
## 
## # Contrasts
## coeff_name = "AAA_Dose"
## 
## # Perform count data normalization to apply Limma
## data_normalized <- dea_limma_normalization(eset_object = eset_count_filtered,
##                                            model = model,
##                                            random = FALSE)
## 
## # Perform differential expression analysis with Limma
## DEA_res <-
##   dea_limma(eset_object = eset_count_filtered,
##             data_norm = data_normalized,
##             model = model,
##             coeff_name = coeff_name,
##             feature_colname = "Genes")
## 
## # Add the tables of DEGs to the results
## DEA_res <- DEA_res %>% dplyr::mutate(DEgenes = purrr::map(data,
##                        ~ .x %>% subset_deg(., "FC", t_FC, "FDR", t_pvalue, order = "FC")))
## 


## ----categorical-vs-categorical, eval = F---------------------------------------------------------------------------------------------
## # Model specification
## model <- model.matrix( ~ 0 + Donor + Treatment  , data = eset_count_filtered)
## 
## colnames(model) <-
##   c("Placebo",
##     "AAA",
##     "BBB",
##     "CCC",
##     "Dose")
## 
## # Contrasts
## coeff_name = "Dose"
## 
## # Perform count data normalization to apply Limma
## data_normalized <- dea_limma_normalization(eset_object = eset_count_filtered,
##                                            model = model,
##                                            random = FALSE)
## 
## # Perform differential expression analysis with Limma
## DEA_res <-
##   dea_limma(eset_object = eset_count_filtered,
##             data_norm = data_normalized,
##             model = model,
##             coeff_name = coeff_name,
##             feature_colname = "Genes")
## 
## # Add the tables of DEGs to the results
## DEA_res <- DEA_res %>% dplyr::mutate(DEgenes = purrr::map(data,
##                        ~ .x %>% subset_deg(., "FC", t_FC, "FDR", t_pvalue, order = "FC")))
## 


## ---- results= 'asis', echo = FALSE, eval= TRUE---------------------------------------------------------------------------------------
cat("\n")
purrr::pwalk(DEA_res,
             function(data, contrast,  ...) {
               
               cat(paste0("  \n#### ", contrast, " \n"))
               cat("\n")
               cat(knitr::knit_print(
                 table_deg_significant(
                   result_table = data,
                   foldchange_column = "FC",
                   foldchange_threshold = c(1.3, 1.5, 2),
                   padj_column = "FDR",
                   padj_threshold = c(0.05, 0.1, 0.2),
                   interactive = interactive
                 )))
              })


## ---- fig.width=15, fig.height=5, fig.fullwidth=TRUE----------------------------------------------------------------------------------
pmap(DEA_res, function(data, contrast, ...) {
  draw_pvalue_histogram(
    data,
    "pvalue",
    pi0_estimate = TRUE,
    alpha = 0.05
  )  
  
})  %>% ggpubr::ggarrange(
  plotlist = .,
  nrow = 1,
  labels = paste(DEA_res$contrast),
  common.legend = TRUE,
  legend = "bottom",
  vjust = 1.2
) %>%   print()



## ---- results='asis', fig.width= 8, fig.height=6--------------------------------------------------------------------------------------
purrr::pwalk(DEA_res,
             function(data, contrast, ...) {
               cat(paste0("  \n#### ", contrast , "  \n"))
               draw_volcano(
                 data = data,
                 col_log2fc = "log2FC",
                 col_signif = "pvalue",
                 col_label = "Genes",
                 fc_cutoff = 1.3,
                 signif_cutoff = 0.05,
                 type = c("pv"),
                 labelled = data %>% pull(Genes) %>% head(100),
                 #focused = c("C1QA", "C1QB", "C1QC"),
                 shape = c(17, 1, 16),
                 alpha = c(0.9, 0.4, 0.7),
                 size = c(5, 2, 2),
                 legend = TRUE
                 )   %>% print()
               cat("\n")
               
             })

  


## ---- results='asis', fig.width = 9, fig.height = 7-----------------------------------------------------------------------------------

# Creates a list of all unique pairs of contrasts
list_contrast <- combn(DEA_res$contrast, 2) %>% t()

# Creates a tab with the plot for each contrast comparison
walk2(list_contrast[, 1], list_contrast[, 2],
      function(.x, .y) {
        cat(paste0(" \n### ", .x, "_VS_", .y, " \n"))
        
        draw_concord(
  fulltable,
  col_fc_x = paste("log2FC",.x , sep =  "_"),
  col_fc_y = paste("log2FC",.y , sep =  "_"),
  col_pv_x = paste("FDR",.x , sep =  "_"),
  col_pv_y = paste("FDR",.y , sep =  "_"),
  th_fc = log2(t_FC),
  th_pv = t_pvalue,
  smooth = "lm",
  bisector = TRUE,
  stats = c("rho", "R2", "MAPE"),
  col_label = "Genes", 
  labs_x = .x, 
  labs_y = .y, 
  stat_to_display = "fc",
  labelled = head(fulltable$Genes, 20)
)  %>% print
               cat("\n")
      })


## ---- fig.width= 8, fig.height=7------------------------------------------------------------------------------------------------------
# Select gene list for heatmap
gene_list = DEA_res$DEgenes %>% purrr::map( ~ .x) %>% bind_rows()  %>%  pull(Genes)

# Plot
draw_heatmap(
  filter_features_from_values(
    eset_count_filtered_log2,
    values = gene_list,
    column = "GeneID"
  ),
  column_annotation_names = c("Treatment", "Donor"),
  column_palette_selection = c("dark2", "nejm"),
  show_row_names = FALSE,
  name = "Exprs",
  scale_data = TRUE,
  column_title = "Heatmap of expression for genes DE in at least one contrast", 
  # cluster_rows = FALSE,
  # cluster_columns = FALSE,
  interactive = FALSE
)


## ---- fig.width= 8, fig.height=7, warnings = FALSE------------------------------------------------------------------------------------

# Creates a bogus eSet with contrast FC as expression
eset_forheatmap <- eset_object(
  expression_data = DEA_res %>% dplyr::select(contrast, data) %>% unnest() %>% pivot_wider(
    id_cols = Genes,
    names_from = contrast,
    values_from = log2FC
  ) %>% drop_na() %>%  data.frame(row.names = "Genes") %>% as.matrix() ,
  phenotype_data = data.frame(row.names = DEA_res$contrast, contrast = DEA_res$contrast),
  feature_data = fData(eset_count_filtered)
)

# Heatmap
draw_heatmap(
  eset_forheatmap %>% filter_features_from_values(
    values = DEA_res$DEgenes  %>% bind_rows()  %>%  pull(Genes),
    column = "GeneID"
  ),
  column_annotation_names = c("contrast"),
  column_palette_selection = c("dark2"),
  # cluster_column = F,
  row_annotation_names = NULL,
  show_row_names = FALSE, 
  name="log2FC", 
  column_title = "Heatmap of log2FC of genes DE in at least one contrast", 
  scale_data = FALSE # Do not scale log2FC
) %>% print()



## ---- fig.width=7, fig.height=5-------------------------------------------------------------------------------------------------------
draw_profile(
  eset_tpm,
  col_y = gene_interest,
  col_x = "Treatment",
  col_group = "Donor",
  col_color = "Donor",
  palette  = "accent",
  point_size = 1,
  hline_threshold  = log2(esti_tpm + 1),
  y_transformation  = ~ log2(. + 1)) + 
  labs(title = "CSTA gene")



## ---- fig.width=7, fig.height=5-------------------------------------------------------------------------------------------------------
# draw_profileplot(eset_tpm,
#                  x = "Treatment",
#                  y_filter ="A1BG",
#                  groupby = "Donor",
#                  y_unit = "tpm",
#                  exprs_trans = ~log2(. + 1),
#                  palette = "Dark2",
#                  background_threshold = NULL,
#                  interactive = TRUE)


## ---- fig.width = 4, fig.height = 4---------------------------------------------------------------------------------------------------

# Creates a named list with all contrasts
list_venn <- pmap(DEA_res, function(DEgenes, ...) {
  DEgenes %>% pull(Genes)
})

names(list_venn) <- DEA_res$contrast

# Ven diagram
NBC = length(list_venn) %>% viridis()
venn::venn(list_venn, zcolor = NBC)


## ---- fig.width = 7, fig.height = 5---------------------------------------------------------------------------------------------------
dataplot <- binary_membership(list_venn, elt_name = "Genes")

upset(dataplot, colnames(dataplot)[-1])



## ----load-modules---------------------------------------------------------------------------------------------------------------------
# mSigDB
m_df = msigdbr(species = "Homo sapiens", category = "H")
m_list = m_df %>% split(x = .$gene_symbol, f = .$gs_name)

enrichDB <- m_df %>%
  dplyr::select(gs_name, gene_symbol)

# Reactome
# enrichDB <- msigdbr::msigdbr(species = "Homo sapiens",
#                  category = "C2",
#                  subcategory = "REACTOME") %>% select(gs_name, gene_symbol)

# KEGG
# enrichDB <- msigdbr::msigdbr(species = "Homo sapiens",
#                  category = "C2",
#                  subcategory = "KEGG") %>% select(gs_name, gene_symbol)



## ----fisher-computation---------------------------------------------------------------------------------------------------------------
# Running Fisher test
DEA_res <- DEA_res %>% mutate(Fisher_results = purrr::map(
    DEgenes,
    ~ .x %>% pull(Genes) %>% fisher_enrichment(
    significative_genes_list = .  ,
    background_genes = rownames(eset_count),
    pathways_database = enrichDB
  )
  ))

# Calculation of the completeness of each pathway for graphical representation
completnesse_data <- completeness_calculation(
  analyzed_data = fData(eset_count),
  gene_colname_data = "GeneID",
  reference_database = enrichDB,
  pathway_colname_reference = "gs_name",
  gene_colname_reference =  "gene_symbol"
)


# adding formatted Fisher results object to DEA_res

DEA_res <- DEA_res %>% mutate(
  Fisher_results_4graph = purrr::map(
    Fisher_results,
    function(.x){
      if(is.null(.x)){
        NULL
      }else{
        .x %>% res_format_enrichment(., type = "fisher")
      }
    } 
  )
)

DEA_res <- DEA_res %>% mutate(
  Fisher_results_4graph = purrr::map(
    Fisher_results_4graph,
    function(.x){
      if(is.null(.x)){
        NULL
      }else{
        .x %>% left_join(., completnesse_data)
      }
    } 
  )
)



## ----compute-fGSEA-scores-------------------------------------------------------------------------------------------------------------
# Adding the fgseascores and results to the DEA_res object

# fGSEA scores
DEA_res <- DEA_res %>% mutate(fGSEA_scores = purrr::map(
  data,
  ~ .x %>% fgsea_ranked_score(
    deg_results = .,
    fc_label = "FC",
    pvalue_label = "pvalue",
    gene_label = "Genes"
  )
))


# Calculation of the completeness of each pathway for graphical representation
completnesse_data <- completeness_calculation(
  analyzed_data = fData(eset_count_filtered),
  gene_colname_data = "GeneID",
  reference_database = enrichDB,
  pathway_colname_reference = "gs_name",
  gene_colname_reference =  "gene_symbol"
)

# Running fGSEA

# IMPORTANT : fGSEA uses random walk and permutation procedures to estimate pvalues. It is necessary to set a seed to have reproducible results when runing several time this code

set.seed(1234)

# adding fgsea results object to DEA_res
DEA_res <- DEA_res %>% mutate(fGSEA_results = purrr::map(
  fGSEA_scores,
  ~ .x %>% fgsea_enrichment(
    statistics_score = .,
    pathways_database = enrichDB,
    pvalue_cutoff = 1
  ) 
))

# adding formatted fgsea results object to DEA_res

DEA_res <- 
  DEA_res %>% 
  mutate(fGSEA_results_4graph = purrr::map(
    fGSEA_results,
    ~ .x %>% res_format_enrichment(., type = "fgsea")
  )
)

DEA_res <- 
  DEA_res %>% 
  mutate(fGSEA_results_4graph = purrr::map(
    fGSEA_results_4graph,
    ~ .x %>% left_join(., completnesse_data)
  )
)



## ----contrast-table-fgsea-------------------------------------------------------------------------------------------------------------
# Full table of pathways for all contrasts

fulltable_pathways <-
  DEA_res %>% dplyr::select(contrast, fGSEA_results_4graph) %>% tidyr::unnest() %>%
  tidyr::pivot_wider(
    id_cols = c("ID", "pathway_size"),
    names_from = "contrast",
    values_from = c("NES", "FDR")
  )

# Display table
display_table(
  fulltable_pathways,
  interactive = TRUE,
  # color_shade_column = "FDR",
  bar_column = paste0("NES_", DEA_res$contrast),
  highlighted_color = "lightblue",
  nmax_display = 10
)


## ---- results='asis', fig.width= 8, fig.height=6--------------------------------------------------------------------------------------
purrr::pwalk(DEA_res,
             function(fGSEA_results_4graph, contrast, ...) {
               cat(paste0("  \n#### ", contrast , "  \n"))
               draw_dotplot_enrichment(
                 fGSEA_results_4graph,
                 xlabel = "NES",
                 pathway_colname = "ID",
                 pvalue_colname = "FDR",
                 display_dot_value = "Completeness",
                 pathway_size_colname = "pathway_size",
                 number_display = 10,
                 up_down_fgsea = TRUE,
                 type = "fgsea",
                 log_pvalue = TRUE
               ) %>% print()
               cat("\n")

             })


## ---- fig.width = 5, fig.height = 5---------------------------------------------------------------------------------------------------

draw_running_score (DEA_res$fGSEA_results[[1]], pathway_name = 1:3)



## ---- results='asis', fig.width= 8, fig.height=6--------------------------------------------------------------------------------------
purrr::pwalk(DEA_res ,
             function(Fisher_results_4graph, contrast, ...) {
               cat(paste0("  \n### ", contrast , "  \n"))
               if (!is.null(Fisher_results_4graph)) {
               draw_barplot_score (
                 Fisher_results_4graph,
                 xlabel = "GeneRatio",
                 pathway_colname = "ID",
                 pvalue_colname = "FDR",
                 display_bar_value = "Completeness",
                 number_display = 15,
                 up_down_fgsea = FALSE,
                 type = "fisher",
                 bar_width = 0.6, 
                 log_pvalue =  TRUE
               ) %>% print()
               }
               cat("\n")
               
             })

