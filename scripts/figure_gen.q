
###################

library(ggplot2)
library(monocle)

exprs.experiment.1 = readMM("experiment.1.exprs")
exprs.experiment.1 = as(exprs.experiment.1, "dgCMatrix")

exprs.experiment.2 = readMM("experiment.2.exprs")
exprs.experiment.2 = as(exprs.experiment.2, "dgCMatrix")

fData.experiment.1 = read.table(
    "experiment.1.fData",
    header = T, row.names = 1, stringsAsFactors = F)

fData.experiment.2 = read.table(
    "experiment.2.fData",
    header = T, row.names = 1, stringsAsFactors = F)

pData.experiment.1 = read.table(
    "experiment.1.pData",
    header = T, row.names = 1, stringsAsFactors = F, sep = "\t")

pData.experiment.2 = read.table(
    "experiment.2.pData",
    header = T, row.names = 1, stringsAsFactors = F, sep = "\t")

cds.experiment.1 = newCellDataSet(
    exprs.experiment.1,
    featureData = new("AnnotatedDataFrame", fData.experiment.1),
    phenoData = new("AnnotatedDataFrame", pData.experiment.1),
    expressionFamily = negbinomial.size())

dim(cds.experiment.1)

cds.experiment.2 = newCellDataSet(
    exprs.experiment.2,
    featureData = new("AnnotatedDataFrame", fData.experiment.2),
    phenoData = new("AnnotatedDataFrame", pData.experiment.2),
    expressionFamily = negbinomial.size())

dim(cds.experiment.2)

rm(list = c(
    "exprs.experiment.1", "exprs.experiment.2",
    "fData.experiment.1", "fData.experiment.2",
    "pData.experiment.1", "pData.experiment.2"))

# Figure 3a
png("worm_t-sne.png", units="in", height=2, width=2, res=600)
ggplot(pData(cds.experiment.1), aes(x = tsne_1, y = tsne_2, color = cluster.name)) +
    geom_point(size = 0.002) +
    xlab("") + ylab("") +
    guides(color = guide_legend(
        title = "Cell type\n(aggregated)",
        override.aes = list(size = 4))) +
    theme_void() +
    monocle:::monocle_theme_opts() +
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.margin = margin(0, -10, 0, 10),
          legend.key.width=unit(0.15, "in"),
          legend.key.height=unit(0.15, "in"),
          legend.position = "none")
dev.off()


###################
png("HSMM_trajectory.png", width=2, height=2, units="in", res=600)
plot_spanning_tree(BJ_MYO_selected, color_by="as.factor(Time)", cell_size=0.5, show_branch_points=FALSE) + 
  theme(legend.position = "none") + 
   stat_density2d(color="gray30", size=I(0.15)) +
   theme_void() +
    theme(legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.margin = margin(0, -10, 0, 10),
          legend.key.width=unit(0.15, "in"),
          legend.key.height=unit(0.15, "in"),
          legend.position = "none") +
  scale_color_manual(name="Time from Serum Switch", 
                       values = c("#EF5B5B", "#FFBA49",  "#8ABF69", "#0FA3B1"), 
                       labels=c("0 Hours", "24 Hours", "48 Hours", "72 Hours")) + 
  guides(colour = guide_legend(title.position = "top", title.hjust =0.5, keywidth=0.1, keyheight=.1)) 
dev.off()
