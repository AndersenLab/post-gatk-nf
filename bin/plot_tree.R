# plot tree
library(ggplot2)
library(ape)
library(ggtree)

args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(args[1])

# how long should the tree be?
id <- length(tree$tip.label)

ggtree::ggtree(tree)+
  geom_tiplab(size = 2)
ggsave(glue::glue("{args[1]}.pdf"), height = id/8.5, width = 20, limitsize = FALSE)
