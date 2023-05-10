#!/usr/bin/env Rscript --vanilla
#' Title this will generate all the seq_logs, seq_plot for the
#' multi fasta file of the selected metagenomics or the genomics
#' sample and will also check for the model test and phylogeny
#' inference
#' comment this line if you want to install the packages
# if (require(argparse)) {
#  install.packages("argparse")
# }
# if (require(methods)) {
#  install.packages("methods")
# }
# if (require(ggmsa)) {
#  install.packages("ggmsa")
# }
# if (require(phangorm)) {
#  install.packages("phangorm")
# }
# if (require(ape)) {
#  install.packages("ape")
# }
suppressPackageStartupMessages(library(argparser, pos = "package:base"))
suppressPackageStartupMessages(library(methods, pos = "package:base"))
suppressPackageStartupMessages(library(ggmsa, pos = "package:base"))
suppressPackageStartupMessages(library(phangorn, pos = "package:base"))
suppressPackageStartupMessages(library(ape, pos = "package:base"))

p <- arg_parser("genomic_alignment_map")
p <-
 add_argument(p, "--alignment", help = "Name of the alignment file")
p <- add_argument(p, "--NT", help = "NT alignment in file")
p <- add_argument(p, "--RNA", help = "Protein alignment in file")
p <- add_argument(p, "--start", help = "start of alignment")
p <- add_argument(p, "--end", help = "end of alignment")
p <- add_argument(p, "--tree", help = "provide a calibrated tree")
argv <- parse_args(p)

alignment_file <- argv$alignment
NT <- argv$NT
RNA <- argv$RNA
start <- argv$start
end <- argv$end
tree <- argv$tree

if (is.na(alignment_file) & is.na(NT)) {
 stop("No alignment file provided and no type selected")
}

if (file.access(alignment_file) == -1)
{
 stop("No alignment file provided and ")
}

if (isFALSE(c("NT", "RNA", "start", "end"))) {
 stop("No parameters for the NT, RNA, start, end
      are provided")
}

if (alignment_file & NT) {
 alignment_plot <- ggmsa(
  alignment_file,
  start,
  end,
  color = "Shapely_NT",
  font = "DroidSansMono",
  char_width = 0.5,
  seq_name = TRUE
 ) +
  geom_seqlogo(color = "Shapely_NT") + geom_msaBar() + geom_GC()
}
save.image(file = "alignment_plot", version = 1)
model_evaluation <- modelTest(
 alignment_file,
 tree,
 model = "c("JC", "F81", "K80", "HKY","SYM", "GTR"))
 model_all <- modelTest(alignment_file,
                  model = "all", multicore = TRUE,mc.cores = 2)

 if (!file.exists("alignment_model_optimization")){
  file.create("model_optimization.txt")}

 write(model_evaluation, file = "model_optimization.txt")
 if (!file.exists("alignment_model_optimization_all")){
  file.create("model_optimization_all")}

 write(model_all, file = "model_optimization_all")

if (alignment_file & RNA) {
alignment_plot <- ggmsa(
  alignment_file,
  start,
  end,
  color = "Shapely_NT",
  font = "DroidSansMono",
  char_width = 0.5,
  seq_name = TRUE
) + geom_seqlogo(color = "Shapely_NT" ) + geom_msaBar() + geom_GC()
}
 save.image(file = "alignment_plot", version = 1)
 model_evaluation <- modelTest(alignment_file,
            tree, model = "c("JC", "F81", "K80", "HKY", "SYM", "GTR")
)
model_all <- modelTest(
 alignment_file,
 model = "all",
 multicore = TRUE,
 mc.cores = 2
)

if (!file.exists("alignment_model_optimization")) {
 file.create("model_optimization.txt")
}

write(model_evaluation, file = "model_optimization.txt")
if (!file.exists("alignment_model_optimization_all")) {
 file.create("model_optimization_all")
}
write(model_all, file = "model_optimization_all")