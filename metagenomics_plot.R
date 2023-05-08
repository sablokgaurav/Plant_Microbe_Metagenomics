#' Title: ggplot function to plot all metagenomics experimental plots
#' # roxygen document parsing
#' @param file : file of the metagenomics classified
#' @param column1 : column of the variable
#' @param column2 : column of the variable
#' @param column3 : column of the facet
#' @param point
#' @param bar
#'
#' @return
#' @export
#'
#' @examples
#'
#' if (require("wdm", quietly = TRUE)) {
#' cor_test(iris, "Sepal.Length", "Sepal.Width",
#'       method = "blomqvist")⁄}
#'
file <- list.files(path = ".", pattern = "*.csv")
read_file <- lapply(file, read.csv)
for (i in read_file)
 print(i)
plot_metagenomics <-
 function(file, column1, column2,
          column3) {
  dataset <- read.csv("file", stringsAsFactors = FALSE)
  options <- c("point", "bar", "col")
  for (i in options) {
   if (i == "point") {
    ggplot(dataset,
           aes(dataset$column1,
               dataset$column2),
           fill = dataset$column1) +
     geom_point() +
     coord_flip() +
     facet_grid(. ~ dataset$column3) +
     theme_light()
    write.csv(dataset$column1, dataset$column2,
              dataset$column3, file = "metagenomics_columns_uses")
   } else if (i == "bar") {
    ggplot(dataset,
           aes(dataset$column1,
               dataset$column2),
           fill = dataset$column1) +
     geom_bar() +
     coord_flip() +
     facet_grid(. ~ column3) +
     theme_light()
    write.csv(read$column1, read$column2,
              read$column3, file = "metagenomics_columns_uses")
   } else if (i == "col") {
    ggplot(dataset, aes(column1, column2), fill = column1) +
     geom_col() +
     coord_flip() +
     facet_grid(. ~ column3) +
     theme_light()
    write.csv(read$column1, read$column2,
              read$column3, file = "metagenomics_columns_uses")
   }
  }
 }
#' Title this is a subset of the function which will prepare the
#'correlation matrix of the metagenomics samples according to the
#'studied variables in the experiment
#'
#' @param dataframe
#' @return
#' @export
#'
#' @examples
#' # if you have multiple files then use this
file <- list.files(path = ".", pattern = "*.csv")
read_file <- lapply(file, read.csv)
for (i in read_file)
 print(i)
# here read_file[1] is the first csv file, you can select
# any number
colnames(read_file[1])
str_replace_all(colnames(read_file[1]), "test_dataset_1.csv.", "")
colnames(read_file[1]) <- str_replace_all(colnames(read_file[1]),
                                          "test_dataset_1.csv.", "")
# you can select the columns of the species
rownames(read_file[1]) <- c(as.character(unlist
                                         (colnames(read_file[1][1]))))
head(read_file[1])
cor(read_file[1])
cor.test.p <- function(x){
 FUN <- function(x, y) cor.test(x, y)[["p.value"]]
 z <- outer(
  colnames(x), 
  colnames(x), 
  Vectorize(function(i,j) FUN(x[,i], x[,j]))
 )
 dimnames(z) <- list(colnames(x), colnames(x))
 z
}
options = c("scatter", "heatmap")
for (i in options) {
 heatmaply_cor(
  read_file[1],
  node_type = i,
  point_size_mat = -log10(p),
  point_size_name = "-log10(p-value)",
  label_names = c("x", "y", "Correlation")
 )
}

#' Title This function will make the Manhattan plot for the metagenomics
#' experiment from a different view point.
#
#' @param dataset
#'
#' @return
#' @export
#'
#' @examples
file <- list.files(path = ".", pattern = "*.csv")
read_file <- lapply(file, read.csv)
for (i in read_file)
 print(i)
select_dataset <- read_file[1][c("bacterial_class", "variable1",
                                 "variable2", "variable3")]
# you can tag as many variable as you want
options = c("c", "b", "cunique", "bunique")
for (i in options) {
 if (i == "c") {
  CMplot(
   select_dataset,
   plot.type = "i",
   threshold =
    c(0.01 / 0.05) / length(unlist(
     select_dataset$bacterial_class,
     use.names = FALSE
    )),
   threshold.col = c('red', 'orange', 'blue'),
   multracks = TRUE,
   pch = 10,
   box = TRUE,
   file.output = TRUE,
   file.name = "C_CMplot",
   file = "jpg",
   dpi = 300,
   height = 4.3,
   width = 2.2
  )
 } else if (i == "b") {
  CMplot(
   select_dataset,
   plot.type = "i",
   threshold =
    c(0.01 / 0.05) / length(unlist(
     select_dataset$bacterial_class,
     use.names = FALSE
    )),
   threshold.col = c('red', 'orange', 'blue'),
   multracks = TRUE,
   pch = 10,
   box = TRUE,
   file.output = TRUE,
   file.name = "C_CMplot",
   file = "jpg",
   dpi = 300,
   height = 4.3,
   width = 2.2
  )
 } else if (i == "cunique") {
  CMplot(
   select_dataset,
   plot.type = "i",
   threshold =
    c(0.01 / 0.05) / length(unique(
     unlist(select_dataset$bacterial_class,
            use.names = FALSE)
    )),
   threshold.col = c('red', 'orange', 'blue'),
   multracks = TRUE,
   pch = 10,
   box = TRUE,
   file.output = TRUE,
   file.name = "C_CMplot",
   file = "jpg",
   dpi = 300,
   height = 4.3,
   width = 2.2
  )
 } else if (i == "bunique") {
  CMplot(
   select_dataset,
   plot.type = "i",
   threshold =
    c(0.01 / 0.05) / length(unique(
     unlist(select_dataset$bacterial_class,
            use.names = FALSE)
    )),
   threshold.col = c('red', 'orange', 'blue'),
   multracks = TRUE,
   pch = 10,
   box = TRUE,
   file.output = TRUE,
   file.name = "C_CMplot",
   file = "jpg",
   dpi = 300,
   height = 4.3,
   width = 2.2
  )
 }
}