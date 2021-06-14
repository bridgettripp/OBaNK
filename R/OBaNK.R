#' @title package
#'
#' @desciption Learns interactions from omics data
#'
#' @param Matrix of omics data (feature x sample)
#' @param Optional: 2) Seed value (default = 378) 3) Output folder (default = results)
#'
#' @return 1) Graphical rendering, 2) Learned interactions, 3) Summary
#'
#' @example OBaNK(input_df = "data_matrix", seed_val = 1132, output_folder = "results_folder")
#'
#' @export


OBaNK <- function(input_df, seed_val, output_folder){

  if(missing(seed_val)) {
    myseed <- 378
  } else {
    myseed <- seed_val
  }
  if(missing(output_folder)) {
    results_folder <- "results"
    if (file.exists(results_folder)){
      results_folder <- results_folder
    }
    else {
      dir.create(results_folder)
      results_folder <- results_folder
    }

  } else {
      if (file.exists(output_folder)){
        results_folder <- output_folder
      }
      else {
        dir.create(output_folder)
        results_folder <- output_folder
      }
    }

  suppressMessages(library(bnlearn))
  suppressMessages(library(parallel))
  suppressPackageStartupMessages(library(Rgraphviz))
  suppressMessages(library(plyr))

  load("EvidenceData.Rda")
  source("confidenceValue.R")
  source("foldChangeCalc.R")

  options(warn = -1)

  start.time <- as.character(date())

  d1 <- read.table(input_df, header = TRUE,stringsAsFactors = FALSE, sep = ",")
  rownames(d1) <- d1[,2]
  name.id <- d1[-1,1:2]
  idcol <- as.character(d1[,2])
  groups <- as.factor(d1[1,-c(1,2)])

  if(length(unique(groups)) > 2){
    stop("Too many experimental groups. Only two are allowed.")}
  if(length(levels(groups)) < 2){
    stop("Not enough experimental groups. Two groups required")}

  fc <- linFold(d1[,-(1:2)]) #foldchange x feature
  colnames(fc) <- idcol[-1]

  dsc3 <- bnlearn::discretize(fc, breaks = 3, method = 'hartemink')

  dt <- dsc3
  dt2 <- "Hartemink"
  al1 <- "mmhc"
  al2 <- mmhc

  output.file <- file(paste0(results_folder,"/", "summary_generated_BN.txt"), open="wt")
  sink(output.file, append = TRUE)
  print(paste("STARTED:", start.time))

  cl = makeCluster(4)
  z <- round(.85*(dim(dt)[1]))
  set.seed(myseed)
  strength <- boot.strength(dt, R = 1000, m = z, algorithm = al1, cpdag = TRUE, cluster = cl) ##R = number of bootstraps m = size of each bootstrap replicate
  stopCluster(cl)

  strg <- (confidenceValue(strength$strength)*.75)
  print(paste0("strength = ", strg))

  name.id <- as.data.frame(name.id)

  learned.arcs <- as.data.frame(cbind(strength$from, strength$to, strength$strength))
  colnames(learned.arcs) <- c("from","to", "streng")

  strength.length <- length(strength$strength)
  remove.perc = .05
  sub.perc = .75
  sorted.strength <- as.matrix(sort(strength$strength, decreasing = TRUE))
  subset.strength <- sorted.strength[! sorted.strength %in% 0] ## remove zero values

  set.seed(myseed)
  norm.bank <- rnorm(length(subset.strength),
                     mean=subset.strength[length(subset.strength)*remove.perc],
                     sd = ((subset.strength[length(subset.strength)*remove.perc])/4))
  sub.bank <- rnorm(length(subset.strength),
                    mean=subset.strength[length(subset.strength)*sub.perc],
                    sd = ((subset.strength[length(subset.strength)*sub.perc]))/4)

  adds <- data.frame()
  subs <- data.frame()

  for (i in 1:dim(learned.arcs)[1]){
    if (length(which(learned.arcs$from[i]==evidence[,1] & learned.arcs$to[i]==evidence[,2])) == 1){
      mch <- which(learned.arcs$from[i]==strength$from & learned.arcs$to[i]==strength$to)
      strength$strength[mch] <- round((as.numeric(strength$strength[mch]) + sample(norm.bank, 1, replace=TRUE)),2)
      if (strength$direction[mch] < 0.5){
        strength$direction[mch] <- 0.5
        }
      adds <- rbind(adds, strength[mch,])
      if(strength$strength[mch]>=1){
        strength$strength[mch] <- 1}
    }
    else {
      if (as.numeric(as.character(learned.arcs$streng[i])) >= strg){
        strength$strength[i] <- round((as.numeric(strength$strength[i]) - sample(sub.bank, 1, replace=TRUE)),2)
        subs <- rbind(subs,strength[i,] )
      }
    }
  }

  write.csv(strength, paste0(results_folder,"/","Strength.csv"))

  avg.boot2 <- averaged.network(strength, threshold = (strg))
  print(avg.boot2)

  pdf(paste0(results_folder,"/","graphical_rendering.pdf"))
  we <- as.data.frame(avg.boot2$arcs)
  netwe <- ugraph(ftM2graphNEL(as.matrix(we), edgemode = "directed"))
  g1 <- layoutGraph(netwe)
  graph.par(list(nodes=list(fill="gray", textCol="blue", fontsize=10), arrowhead = "normal"))
  graphviz.plot(skeleton(as.bn(g1)), shape = "ellipse", layout = "neato",
                main = paste("Equivalency Class: IDs"), render = TRUE)

  for (i in seq(length(strength$strength))){
    strength$from[i] <- name.id$Symbol[BiocGenerics::which(strength$from[i] == name.id$ID)]
    strength$to[i] <- name.id$Symbol[BiocGenerics::which(strength$to[i] == name.id$ID)]}

  avg.boot2 <- averaged.network(strength, threshold = (strg))

  we <- as.data.frame(avg.boot2$arcs)
  netwe <- ugraph(ftM2graphNEL(as.matrix(we), edgemode = "directed"))
  g1 <- layoutGraph(netwe)
  graph.par(list(nodes=list(fill="gray", textCol="blue", fontsize=10), arrowhead = "normal"))
  graphviz.plot(skeleton(as.bn(g1)), shape = "ellipse", layout = "neato",
                main = paste("Equivalency Class"), render = TRUE)

  print(paste("COMPLETED:", date()))

  dev.off()
  sink()
  close(output.file)

  options(warn = 0)
}





