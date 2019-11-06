---
title: "MSMB-Chapter5-Clustering"
author: "Aleeza Gerstein"
date: '2019-11-06'
output:
  rmarkdown::pdf_document:
  fig_caption: yes
includes:
  in_header: preamble-latex.tex
header-includes:
  - \usepackage{color}
  - \usepackage{mathtools}
---

```{r global_options, include=FALSE}
knitr::opts_chunk$set(fig.pos = 'H')
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("mixtools"))
suppressPackageStartupMessages(library("flowCore"))
#BiocManager::install("ggcyto")
#BiocManager::install("labeling")
#BiocManager::install("dbscan")
#BiocManager::install("pheatmap")
suppressPackageStartupMessages(library("ggcyto"))
suppressPackageStartupMessages(library("labeling"))
suppressPackageStartupMessages(library("dbscan"))
suppressPackageStartupMessages(library("pheatmap"))
```

```{r}
v1 <- seq(0, 1, length.out=100)
plot(log(v1), asinh(v1), type="l")
plot(v1, asinh(v1), type= "l")
v3 <- seq(30, 3000, length = 100)
plot(log(v3), asinh(v3), type="l")
````

#5.5 Clustering examples: Flow cytometry

```{r}
library("flowCore")
library("flowViz")
fcsB <- read.FCS("data/Bendall_2011.fcs")
markersB <- readr::read_csv("data/Bendall_2011_markers.csv")
mt <- match(markersB$isotope, colnames(fcsB))
colnames(fcsB)[mt] <- markersB$marker

asinhtrsf <- arcsinhTransform(a = 0.1, b = 1)
fcsBT <- transform(fcsB, transformList(colnames(fcsB)[-c(1, 2, 4)], asinhtrsf))

kf <- kmeansFilter("CD3all" = c("Pop1", "Pop2"), filterID = "myKmFilter")
fres <- flowCore::filter(fcsBT, kf)
summary(fres)
fcsBT1 <- flowCore::split(fcsBT, fres, population = "Pop1")
fcsBT2 <- flowCore::split(fcsBT, fres, population = "Pop2")
```

```{r}
library("ggcyto")
ggcd4cd8 <- ggcyto(fcsB, aes(x = CD4, y = CD8))
ggcd4 <- ggcyto(fcsB, aes(x = CD4))
ggcd8 <- ggcyto(fcsB, aes(x = CD8))
p1 <- ggcd4 +
  geom_histogram(bins = 60)
p1b <- ggcd8 +
  geom_histogram(bins = 60)
asinht <- arcsinhTransform(a = 0, b =1)
trans1 <- transformList(colnames(fcsB)[-c(1, 2, 4)], asinht)
fcsBT <- transform(fcsB, trans1)
p1t <- ggcyto(fcsBT, aes(x=CD4)) +
  geom_histogram(bins = 90)
p1t

p2t <- ggcyto(fcsBT, aes(x=CD4, y = CD8)) +
  geom_density2d(colour="black")
p2t

p3t <- ggcyto(fcsBT, aes(x = CD45RA, y = CD20)) +
  geom_density2d(colour = "black")
p3t
```

## 5.5.3 Density-based clustering
```{r}
mc5 <- exprs(fcsBT)[,c(15, 16, 19, 40, 33)] #why this order?
res5 <- dbscan::dbscan(mc5, eps = 0.65, minPts = 30)
mc5df <- data.frame(mc5, cluster = as.factor(res5$cluster))
table(mc5df$cluster)
```

```{r}
mc6 <- exprs(fcsBT)[,c(15, 16, 19, 40, 25, 33)] #why this order?
res6 <- dbscan::dbscan(mc6, eps = 0.65, minPts = 20)
mc6df <- data.frame(mc6, cluster = as.factor(res6$cluster))
table(mc6df$cluster)

res6_b <- dbscan::dbscan(mc6, eps = 0.45, minPts = 20)
mc6df_b <- data.frame(mc6, cluster = as.factor(res6_b$cluster))
table(mc6df_b$cluster)

res6_c <- dbscan::dbscan(mc6, eps = 0.75, minPts = 20)
mc6df_c <- data.frame(mc6, cluster = as.factor(res6_c$cluster))
table(mc6df_c$cluster)
```

\textcolor{red}{Question 5.8}
```{r}
load("data/Morder.RData")
#without dendogram or reordering, Euclidean and Manhattan distances

library(gridExtra)
library(grid)
#I want to plot both on the same panel
#two ways of doing it
#https://www.biostars.org/p/128229/
#https://stackoverflow.com/questions/39590849/using-a-pheatmap-in-arrangegrob

plot_list <- list()
ph1 <- pheatmap(Morder)
ph2 <- pheatmap(Morder, clustering_distance_rows = "manhattan")
plot_list[[1]] <- ph1[[4]]
plot_list[[2]] <- ph2[[4]]

g <- grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
#g<-do.call(grid.arrange,plot_list)

#Question 5.9: which orderings do not match
#https://www.biostars.org/p/170614/
res1 <- Morder[c(ph1$tree_row[["order"]]),ph1$tree_col[["order"]]]
res2 <- Morder[c(ph2$tree_row[["order"]]),ph2$tree_col[["order"]]]
row.names(res1) == row.names(res2)
```

#here
\textcolor{blue}{Question 5.12}
```{r}
library(tidyverse)
simdat <- lapply(c(0, 8), function(mx){
  lapply(c(0, 8), function(my){
    tibble(x = rnorm(100, mean = mx, sd = 2),
           y = rnorm(100, mean = my, sd = 2),
           class = paste(mx, my, sep=":"))
  }) %>% bind_rows
}) %>% bind_rows
simdat
simdatxy <- simdat[, c("x", "y")]

rsim1 <- range(c(simdat$x[1:100], simdat$y[1:100]))
rsim2 <- range(c(simdat$x[101:200], simdat$y[101:200]))
rsim3 <- range(c(simdat$x[201:300], simdat$y[201:300]))
rsim4 <- range(c(simdat$x[301:400], simdat$y[301:400]))

simdat_un <- lapply(c(0, 8), function(mx){
  lapply(c(0, 8), function(my){
    tibble(x = runif(100, rsim[1], rsim[2]),
           y = runif(100, rsim[1], rsim[2]),
           class = paste(mx, my, sep=":"))
  }) %>% bind_rows
}) %>% bind_rows
simdat_un
```

```{r}
library(cluster)

pamfun <- function(x, k)list(cluster=pam(x, k, cluster.only = TRUE))

gss <- clusGap(simdatxy, FUN = pamfun, K.max = 8, B = 50, verbose = FALSE)

plot_gap <- function(x){
  gstab <- data.frame(x$Tab, k = seq_len(nrow(x$Tab)))
  ggplot(gstab, aes(k, gap)) +
    geom_line()+
    geom_errorbar(aes(ymax = gap + SE.sim,
                      ymin = gap - SE.sim), width = 0.1) +
    geom_point(size = 3, col = "red")
}

plot_gap(gss)

library("Hiiragi2013")
data("x")

selFeats = order(rowVars(Biobase::exprs(x)), decreasing = TRUE)[1:50]
embmat = t(Biobase::exprs(x)[selFeats, ])
embgap = clusGap(embmat, FUN = pamfun, K.max = 24, verbose = FALSE)
k1 = maxSE(embgap$Tab[, "gap"], embgap$Tab[, "SE.sim"])
k2 = maxSE(embgap$Tab[, "gap"], embgap$Tab[, "SE.sim"],
           method = "Tibs2001SEmax")
c(k1, k2)
```

\textcolor{red}{Question 5.15 Use all features in x}

```{r}
````
