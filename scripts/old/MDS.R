library(ggplot2)
library(MASS)
library(vegan)
library(dplyr)
library(magrittr)
library(data.table)
library(readr)

#read data into variable data
data <- read_csv("delims/pCERC4_plasmid_coverage.csv")

data <- as.data.frame(data)

nems <- data$X1

data <- data[,2:ncol(data)]

rownames(data) <- nems

#data <- t(data[,2:ncol(data)])

data <- as.matrix(data)

data <- unique(data)


virdist <- vegdist(wisconsin(data), method = "jaccard")
# Get the baseline solution: start with cmdscale
mds.null.vir <- isoMDS(virdist, tol=1e-7, maxit = 200)
## See if you can get any better.
repeat{
        mds.1.vir <- isoMDS(virdist, k=2, initMDS(virdist, k = 2), maxit=200, trace=FALSE, tol=1e-7)
        print("another iteration")
        print(mds.1.vir$stress)
        print(mds.null.vir$stress)
        if(mds.1.vir$stress < mds.null.vir$stress) break
}








# Scale solutions ("fix translation, rotation and scale")
mds.null.vir <- postMDS(mds.null.vir, virdist)
mds.1.vir <- postMDS(mds.1.vir, virdist)
# Compare solutions
plot(procrustes(mds.1.vir, mds.null.vir))
#convert to dataframe
mds.null.vir.df <- as.data.frame(mds.null.vir,
                                 attr(x = mds.null.vir,
                                      which = "dimnames"))
mds.1.vir.df <- as.data.frame(mds.1.vir, 
                              attr(x = mds.1.vir, 
                                   which = "dimnames"))

rownames(mds.1.vir.df) -> mds.1.vir.df$name
mds.1.vir.df <- left_join(mds.1.vir.df, mds.metadata)

mds.1.vir.df$ST <- as.factor(mds.1.vir.df$ST)
