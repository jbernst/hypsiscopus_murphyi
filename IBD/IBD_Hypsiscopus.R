getwd()
setwd("C:/Users/Justin/Documents/Publications/Hypsiscopus_murphyi/IBD/")


###Get matrix of geographic distances
library(geodist)
geo <- read.csv(file = "cyt-b_geo.csv", header = TRUE)
geo.dist <- geodist(geo, measure = "geodesic")

###Get matrix of genetic distances
library(ape)
library(phytools)
library(adegenet)

dna <- read.FASTA(file = "CYTB_NoOG.fasta", type = "DNA")
dna <- read.dna(file = "CYTB_NoOG.fasta", format = "fasta", as.matrix = TRUE)
class(dna)
View(dna)

dna.distances <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
#There will be one Na due to two small sequences, so turn it into a zero
dna.distances[is.na(dna.distances)] <- 0

###Mantel test
library(ggplot2)

DNA <- as.data.frame(dna.distances)
GEO <- as.data.frame(geo.dist)
View(DNA)

mantel.test(DNA, GEO)

lm <- lm(dna.distances~geo.dist)
summary(lm)
plot.lm
plot(dna.distances~geo.dist)
abline(lm(dna.distances~geo.dist))

#This shows there is clear IBD going on, but we should limit it to two species
ibd <- mantel.randtest(as.dist(DNA), as.dist(GEO))
plot(ibd)

#Another type of plot
plot(geo.dist, dna.distances)
abline(lm(geo.dist~dna.distances), col="red",lty=2)






#####
#####Let's repeat for just murphi and just plumbea
getwd()
setwd("C:/Users/Justin/Documents/Publications/Hypsiscopus_murphyi/IBD/")


###Get matrix of geographic distances
library(geodist)
geo.murphi <- read.csv(file = "murphi_geo.csv", header = TRUE)
geo.plumbea <- read.csv(file = "plumbea_geo.csv", header = TRUE)
geo.dist.m <- geodist(geo.murphi, measure = "geodesic")
geo.dist.p <- geodist(geo.plumbea, measure = "geodesic")

###Get matrix of genetic distances
library(ape)
library(phytools)
library(adegenet)

dna.m <- read.FASTA(file = "CYTB_NoOG_murphi.fasta", type = "DNA")
dna.p <- read.FASTA(file = "CYTB_NoOG_plumbea.fasta", type = "DNA")
class(dna)
View(dna)

dna.distances.m <- dist.dna(dna.m, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
dna.distances.p <- dist.dna(dna.p, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)

#There will be one Na due to two small sequences, so turn it into a zero
dna.distances.p[is.na(dna.distances.p)] <- 0
dna.distances.m[is.na(dna.distances.m)] <- 0

###Mantel test
library(ggplot2)

DNA.m <- as.data.frame(dna.distances.m)
GEO.m <- as.data.frame(geo.dist.m)
DNA.p <- as.data.frame(dna.distances.p)
GEO.p <- as.data.frame(geo.dist.p)
View(DNA)

#Run the mantel tests
ibd.m <- mantel.randtest(as.dist(DNA.m), as.dist(GEO.m))
ibd.p <- mantel.randtest(as.dist(DNA.p), as.dist(GEO.p))
plot(ibd.m)
plot(ibd.p)

#Another type of plot
plot(geo.dist.m, dna.distances.m)
plot(geo.dist.p, dna.distances.p)
abline(lm(as.numeric(geo.dist.m)~as.numeric(dna.distances.m)), col="red",lty=2)


###
###And now for plumbea and murphi together

###Get matrix of geographic distances
geo.plumbea.murphi <- read.csv(file = "plumbea_murphi_geo.csv", header = TRUE)
geo.dist.pm <- geodist(geo.plumbea.murphi, measure = "geodesic")


###Get matrix of genetic distances
dna.pm <- read.FASTA(file = "CYTB_NoOG_plumbea_murphi.fasta", type = "DNA")

dna.distances.pm <- dist.dna(dna.pm, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)

#There will be one Na due to two small sequences, so turn it into a zero



###Mantel test
DNA.pm <- as.data.frame(dna.distances.pm)
GEO.pm <- as.data.frame(geo.dist.pm)

View(DNA)

#Run the mantel tests
ibd.pm <- mantel.randtest(as.dist(DNA.pm), as.dist(GEO.pm))
plot(ibd.pm)

#Another type of plot
plot(geo.dist.pm, dna.distances.pm)
abline(lm(geo.dist.pm~dna.distances.pm), col="red",lty=2)

?mantel.randtest
citation("ade4")
citation("ape")
citation("geodist")
packageVersion("ape")
