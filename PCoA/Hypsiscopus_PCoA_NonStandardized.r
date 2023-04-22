# Set working directory
getwd()
setwd("C:/Users/Justin/Documents/Publications/Hypsiscopus_murphi/Quant_Stats/")

# Load packages
library(vegan)
library(labdsv)
library(ggplot2)
library(dplyr)

Hypsiscopus <- read.csv("Hypsiscopus_MASTER_Morph-Data_PCOA_Adults.csv", header = TRUE, na.strings = c("","NA"))
View(Hypsiscopus)

testpcoa <- Hypsiscopus[c(11:13,15,17:43)]
View(testpcoa)

#calculate distance matrix with Gower transform in Vegan
Gowerdist <- vegdist(testpcoa,method="gower")

#spat out error about missing values, fixed with:
Gowerdist <- vegdist(testpcoa,method="gower", na.rm= TRUE)

#error, must be numeric
sapply(Hypsiscopus, class)  #lets you know which columns are factors (a no-no)
sapply(testpcoa, class)  #lets you know which columns are factors (a no-no)
View(testpcoa)
Gowerdist

#can remove all rows with NAs as follows:
testpcoa.full <- na.omit(testpcoa)
View(testpcoa.full)
str(testpcoa)
str(testpcoa.full)
View(testpcoa.full)

#Now make a new Gowerdist object with the new testpcoa.full object
Gowerdist.full <- vegdist(testpcoa.full,method="gower", na.rm= TRUE)
View(Gowerdist.full)

# Documentation is a little sparse, From doc re: na.rm "logical. Should missing values (including NaN) be omitted from the calculations?"
# Need to check and see if this is solved by ignoring columns. Can do by making dataset omitting columns with missing data to see if identical results.
# Note: na.rm does not remove columns, only cells

#Run PCoA in labdsv saving first four dimensions
pcotest <- pco(Gowerdist.full,k=2)
View(pcotest)
#Save output as .csv
#Commented out since we already have this and don't want to overwrite it
#write.csv(pcotest$points,'pcotestPOINTSAdults.csv')

#Visualize plot

#Don't forget to append a species column to the outputted .csv (also, can turn the species into numbers to make a gradient color scheme)
pcoaTESTplot <- read.csv("pcotestPOINTS.csv")       

#pcoaTESTplot.4sp <- read.csv("pcotestPOINTS_4sp.csv") #3 species; can also substitute with the file 'pcotestPOINTS_3sp.csv' (this latter file is not made yet)
ggplot(pcoaTESTplot, aes(x=V1, y=V2, color=species)) + geom_point()

#("#7ECEFA", "#93E07F", "#B971F0", "#E0A618"
ggplot(pcoaTESTplot, aes(x=V1, y=V2, color=species)) + geom_point() + geom_text(aes(label=species),hjust=0, vjust=0) +
  scale_color_manual(values=c("#7ECEFA", "#93E07F")) 

#Code block for experimenting with drawing hulls
#calculating a convex polygon hull for each data group
hull_fr<-pcoaTESTplot%>%
  group_by(species)%>%
  slice(chull(V1,V2))
#now create base scatterplot
fr.pco.plot<-ggplot(pcoaTESTplot,aes(x=V1,y=V2,colour=species))+
  geom_point()  +
  scale_color_manual(values=c("#7ECEFA", "#93E07F"))
#put in a fill group on the plot and overlay the hulls
fr.pco.plot+aes(fill=factor(species))+geom_polygon(data=hull_fr,alpha=0.3)
#THIS BLOCK WORKS FOR CONVEX HULLS
fr.pco.plot+aes(fill=factor(species))+geom_polygon(data=hull_fr,alpha=0.3) +
  scale_fill_manual(values=c("#7ECEFA", "#93E07F"))

#If you want numbers of the specimens, run this
pcoaTESTplot <- read.csv("pcotestPOINTSnumbers.csv")       

#pcoaTESTplot.4sp <- read.csv("pcotestPOINTS_4sp.csv") #3 species; can also substitute with the file 'pcotestPOINTS_3sp.csv' (this latter file is not made yet)
ggplot(pcoaTESTplot, aes(x=V1, y=V2, color=x)) + geom_point()

#("#7ECEFA", "#93E07F", "#B971F0", "#E0A618"
ggplot(pcoaTESTplot, aes(x=V1, y=V2, color=x)) + geom_point() + geom_text(aes(label=x),hjust=0, vjust=0) 

#Code block for experimenting with drawing hulls
#calculating a convex polygon hull for each data group
hull_fr<-pcoaTESTplot%>%
  group_by(species)%>%
  slice(chull(V1,V2))
#now create base scatterplot
fr.pco.plot<-ggplot(pcoaTESTplot,aes(x=V1,y=V2,colour=species))+
  geom_point()  +
  scale_color_manual(values=c("#7ECEFA", "#93E07F"))
#put in a fill group on the plot and overlay the hulls
fr.pco.plot+aes(fill=factor(species))+geom_polygon(data=hull_fr,alpha=0.3)
#convex hulls
fr.pco.plot+aes(fill=factor(species))+geom_polygon(data=hull_fr,alpha=0.3) +
  scale_fill_manual(values=c("#7ECEFA", "#93E07F"))

#Making a 3D, rotating video of your PCA
library(scatterplot3d)
library(rgl)
library(ggplot2)
library(ggfortify)
#library(ggbiplot)
library(readr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(heplots)
library(candisc)
library(rpart)
library(rpart.plot)
library(corrplot)
library(magrittr)
library(pca3d)
library(scatterplot3d)
library(plot3D)
library(cluster)

# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1_3<-c("#E69F00", "#56B4E9", "#009E73")
cbp1_4<-c("#E69F00", "#56B4E9", "#009E73", "#999999")
cbp1_5<-c("#E69F00", "#56B4E9", "#009E73", "#999999", "#CC79A7")
cbp1_4b<-c("#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1_6<-c("#E69F00", "#56B4E9", "#009E73", "#999999", "#D55E00", "#CC79A7")

scatter3d(V1~V2+V3|as.factor(species),data=pcoaTESTplot,
          surface=FALSE, sphere.size=0.8,ellipsoid=TRUE, level=0.5, 
          surface.col=c("#E69F00", "#56B4E9", "#009E73"))

#No ellipsoids
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot$V2, z = pcoaTESTplot$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.3sp$Species), grid = FALSE, surface = FALSE)


#With ellipsoids : this isn't working, but may be due to how many points per grouping we have
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot$V2, z = pcoaTESTplot$V3,
          point.col = "blue", groups = as.factpr(pcoaTESTplot.3sp$Species), ellipsoid = TRUE, grid = TRUE, surface = TRUE)


#No ellipsoids 3 species version
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot.3sp$V2, z = pcoaTESTplot.3sp$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.3sp$Species), grid = FALSE, surface = FALSE)

pcoaTESTplot.3sp <- read.csv("pcotestPOINTS_3sp.csv")
#With ellipsoids 3 species version
scatter3d(x = pcoaTESTplot.3sp$V1, y = pcoaTESTplot.3sp$V2, z = pcoaTESTplot.3sp$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.3sp$Species), ellipsoid = TRUE, grid = TRUE, surface = FALSE)

#With ellipsoids 4 species version: This version works because before cf. geminaulis and cf. calligaster had too few points to plot as an ellipse

pcoaTESTplot.4sp <- read.csv("pcotestPOINTS_4sp.csv")
scatter3d(x = pcoaTESTplot.4sp$V1, y = pcoaTESTplot.4sp$V2, z = pcoaTESTplot.4sp$V3,
          point.col = "blue", groups = as.factor(pcoaTESTplot.4sp$Species), ellipsoid = TRUE, grid = TRUE, surface = FALSE)

#Rotation Gif
library(tidyverse)
library(dplyr)
library(rgl)
library(ggpubr)
library(heplots)
library(candisc)
library(rpart)
library(rpart.plot)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(magrittr)
library(pca3d)
library(scatterplot3d)
library(plot3D)
library(cluster)
library(sf)
library(sp)
library(GISTools)

scatter3d(V2~V1+V3|as.factor(Species),data=pcoaTESTplot.4sp, fogtype ="none", fov=20,  box=TRUE, xlab="V1",ylab="V2",zlab="V3",
          bg.col= "black",
          axis.col=c("white", "white", "white"), axis.scales=FALSE, axis.ticks=TRUE, text.col="pink",
          surface=FALSE, fill=TRUE, grid=TRUE, sphere.size=0.8, ellipsoid=TRUE, level=0.5,
          surface.col=c("#7ECEFA", "#93E07F", "#B971F0", "#E0A618"), font="10")


grid3d(c("x", "y", "z"), at = NULL, col = "gray", lwd = 1, lty = 1, n = 5)

#You must leave the temporary window from line 165 (scatter3d) open and the images WILL be the size of the window; maximize the window for best results and THEN run the export code below
movie3d(spin3d(axis = c(0, 1, 0)), duration = 12,
        dir = "C:/Users/Justin/Documents/Publications/Hemibungarus_NewRecord_2019/Quant_Stats/Hemibungarus_PCoA/", convert=TRUE)


##duration = 12 gave me a full rotation for my data


