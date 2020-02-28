# Load packages to the library
library(vegan)
library(factoextra)
library(lavaan)
library(semPlot)
library(semTools)
library(dplyr)
library(scales)
library(gridExtra)

#Set the working directory
setwd("V:/LEWA Synthesis/Analyses")

# Load the data
new_data <- read.csv(file="LEWADataPath3.csv", header=TRUE)
names(new_data)

# Subset the variables
data_ten <- subset(new_data, Year==2010)
data_ten <- data_ten[,-cbind(19,27)]
data_nine <- subset(new_data, Year==2009)
data_nine <- data_nine[,-cbind(19,27)]
#data_nine <- na.omit(data_nine)
#data_ten <- na.omit(data_ten)

ten <- data_ten[,9:20]
ten$pH <- data_ten[,5]
#ten$Nutrients <- data_ten[,22]
ten <- ten[,-cbind(3,5,7,8)]

nine <- data_nine[,9:20]
nine$pH <- data_nine[,5]
#nine$Nutrients <- data_nine[,23]
nine <- nine[,-cbind(3,7,5,8)]

# Perform PCA
nine_pca <- prcomp(nine, scale = TRUE)
ten_pca <- prcomp(ten, scale = TRUE)

#Add PCA to the datafiles for analyses
ynine <-cbind(data_nine[,2], data_nine[,7:8], data_nine[,23:27],data_nine[,21], nine_pca$x[,1], nine_pca$x[,2]*-1)
colnames(ynine) <- c("Glyphosate","Min.Depth", "Max.Depth", "Richness", "Abundance", "Surface.Area", "pr.predators", "chiro.emerge", "Cover","PC1.Chemicals","PC2.Chemicals")

yten <-cbind(data_ten[,2], data_ten[,7:8], data_ten[,23:27], data_ten[,21], ten_pca$x[,1:2])
colnames(yten) <- c("Glyphosate","Min.Depth", "Max.Depth", "Richness", "Abundance", "Surface.Area", "pr.predators", "chiro.emerge", "Cover","PC1.Chemicals","PC2.Chemicals")

#Making PCA figures and calculating eigenvalues to report
#2009 figure

eig.val <- as.data.frame(get_eigenvalue(nine_pca)) #get the Eigenvalues
eig.val

res.ind <- get_pca_ind(nine_pca) #Get the coords for individual points
vals <- as.data.frame(res.ind$coord[,1:2]) #we only want the first two axis

res.var <- get_pca_var(nine_pca) #Get the coords for end of arrows
res.var$contrib #contributions of each variable to the first two axis
vars <- as.data.frame(res.var$coord[,1:2]) #Coordinates of the variables on the first two axis

vars[,1] <- rescale(vars[,1], to = c(min(vals[,1]),max(vals[,1]))) #rescale to the max and min of the coords of individual points
vars[,2] <- rescale(vars[,2], to = c(min(vals[,2]),max(vals[,2]))) #rescale to the max and min of the coords of individual points

veclab <- c("Cl", "Ca", "Fe", "Mg", "Na", "NH4N", "NO3N", "P", "pH")
veclabloc <- vars

veclabloc[1,1] <- veclabloc[1,1] + .4
veclabloc[1,2] <- veclabloc[1,2] - .05
veclabloc[2,1] <- veclabloc[2,1] + .6
veclabloc[2,2] <- veclabloc[2,2]
veclabloc[3,1] <- veclabloc[3,1] - .5
veclabloc[3,2] <- veclabloc[3,2] + .1
veclabloc[4,1] <- veclabloc[4,1] + .6
veclabloc[4,2] <- veclabloc[4,2]
veclabloc[5,1] <- veclabloc[5,1] + .6
veclabloc[5,2] <- veclabloc[5,2]
veclabloc[6,1] <- veclabloc[6,1]
veclabloc[6,2] <- veclabloc[6,2] - .4
veclabloc[7,1] <- veclabloc[7,1]
veclabloc[7,2] <- veclabloc[7,2] - .2
veclabloc[8,1] <- veclabloc[8,1] - .1
veclabloc[8,2] <- veclabloc[8,2] - .2
veclabloc[9,1] <- veclabloc[9,1] + .6
veclabloc[9,2] <- veclabloc[9,2] + .3

#Make the basic 2009 plot
biplot2009 <- ggplot() +
  geom_hline(yintercept = 0, size =1) + #adds horizontal line
  geom_vline(xintercept = 0, size =1) + #adds verical line
  geom_point(data=vals, aes(x=Dim.1, y=Dim.2), colour="grey") + #plots the points
  geom_segment(data=vars, aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1) + #plots the arrows
  geom_text(data=veclabloc, label = veclab, aes(x=Dim.1, y=Dim.2), alpha=1, size=3)+ #lables the arrows
  scale_x_continuous(sprintf("PC1 (%s%%)", round(eig.val[1,2],digits=2)), limits = c(-7,7)) +
  scale_y_continuous(sprintf("PC2 (%s%%)", round(eig.val[2,2],digits=2)), limits = c(-5.5,5.5)) +
  annotate("text", label = "A", x = -6.4, y = 5.5, size = 4) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#Making the 2010 plot

eig.val.2010 <- as.data.frame(get_eigenvalue(ten_pca)) #get the Eigenvalues
eig.val.2010

res.ind.2010 <- get_pca_ind(ten_pca) #Get the coords for individual points
vals.2010 <- as.data.frame(res.ind.2010$coord[,1:2]) #we only want the first two axis

res.var.2010 <- get_pca_var(ten_pca) #Get the coords for end of arrows
res.var.2010$contrib #contributions of each variable to the first two axis
vars.2010 <- as.data.frame(res.var.2010$coord[,1:2]) #we only want the first two axis

vars.2010[,1] <- rescale(vars.2010[,1], to = c(min(vals.2010[,1]),max(vals.2010[,1]))) #rescale to the max and min of the coords of individual points
vars.2010[,2] <- rescale(vars.2010[,2], to = c(min(vals.2010[,2]),max(vals.2010[,2]))) #rescale to the max and min of the coords of individual points

veclab.2010 <- c("Cl", "Ca", "Fe", "Mg", "Na", "NH4N", "NO3N", "P", "pH")
veclabloc.2010 <- vars.2010

veclabloc.2010[1,1] <- veclabloc.2010[1,1] + .55
veclabloc.2010[1,2] <- veclabloc.2010[1,2]
veclabloc.2010[2,1] <- veclabloc.2010[2,1] + .6
veclabloc.2010[2,2] <- veclabloc.2010[2,2] + .05
veclabloc.2010[3,1] <- veclabloc.2010[3,1] - .1
veclabloc.2010[3,2] <- veclabloc.2010[3,2] - .2
veclabloc.2010[4,1] <- veclabloc.2010[4,1] + .2
veclabloc.2010[4,2] <- veclabloc.2010[4,2] - .2
veclabloc.2010[5,1] <- veclabloc.2010[5,1] + .6
veclabloc.2010[5,2] <- veclabloc.2010[5,2]
veclabloc.2010[6,1] <- veclabloc.2010[6,1] - 1.6
veclabloc.2010[6,2] <- veclabloc.2010[6,2] + .3
veclabloc.2010[7,1] <- veclabloc.2010[7,1] + 1.7
veclabloc.2010[7,2] <- veclabloc.2010[7,2] + .1
veclabloc.2010[8,1] <- veclabloc.2010[8,1]
veclabloc.2010[8,2] <- veclabloc.2010[8,2] + .3
veclabloc.2010[9,1] <- veclabloc.2010[9,1] + .55
veclabloc.2010[9,2] <- veclabloc.2010[9,2] - .2

#Make the basic plot
biplot2010 <- ggplot() +
  geom_hline(yintercept = 0, size =1) + #adds horizontal line
  geom_vline(xintercept = 0, size =1) + #adds verical line
  geom_point(data=vals.2010, aes(x=Dim.1, y=Dim.2), colour="grey") + #plots the points
  geom_segment(data=vars.2010, aes(x=0, y=0, xend=Dim.1, yend=Dim.2), arrow=arrow(length=unit(0.2,"cm")), alpha=1) + #plots the arrows
  geom_text(data=veclabloc.2010, label = veclab.2010, aes(x=Dim.1, y=Dim.2), alpha=1, size=3)+ #lables the arrows
  scale_x_continuous(sprintf("PC1 (%s%%)", round(eig.val.2010[1,2],digits=2)), limits = c(-8.5,8.5)) +
  scale_y_continuous(sprintf("PC2 (%s%%)", round(eig.val.2010[2,2],digits=2)), limits = c(-4,4)) +
  annotate("text", label = "B", x = -8, y = 4, size = 4) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

#Create the combined PCA figure
jpeg("PCABiplots_Revised.jpeg", width = 14, height = 7, units="cm", res=300)
grid.arrange(biplot2009, biplot2010, ncol=2)
dev.off()



####END OF PCA----------------------------------------------------####


###SEM ANALYSES####

#centre and scale the data
#This really isn't necessary if we use the standardized Beta's provided
ynine_mc <- as.data.frame(scale(ynine, center = TRUE, scale = TRUE))
yten_mc <- as.data.frame(scale(yten, center = TRUE, scale = TRUE))

#Old version of the analyses used in the first draft of the manuscript
All2009_SEM <- '
  Cover ~ Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals
  Richness ~ Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals + Cover
  Abundance ~ Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals + Cover + Richness
'
fit_All2009 <- sem(All2009_SEM, data = ynine_mc) 
summary(fit_All2009, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)

#Including the direct effect of glyphosate on Abundance and the indirect effect through changing cover
All2009_SEM_2 <- '
  Cover ~ a*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals
  Richness ~ Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals + Cover
  Abundance ~ c*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals + b*Cover + Richness
  
  ind_abun := a*b
  dir_abun := c
  tot_abun := c + (a*b)
'

fit_All2009_2 <- sem(All2009_SEM_2, data = ynine) 
summary(fit_All2009_2, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)

#Including the direct effect of glyphosate on Abundance and the indirect effect through changing cover
#Including the direct effect of glyphosate on richness and the indirect effect through changing cover
All2009_SEM_3 <- '
  Cover ~ a*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals
  Richness ~ e*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals + d*Cover
  Abundance ~ c*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals + b*Cover + Richness

  ind_rich := a*d
  dir_rich := d
  tot_rich := e + (a*d)

  ind_abun := a*b
  dir_abun := c
  tot_abun := c + (a*b)
'

fit_All2009_3 <- sem(All2009_SEM_3, data = ynine_mc) 
summary(fit_All2009_3, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)

#Specifying a model that only contains the Paths that we defined and includes both indirect pathways
SEM_4 <- '
  Cover ~ a*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals
  Richness ~ e*Glyphosate + Min.Depth + Max.Depth + Surface.Area + d*Cover
  Abundance ~ c*Glyphosate + Min.Depth + Max.Depth + Surface.Area + b*Cover + Richness

  ind_rich := a*d
  dir_rich := e
  tot_rich := e + (a*d)

  ind_abun := a*b
  dir_abun := c
  tot_abun := c + (a*b)
'

fit_All2009_4 <- sem(SEM_4, data = ynine) 
summary(fit_All2009_4, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)

fit_All2010_4 <- sem(SEM_4, data = yten) 
summary(fit_All2010_4, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)


fit_All2009_5 <- sem(SEM_5, data = ynine_mc) #need to mean center and scale because variances are very different
summary(fit_All2009_5, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)


#Model that fits benthic predators and includes the indirect effect of glyphosate > cover > predators > amphibian
ynine_mc <- as.data.frame(scale(ynine, center = TRUE, scale = TRUE)) # need to mean centre and scale
SEM_6 <- '
  Cover ~ a*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals + PC2.Chemicals
  pr.predators ~ g*Glyphosate + Min.Depth + Max.Depth + Surface.Area + f*Cover
  Richness ~ e*Glyphosate + Min.Depth + Max.Depth + Surface.Area + d*Cover + pr.predators
  Abundance ~ c*Glyphosate + Min.Depth + Max.Depth + Surface.Area + b*Cover + Richness + i*pr.predators


  ind_rich := a*d
  dir_rich := e
  tot_rich := e + (a*d)

  ind_abun := a*b
  dir_abun := c
  tot_abun := c + (a*b)

  ind_pabun := g*(a*f)
  dir_pabun := i
  tot_pabun := i + (g*(a*f))

  ind_pred := a*f
  dir_pred := g
  tot_pred := g + (a*f)
'

fit_All2009_6 <- sem(SEM_6, data = ynine_mc) #need to mean center and scale because variances are very different
summary(fit_All2009_6, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)
#AIC = 122.215 chi square p = 0

#Predator model does not fit the data very well, p from chi square < 0.05.
#I reduce the number of factors by removing the smallest non-significant effects for each response until the model fits
#I check AIC to make sure the model is explaining more as I go along
#All variables are left for amphibian abundance because they are all boarderline significant
SEM_7 <- '
  Cover ~ a*Glyphosate + Min.Depth + Max.Depth + Surface.Area + PC1.Chemicals
  pr.predators ~ g*Glyphosate + Min.Depth + Surface.Area + f*Cover
  Richness ~ e*Glyphosate + Min.Depth + Max.Depth + Surface.Area + d*Cover
  Abundance ~ c*Glyphosate + Min.Depth + Max.Depth + Surface.Area + b*Cover + Richness + i*pr.predators

  ind_rich := a*d
  dir_rich := e
  tot_rich := e + (a*d)

  ind_abun := a*b
  dir_abun := c
  tot_abun := c + (a*b)

  ind_pabun := g*(a*f)
  dir_pabun := i
  tot_pabun := i + (g*(a*f))

  ind_pred := a*f
  dir_pred := g
  tot_pred := g + (a*f)
'

fit_All2009_7 <- sem(SEM_7, data = ynine_mc)
summary(fit_All2009_7, standardized=TRUE, fit.measures = TRUE, rsq = TRUE)
#chi square p = 0.308
#AIC = 116.356
