#### Bernoulli GAM WITH SPATIAL CORRELATION
#############################################################################





######################################################



######################################################
#Load packages and support files
source("~/CAMBRIDGE/6. Data_analysis/INLA_kenya/Lattice_data/HighstatLibV13.R")
library(lattice)
library(ggplot2)
library(ggmap)
library(mgcv)
library(lme4)
library(rgdal)
library(sp)
library(gstat)
library(raster)
library(plyr)
library(reshape)
library(fields)
library(maps)
library(maptools)
library(mapdata)
library(plyr)
data("worldHiresMapEnv")

# And load extra packages:
library(INLA)
library(brinla)
library(MASS)
library("RColorBrewer")

library(ggregplot)

## edit function
Efxplot.Val <- function (ModelList, Sig = TRUE, StarLoc = NULL, Alpha1 = 1, 
                         Alpha2 = 1, PointOutline = F, ModelNames = NULL, VarNames = NULL, 
                         VarOrder = NULL, Intercept = TRUE, Size = 2, tips = 0.2) 
{
  require(dplyr)
  require(ggplot2)
  require(INLA)
  require(MCMCglmm)
  Graphlist <- list()
  if (!class(ModelList) == "list") {
    ModelList <- list(ModelList)
  }
  for (i in 1:length(ModelList)) {
    model <- ModelList[[i]]
    if (class(model) == "inla") {
      Graph <- as.data.frame(summary(model)$fixed)
      Graph <- Graph[1:4,] ## edit Valentina
      colnames(Graph)[which(colnames(Graph) %in% c("0.025quant", 
                                                   "0.975quant"))] <- c("Lower", "Upper")
      colnames(Graph)[which(colnames(Graph) %in% c("0.05quant", 
                                                   "0.95quant"))] <- c("Lower", "Upper")
      colnames(Graph)[which(colnames(Graph) %in% c("mean"))] <- c("Estimate")
    }
    if (class(model) == "MCMCglmm") {
      Graph <- as.data.frame(summary(model)$solutions)
      colnames(Graph)[1:3] <- c("Estimate", "Lower", 
                                "Upper")
    }
    Graph$Model <- i
    Graph$Factor <- rownames(Graph)
    Graphlist[[i]] <- Graph
  }
  Graph <- bind_rows(Graphlist)
  Graph$Sig <- with(Graph, ifelse(Lower * Upper > 0, "*", 
                                  ""))
  Graph$Model <- as.factor(Graph$Model)
  if (!is.null(ModelNames)) {
    levels(Graph$Model) <- ModelNames
  }
  position <- ifelse(length(unique(Graph$Model)) == 1, "none", 
                     "right")
  if (is.null(VarOrder)) 
    VarOrder <- rev(unique(Graph$Factor))
  if (is.null(VarNames)) 
    VarNames <- VarOrder
  Graph$Factor <- factor(Graph$Factor, levels = VarOrder)
  levels(Graph$Factor) <- VarNames
  Graph %<>% as.data.frame %>% filter(!is.na(Factor))
  if (!Intercept) {
    VarNames <- VarNames[!str_detect(VarNames, "ntercept")]
    Graph <- Graph %>% filter(Factor %in% VarNames)
  }
  Graph$starloc <- NA
  min <- min(Graph$Lower, na.rm = T)
  max <- max(Graph$Upper, na.rm = T)
  if (Sig == TRUE) {
    Graph$starloc <- max + (max - min)/10
  }
  if (!is.null(StarLoc)) {
    Graph$starloc <- StarLoc
  }
  Graph$Alpha <- with(Graph, ifelse(Lower * Upper > 0, Alpha1, 
                                    Alpha2))
  Graph <- Graph %>% mutate(SigAlpha = factor(as.numeric(Lower * 
                                                           Upper > 0), levels = c(1, 0)))
  if (PointOutline) {
    PointOutlineAlpha <- Alpha1
  }
  else {
    PointOutlineAlpha <- 0
  }
  ggplot(Graph, aes(x = as.factor(Factor), y = Estimate, group = Model, 
                    colour = Model, alpha = SigAlpha)) + 
    theme(text = element_text(size = 18)) +
    geom_point(position = position_dodge(w = 0.5), 
               size = Size,
               col = "navyblue") + geom_errorbar(position = position_dodge(w = 0.5), 
                                                 aes(ymin = Lower, ymax = Upper), 
                                                 size = 1.2, 
                                                 width = tips, col = "navyblue") + 
    geom_hline(aes(yintercept = 0), lty = 2) + labs(x = NULL) + 
    coord_flip() + theme(legend.position = position) + geom_text(aes(label = Sig, 
                                                                     y = starloc), position = position_dodge(w = 0.5), show.legend = F, col="blue") + 
    scale_alpha_manual(values = c(Alpha1, Alpha2)) + guides(alpha = "none") + 
    geom_point(colour = "black", aes(group = Model), 
               position = position_dodge(w = 0.5), size = 4, alpha = PointOutlineAlpha) + 
    geom_errorbar(aes(ymin = Lower, ymax = Upper, group = Model), 
                  width = 0.1, position = position_dodge(w = 0.5), 
                  colour = "black", alpha = PointOutlineAlpha) + 
    geom_point(position = position_dodge(w = 0.5), size = 3, 
               alpha = PointOutlineAlpha)
}



####################################################################
#Set the working directory and load the data
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")

#Import the data 
mydata <- read.csv(file = "PresenceAbsence_100kmBuffer.csv", 
                   header = TRUE,
                   na.strings = "NA",
                   stringsAsFactors = TRUE,
                   dec = ".")
colnames(mydata)
mydata$longitude <- mydata$xcoord
mydata$latitude <- mydata$ycoord

#Inspect the file
names(mydata)
str(mydata)  # Make sure you have num and not factors for 
summary(mydata)
dim(mydata)
###################################################################





##################################################################
# Underlying task. 
# Model total anthrax outbreaks as a function of covariates: 



##################################################################
# DATA EXPLORATION

## 1.1: Outliers 
#**********************
# in the response variable (outbreaks) and the covariates
MyVar <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10",
           "bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio1",
           "elevation","distance_to_water","soil_calcium","soil_pH","soilwater",
           "soil_carbon", "evi" )
Mydotplot(mydata[, MyVar])



## 1.2: Collinearity
#**********************
MyVar <- data.frame(
  Annual.Mean.Temp = mydata$bio1, 
  Mean.Diurnal.Range = mydata$bio2, 
  Isothermality = mydata$bio3, 
  Temp.Seasonality = mydata$bio4, 
  Max.Temp.Warmest.Month = mydata$bio5, 
  Minimum.Temp.Coldest.Month = mydata$bio6,
  Temp.Annual.Range = mydata$bio7, 
  Mean.Temp.Wettest.Quarter = mydata$bio8, 
  Mean.Temp.Driest.Quarter = mydata$bio9, 
  Mean.Temp.Warmest.Quarter = mydata$bio10, 
  Mean.Temp.Coldest.Quarter = mydata$bio11, 
  Annual.Precipitation = mydata$bio12, 
  Precipitation.Wettest.Month = mydata$bio13, 
  Precipitation.Driest.Month = mydata$bio14, 
  Precipitation.Seasonality = mydata$bio15, 
  Precipitation.Wettest.Quarter = mydata$bio16, 
  Precipitation.Driest.Quarter = mydata$bio17, 
  Precipitation.Warmest.Quarter = mydata$bio18,
  Precipitation.Coldest.Quarter = mydata$bio19, 
  Elevation = mydata$elevation, 
  Distance.to.Water = mydata$distance_to_water, 
  Soil.Calcium = mydata$soil_calcium, 
  Soil.Water = mydata$ soilwater,
  Soil.Carbon = mydata$soil_carbon, 
  EVI = mydata$evi,
  Soil.pH = mydata$soil_pH)

###
MyVar1 <- data.frame(
  BIO1 = mydata$bio1, 
  BIO2 = mydata$bio2, 
  BIO3 = mydata$bio3, 
  BIO4 = mydata$bio4, 
  BIO5 = mydata$bio5, 
  BIO6 = mydata$bio6,
  BIO7 = mydata$bio7, 
  BIO8 = mydata$bio8, 
  BIO9 = mydata$bio9, 
  BIO10 = mydata$bio10, 
  BIO11 = mydata$bio11, 
  BIO12 = mydata$bio12, 
  BIO13 = mydata$bio13, 
  BIO14 = mydata$bio14, 
  BIO15 = mydata$bio15, 
  BIO16 = mydata$bio16, 
  BIO17 = mydata$bio17, 
  BIO18 = mydata$bio18,
  BIO19 = mydata$bio19)

MyVar2 <- data.frame(
  BIO4 = mydata$bio4, 
  BIO12 = mydata$bio12, 
  Elevation = mydata$elevation, 
  Distance.to.Water = mydata$distance_to_water, 
  Soil.Calcium = mydata$soil_calcium, 
  Soil.Water = mydata$ soilwater,
  Soil.Carbon = mydata$soil_carbon, 
  EVI = mydata$evi,
  Soil.pH = mydata$soil_pH)

#*## Corrplot
#??corrplot

library(corrplot)
## corrplot 0.92 loaded
M = cor(MyVar1)
corrplot(M, method = 'number') # colorful number

M = cor(MyVar2)
corrplot(M, method = 'number') # colorful number



# Variables were removed gradually until the Variance Inflation Factor (VIF) values 
# were less than 4.
# Used VIF values, Pearson correlations, and scatterplots

#Variance inflation factors
MyVar <- c("bio4","bio12","elevation","distance_to_water","soil_calcium",
           "soilwater","soil_carbon", "evi","soil_pH" )
ke.vif <- corvif(mydata[,MyVar])

# remove bio12, carbon, evi, and soil pH,
MyVar <- c("bio4","distance_to_water","elevation","soil_calcium","soilwater")
ke.vif <- corvif(mydata[,MyVar])

## Use Pearson correlations and scatterplots:
library(GGally)
ggpairs(mydata[,MyVar])


## 1.3: Relationships
#**********************
# Relationships between the response variable and covariates
MyVar <- c("bio4","distance_to_water","elevation","soil_calcium","soilwater")
MyMultipanel.ggp2(Z = mydata, 
                  varx = MyVar, 
                  vary = "outbreaks", 
                  ylab = "Total outbreaks",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)


## Distibution of outbreak counts across locations?
par(mfrow=c(1,1))
plot(table(mydata$outbreaks), type = "h", 
     xlab = "Number of outbreaks", 
     ylab = "Number of locations", 
     main = "Frequency of locations with given outbreak counts")
round(100 * table(mydata$outbreaks)/nrow(mydata), digits = 2)
nrow(mydata[(mydata$outbreaks==0),])/nrow(mydata)

###################################################

##INLA Analysis:

###################################################
# What are the distances between the sites in km ?
mydata$Xkm   <- mydata$xcoord  #/ 1000 
mydata$Ykm   <- mydata$ycoord  #/ 1000
Loc <- mydata[, c("Xkm", "Ykm")] #Spatial locations
Distances <- dist(Loc)
hist(Distances)

# We have lots of sites positioned very close
# to each other. And there is where you would
# expect spatial dependency!


#######################################################
# House keeping
# Because this is a GAM with an exponential relationship, it is better
# to standardize each continuous covariate. 
mydata$bio4c    <- MyStd(mydata$bio4)
mydata$bio12c    <- MyStd(mydata$bio12)
mydata$elevationc    <- MyStd(mydata$elevation)
mydata$distance_to_waterc    <- MyStd(mydata$distance_to_water)
mydata$soil_calciumc    <- MyStd(mydata$soil_calcium)
mydata$soilwaterc    <- MyStd(mydata$soilwater)

# Bernoulli part:
mydata$outbreaks01 <- ifelse(mydata$outbreaks==0, 0, 1) # new response value
######################################################################






######################################################################
### Bernoulli GAM WITH SPATIAL CORRELATION

# Basis functions for the distance to water  data:
Legodistance_to_water <- smoothCon(mgcv::s(distance_to_water, bs = "cr", k = 4, fx = TRUE), ### edit
                      data= mydata,
                      absorb.cons = TRUE)[[1]]
# Extract the actual basis functions.
Xdistance_to_water <- Legodistance_to_water$X 
ncol_b <- ncol(Xdistance_to_water) # all basis functions have similar number of columns
nrow_b <- nrow(Xdistance_to_water)
# Give the basis functions unique names.
colnames(Xdistance_to_water)   <- paste("distance_to_water", 1:ncol_b, sep = "")
lcs.distance_to_water <- inla.make.lincombs(distance_to_water1 = Xdistance_to_water[,"distance_to_water1"],
                                            distance_to_water2 = Xdistance_to_water[,"distance_to_water2"],
                                            distance_to_water3 = Xdistance_to_water[,"distance_to_water3"])
names(lcs.distance_to_water)   <- paste(names(lcs.distance_to_water), "distance_to_water", sep = "")

###################
# Basis functions for the elevation  data:
Legoelevation <- smoothCon(mgcv::s(elevation, bs = "cr", k = 4, fx = TRUE), 
                           data= mydata,
                           absorb.cons = TRUE)[[1]]
Xelevation <- Legoelevation$X 
# Give the basis functions unique names.
colnames(Xelevation)   <- paste("elevation", 1:ncol_b, sep = "")
lcs.elevation <- inla.make.lincombs(elevation1 = Xelevation[,"elevation1"],
                                    elevation2 = Xelevation[,"elevation2"],
                                    elevation3 = Xelevation[,"elevation3"])
names(lcs.elevation)   <- paste(names(lcs.elevation), "elevation", sep = "")

###################
# Basis functions for the soil_calcium  data:
Legosoil_calcium <- smoothCon(mgcv::s(soil_calcium, bs = "cr", k = 4, fx = TRUE), 
                           data= mydata,
                           absorb.cons = TRUE)[[1]]
Xsoil_calcium <- Legosoil_calcium$X 
# Give the basis functions unique names.
colnames(Xsoil_calcium)   <- paste("soil_calcium", 1:ncol_b, sep = "")
lcs.soil_calcium <- inla.make.lincombs(soil_calcium1 = Xsoil_calcium[,"soil_calcium1"],
                                       soil_calcium2 = Xsoil_calcium[,"soil_calcium2"],
                                       soil_calcium3 = Xsoil_calcium[,"soil_calcium3"])
names(lcs.soil_calcium)   <- paste(names(lcs.soil_calcium), "soil_calcium", sep = "")


###################
# Basis functions for the soilwater  data:
Legosoilwater <- smoothCon(mgcv::s(soilwater, bs = "cr", k = 4, fx = TRUE), 
                      data= mydata,
                      absorb.cons = TRUE)[[1]]
Xsoilwater <- Legosoilwater$X 
# Give the basis functions unique names.
colnames(Xsoilwater)   <- paste("soilwater", 1:ncol_b, sep = "")
lcs.soilwater <- inla.make.lincombs(soilwater1 = Xsoilwater[,"soilwater1"],
                                    soilwater2 = Xsoilwater[,"soilwater2"],
                                    soilwater3 = Xsoilwater[,"soilwater3"])
names(lcs.soilwater)   <- paste(names(lcs.soilwater), "soilwater", sep = "")

###################
# Basis functions for the bio4 data:
Legobio4 <- smoothCon(mgcv::s(bio4, bs = "cr", k = 4, fx = TRUE), 
                      data= mydata,
                      absorb.cons = TRUE)[[1]]
Xbio4 <- Legobio4$X 
# Give the basis functions unique names.
colnames(Xbio4)   <- paste("bio4", 1:ncol_b, sep = "")
lcs.bio4 <- inla.make.lincombs(bio41 = Xbio4[,"bio41"],
                               bio42 = Xbio4[,"bio42"],
                               bio43 = Xbio4[,"bio43"])
names(lcs.bio4)   <- paste(names(lcs.bio4), "bio4", sep = "")


###################
## Organize covariates and basis functions in a dataframe:
N=nrow(mydata)
Xm <- data.frame(
  Intercept           = rep(1, N),
  bio4c              = mydata[,"bio4c"],
  distance_to_waterc  = mydata[,"distance_to_waterc"],
  elevationc          = mydata[,"elevationc"],
  soil_calciumc       = mydata[,"soil_calciumc"],
  soilwaterc          = mydata[,"soilwaterc"],
  distance_to_water1 = Xdistance_to_water[,"distance_to_water1"],
  distance_to_water2 = Xdistance_to_water[,"distance_to_water2"],
  distance_to_water3 = Xdistance_to_water[,"distance_to_water3"],
  elevation1 = Xelevation[,"elevation1"],
  elevation2 = Xelevation[,"elevation2"],
  elevation3 = Xelevation[,"elevation3"],
  soil_calcium1 = Xsoil_calcium[,"soil_calcium1"],
  soil_calcium2 = Xsoil_calcium[,"soil_calcium2"],
  soil_calcium3 = Xsoil_calcium[,"soil_calcium3"],
  soilwater1 = Xsoilwater[,"soilwater1"],
  soilwater2 = Xsoilwater[,"soilwater2"],
  soilwater3 = Xsoilwater[,"soilwater3"]
  )

####################################################



####################################################

## mydata
#1. Make a mesh.
#   Get a sense for the distribution of distances between sampling locations. 
Loc <- cbind(mydata$Xkm, mydata$Ykm)
plot(mydata$Xkm, mydata$Ykm)
D   <- dist(Loc)
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Distance between sites (km)",
     ylab = "Frequency")
# Create the mesh
RangeGuess <- 20 #km
MaxEdge  <- RangeGuess / 5
ConvHull <- inla.nonconvex.hull(Loc)		# convex=4
mesh    <- inla.mesh.2d(loc = Loc,
                        boundary = ConvHull,
                        max.edge = c(1, 5) * MaxEdge, 
                        #max.edge = c(2, 6) * MaxEdge,  #Use during the trial
                        cutoff  = MaxEdge / 5)
# plot
# Read data
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/shapefiles")
ug <- readOGR("Uganda_boundary_UTM_km.shp")

par(mfrow = c(1, 1), mar=c(2, 2, 2, 2))

plot(mesh)
plot(ug, col = "lightgrey", add=T)
plot(mesh, add=T)
points(Loc, col = "red", pch = 16, cex = 1.0)
mesh$n # [1] 5716


# 2. Define the weighting factors a_ik (also called 
#    the projector matrix).
A <- inla.spde.make.A(mesh, loc = Loc)

# 3. Define the SPDE.
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(20, 0.5), #, 100, 0.05 
                            prior.sigma = c(0.01,  0.05)) #0.01, 0.05

# 4. Define the spatial field.
w.index <- inla.spde.make.index(
  name    = 'w', 
  n.spde  = spde$n.spde)

# 5. Make a stack. 
# Here is the stack (bernoulli)
N <- nrow(mydata)
Stack.c <- inla.stack(
  tag = "Fit",
  data = list(outbreaks = mydata$outbreaks01),  
  A = list(1, 1, A),                 
  effects = list(  
    Intercept = rep(1, N),
    Xm        = Xm[,-1],
    w         = w.index))

# Run
Mb1.c <- inla(formula = outbreaks ~ -1 + Intercept, 
            family = "binomial", 
            control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
            data = inla.stack.data(Stack.c),
            control.predictor = list(A = inla.stack.A(Stack.c)))

Mb2.c <- inla(formula = outbreaks ~ -1 + Intercept + 
                f(w, model = spde) , 
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))

Mb3.c <- inla(formula = outbreaks ~ -1 + Intercept + bio4c + distance_to_waterc + 
                 elevationc + soil_calciumc + soilwaterc + 
                f(w, model = spde) , 
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))

Mb4.c <- inla(formula = outbreaks ~ -1 + Intercept + bio4c + distance_to_waterc + 
                 elevationc + soil_calciumc + soilwaterc, 
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))

Mb5.c <- inla(formula = outbreaks ~ -1 + Intercept + distance_to_waterc + 
                elevationc + soil_calciumc + soilwaterc, 
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))


Mb6.c <- inla(formula = outbreaks ~ -1 + Intercept + 
                distance_to_water1 + distance_to_water2 + distance_to_water3 +
                elevation1 + elevation2 + elevation3 +
                soil_calcium1 + soil_calcium2 + soil_calcium3 +
                soilwater1 + soilwater2 + soilwater3, 
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))

Mb7.c <- inla(formula = outbreaks ~ -1 + Intercept + 
                distance_to_water1 + distance_to_water2 + distance_to_water3 +
                elevation1 + elevation2 + elevation3 +
                soil_calcium1 + soil_calcium2 + soil_calcium3 +
                soilwater1 + soilwater2 + soilwater3 + 
                f(w, model = spde), 
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))

Mb8.c <- inla(formula = outbreaks ~ -1 + Intercept + distance_to_waterc +
                soil_calciumc + soilwaterc +
                elevation1 + elevation2 + elevation3, 
              lincomb = c(lcs.elevation),
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))

Mb9.c <- inla(formula = outbreaks ~ -1 + Intercept + distance_to_waterc +
                soil_calciumc + soilwaterc +
                elevation1 + elevation2 + elevation3 +
                f(w, model = spde), 
              lincomb = c(lcs.elevation),
              family = "binomial", 
              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),    
              data = inla.stack.data(Stack.c),
              control.predictor = list(A = inla.stack.A(Stack.c)))




## calculate AUC
library(PresenceAbsence)

DATA <- data.frame(
  ID = c(1:nrow(mydata)),
  observed = mydata$outbreaks01,
  predicted1 = Mb1.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted2 = Mb2.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted3 = Mb3.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted4 = Mb4.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted5 = Mb5.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted6 = Mb6.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted7 = Mb7.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted8 = Mb8.c$summary.fitted.values[1:nrow(mydata),"mean"],
  predicted9 = Mb9.c$summary.fitted.values[1:nrow(mydata),"mean"]
  )

auc1 <- auc(DATA, st.dev = TRUE, which.model = 1, na.rm = FALSE)
auc2 <- auc(DATA, st.dev = TRUE, which.model = 2, na.rm = FALSE)
auc3 <- auc(DATA, st.dev = TRUE, which.model = 3, na.rm = FALSE)
auc4 <- auc(DATA, st.dev = TRUE, which.model = 4, na.rm = FALSE)
auc5 <- auc(DATA, st.dev = TRUE, which.model = 5, na.rm = FALSE)
auc6 <- auc(DATA, st.dev = TRUE, which.model = 6, na.rm = FALSE)
auc7 <- auc(DATA, st.dev = TRUE, which.model = 7, na.rm = FALSE)
auc8 <- auc(DATA, st.dev = TRUE, which.model = 8, na.rm = FALSE)
auc9 <- auc(DATA, st.dev = TRUE, which.model = 9, na.rm = FALSE)

auc1[1,1]
# And compare the models with AUCs
dic  <- c(Mb1.c$dic$dic,
          Mb2.c$dic$dic,
          Mb3.c$dic$dic,
          Mb4.c$dic$dic,
          Mb5.c$dic$dic,
          Mb6.c$dic$dic,
          Mb7.c$dic$dic,
          Mb8.c$dic$dic,
          Mb9.c$dic$dic
          )   
waic <- c(Mb1.c$waic$waic,
          Mb2.c$waic$waic,
          Mb3.c$waic$waic,
          Mb4.c$waic$waic,
          Mb5.c$waic$waic,
          Mb6.c$waic$waic,
          Mb7.c$waic$waic,
          Mb8.c$waic$waic,
          Mb9.c$waic$waic
          )
auc  <- c(auc1[1,1], auc2[1,1],auc3[1,1],auc4[1,1],auc5[1,1],auc6[1,1],auc7[1,1],auc8[1,1],auc9[1,1]
          )
auc.sd  <- c(auc1[1,2], auc2[1,2],auc3[1,2],auc4[1,2],auc5[1,2],auc6[1,2],auc7[1,2],auc8[1,2],auc9[1,2]
             )

## more on CPO http://www.math.chalmers.se/~bodavid/GMRF2015/Lectures/Flab4.pdf
# use the logarithmic score calculated as -mean(log(inla.result$cpo))
logscore <- -mean(log(Mb15.c$cpo$cpo))
cpo.score  <- c(-mean(log(Mb1.c$cpo$cpo)),
                -mean(log(Mb2.c$cpo$cpo)),
                -mean(log(Mb3.c$cpo$cpo)),
                -mean(log(Mb4.c$cpo$cpo)),
                -mean(log(Mb5.c$cpo$cpo)),
                -mean(log(Mb6.c$cpo$cpo)),
                -mean(log(Mb7.c$cpo$cpo)),
                -mean(log(Mb8.c$cpo$cpo)),
                -mean(log(Mb9.c$cpo$cpo))
                )

Z.out.b     <- cbind(dic, waic, cpo.score, auc)
rownames(Z.out.b) <- c("linear baseline",
                       "linear baseline + SRF",
                       "linear + bio4 + distw + elev + Ca + Swater + SRF",
                       "linear + bio4 + distw + elev + Ca + Swater",
                       "linear + distw + elev + Ca + Swater",
                       "non-linear + f(distw) + f(elev) + f(Ca) + f(Swater)",
                       "non-linear + f(distw) + f(elev) + f(Ca) + f(Swater) + SRF",
                       "non-linear + distw + f(elev) + Ca + Swater",
                       "non-linear + distw + f(elev) + Ca + Swater + SRF"
                       )
Z.out.b
###########################


###########################
## Save in a data frame
Model <- c("linear baseline",
           "linear baseline + SRF",
           "linear + bio4 + distw + elev + Ca + Swater + SRF",
           "linear + bio4 + distw + elev + Ca + Swater",
           "linear + distw + elev + Ca + Swater",
           "non-linear + f(distw) + f(elev) + f(Ca) + f(Swater)",
           "non-linear + f(distw) + f(elev) + f(Ca) + f(Swater) + SRF",
           "non-linear + distw + f(elev) + Ca + Swater",
           "non-linear + distw + f(elev) + Ca + Swater + SRF")
model <- as.data.frame(Model)
dic <- as.data.frame(dic)
waic <- as.data.frame(waic)
auc <- as.data.frame(auc)
cpo.score <- as.data.frame(cpo.score)
Valid.results <- cbind(model,dic,waic,auc,cpo.score)

setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")
write.csv(Valid.results, "Model_Selection_100kmBuffer.csv")
###################################
