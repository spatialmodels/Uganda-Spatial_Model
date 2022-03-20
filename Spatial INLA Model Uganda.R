#### Bernoulli GAM WITH SPATIAL CORRELATION
#############################################################################



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
mydata <- read.csv(file = "PresenceAbsence_75kmBuffer.csv",header = TRUE,na.strings = "NA",stringsAsFactors = TRUE,
                   dec = ".")
colnames(mydata)
mydata$longitude <- mydata$xcoord
mydata$latitude <- mydata$ycoord
####################################################################
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

# remove distance to water, carbon, evi, and soil pH,
MyVar <- c("bio4","bio12","elevation","soil_calcium","soilwater")
ke.vif <- corvif(mydata[,MyVar])

## Use Pearson correlations and scatterplots:
library(GGally)
ggpairs(mydata[,MyVar])







## 1.3: Relationships
#**********************
# Relationships between the response variable and covariates
MyVar <- c("bio4","bio12","elevation","distance_to_water","soil_calcium","soilwater")
MyMultipanel.ggp2(Z = mydata, 
                  varx = MyVar, 
                  vary = "outbreaks", 
                  ylab = "Total outbreaks",
                  addSmoother = TRUE,
                  addRegressionLine = FALSE,
                  addHorizontalLine = FALSE)


## Distibution of outbreak counts across locations?
#**********************
# 
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
# auc calculates the area under the ROC curve approximated with a Mann-Whitney U statistic, and
# (optionally) the associated standard deviation.

# auc(DATA, st.dev = TRUE, which.model = 1, na.rm = FALSE)
# Arguments:
# DATA a matrix or dataframe of observed and predicted values where 
# each row represents one plot and where columns are:
#  DATA[,1] plot ID text
#  DATA[,2] observed values zero-one values
#  DATA[,3] predicted probabilities from first model numeric (between 0 and 1)

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
write.csv(Valid.results, "Selection_and_validation_results_07Dec.csv")

#AUC ROC plot for best model to determine threshold for maximum sensitivity and specificity
DATA <- data.frame(
  ID = c(1:nrow(mydata)),
  observed = mydata$outbreaks01,
  predicted15 = Mb8.c$summary.fitted.values[1:nrow(mydata),"mean"]
)
auc.roc.plot(DATA, threshold=101, which.model=c(1), model.names=c("Best model"),
             na.rm=TRUE, xlab="1-Specificity (false positives)", ylab="Sensitivity (true positives)",
             main="ROC Plot", color=TRUE, line.type=TRUE, lwd=1, mark=0, mark.numbers=TRUE,
             opt.thresholds=TRUE, opt.methods=c(3), req.sens=0.85, req.spec=0.85, obs.prev=NULL,
             add.legend=TRUE, legend.text=NULL, legend.cex=0.8, add.opt.legend=TRUE, opt.legend.text=NULL,
             opt.legend.cex=0.7, counter.diagonal=TRUE, pch=NULL)


###################
library(ggregplot)
Efxplot(list(Mb8.c),ModelNames = c(""))
## Plot new fixed effects
Efxplot.Val(list(Mb8.c), ModelNames = c(""))

# Covariate effects
round(Mb8.c$summary.fixed[, c("mean", "0.025quant", "0.975quant")], 3)
###################


## Smoothers**** 
Ns <- nrow(mydata)
f.elevation    <- Mb8.c$summary.lincomb.derived[1:Ns + 0 * Ns, "mean"] 
SeLo.elevation <- Mb8.c$summary.lincomb.derived[1:Ns + 0 * Ns,"0.025quant"] 
SeUp.elevation <- Mb8.c$summary.lincomb.derived[1:Ns + 0 * Ns,"0.975quant"]

f.elevation <- exp(f.elevation)/(1 + exp(f.elevation))
SeLo.elevation <- exp(SeLo.elevation)/(1 + exp(SeLo.elevation))
SeUp.elevation <- exp(SeUp.elevation)/(1 + exp(SeUp.elevation))

Is.elevation <- order(mydata$elevation.std)

MyData <- data.frame(
  mu   = c(f.elevation[Is.elevation]), 
  SeUp = c(SeUp.elevation[Is.elevation]), 
  SeLo = c(SeLo.elevation[Is.elevation]), 
  Xaxis = c(sort(mydata$elevation)),
  ID    = factor(rep(c("Elevation"), 
                     each = nrow(mydata))))

# Plot
p <- ggplot()
p <- p + xlab("Covariate") + ylab("Probability of occurrence")
p <- p + theme(text = element_text(size = 15)) 
p <- p + geom_line(data = MyData,aes(x = Xaxis, y = mu), color="dodgerblue3", size = 0.6)
p <- p + geom_ribbon(data = MyData, aes(x = Xaxis,ymax = SeUp,ymin = SeLo),alpha = 0.2)
p <- p + facet_wrap(~ID, scales = "free", ncol = 1)  
p <- p + theme(panel.background = element_rect(fill = "grey96", colour = "grey92"))
p


############################################################
### Probability of occurrence curves for linear predictors 
Xm.p <- data.frame(
  Intercept              = rep(1, N),
  outbreaks              = mydata[,"outbreaks"],
  bio4                   = mydata[,"bio4"],
  distance_to_water      = mydata[,"distance_to_water"],
  elevation              = mydata[,"elevation"],
  soil_calcium           = mydata[,"soil_calcium"],
  soilwater              = mydata[,"soilwater"])


#### Bio4
MyData <- data.frame(bio4 = seq(from = min(Xm.p$bio4),to = max(Xm.p$bio4), length = 100))
Xmat <- model.matrix(~ bio4,data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)
Mb5.p <- inla(outbreaks ~ bio4,lincomb = lcb,
              control.predictor = list(link = 1,compute = TRUE),
              family = "binomial",Ntrials = 1,data = Xm.p)
Pred.pm <- Mb5.p$summary.lincomb.derived[,c("mean","0.025quant", "0.975quant")]
print(head(Mb5.p$summary.lincomb.derived), digits = 2)
Pred.marg <- Mb5.p$marginals.lincomb.derived
Pred.marg[[1]] 
MyFun <- function(x) {exp(x) / (1 + exp(x))}
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(MyFun , Pred.marg[[1]]))
inla.emarginal(MyFun, Pred.marg[[1]])
MyData$mu <- unlist(lapply(Pred.marg,function(x) inla.emarginal(MyFun, x)))
MyData$selo <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.025),inla.tmarginal(MyFun, x))))
MyData$seup <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.975),inla.tmarginal(MyFun, x))))
# Plot v.1:182   
p <- ggplot()
p <- p + geom_point(data = Xm.p,aes(y = outbreaks, x = bio4),shape = 1,size = 1)
p <- p + xlab("Temperature Seasonality (BIO4)") + ylab("Probability of occurrence")
p <- p + theme(text = element_text(size=18)) 
p <- p + geom_line(data = MyData,aes(x = bio4,y = mu),colour = "navyblue",size = 1.0)
p <- p + geom_ribbon(data = MyData,aes(x = bio4,ymax = seup,ymin = selo),alpha = 0.2)
p


#### Distance to water
MyData <- data.frame(distance_to_water = seq(from = min(Xm.p$distance_to_water),to   = 0.3, length = 100))
Xmat <- model.matrix(~ distance_to_water,data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)
Mb5.p <- inla(outbreaks ~ distance_to_water,lincomb = lcb,
              control.predictor = list(link = 1,compute = TRUE),
              family = "binomial",Ntrials = 1,data = Xm.p)
Pred.pm <- Mb5.p$summary.lincomb.derived[,c("mean","0.025quant", "0.975quant")]
print(head(Mb5.p$summary.lincomb.derived), digits = 2)
Pred.marg <- Mb5.p$marginals.lincomb.derived
Pred.marg[[1]] 
MyFun <- function(x) {exp(x) / (1 + exp(x))}
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(MyFun , Pred.marg[[1]]))
inla.emarginal(MyFun, Pred.marg[[1]])
MyData$mu <- unlist(lapply(Pred.marg,function(x) inla.emarginal(MyFun, x)))
MyData$selo <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.025),inla.tmarginal(MyFun, x))))
MyData$seup <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.975),inla.tmarginal(MyFun, x))))
# Plot v.1:182   
p <- ggplot()
p <- p + geom_point(data = Xm.p,aes(y = outbreaks, x = distance_to_water),shape = 1,size = 1)
p <- p + xlab("Distance to Water (km)") + ylab("Probability of occurrence")
p <- p + theme(text = element_text(size=18)) 
p <- p + geom_line(data = MyData,aes(x = distance_to_water,y = mu),colour = "navyblue",size = 1.0)
p <- p + geom_ribbon(data = MyData,aes(x = distance_to_water,ymax = seup,ymin = selo),alpha = 0.2)
p


## Elevation
MyData <- data.frame(elevation = seq(from = min(Xm.p$elevation),to = max(Xm.p$elevation),length = 100))
Xmat <- model.matrix(~ elevation,data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)
Mb5.p <- inla(outbreaks ~ elevation,lincomb = lcb,
              control.predictor = list(link = 1, compute = TRUE),
              family = "binomial",Ntrials = 1,data = Xm.p)
Pred.pm <- Mb5.p$summary.lincomb.derived[,c("mean","0.025quant", "0.975quant")]
print(head(Mb5.p$summary.lincomb.derived), digits = 2)
Pred.marg <- Mb5.p$marginals.lincomb.derived
Pred.marg[[1]] 
MyFun <- function(x) {exp(x) / (1 + exp(x))}
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(MyFun , Pred.marg[[1]]))
inla.emarginal(MyFun, Pred.marg[[1]])
MyData$mu <- unlist(lapply(Pred.marg,function(x) inla.emarginal(MyFun, x)))
MyData$selo <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.025),inla.tmarginal(MyFun, x))))
MyData$seup <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.975),inla.tmarginal(MyFun, x))))
# Plot v.1:182   
p <- ggplot()
p <- p + geom_point(data = Xm.p,aes(y = outbreaks, x = elevation),shape = 1,size = 1)
p <- p + xlab("Elevation (m)") + ylab("Probability of occurrence")
p <- p + theme(text = element_text(size=18)) 
p <- p + geom_line(data = MyData,aes(x = elevation,y = mu),colour = "navyblue",size = 1.0)
p <- p + geom_ribbon(data = MyData,aes(x = elevation,ymax = seup,ymin = selo),alpha = 0.2)
p


## soil calcium
MyData <- data.frame(soil_calcium = seq(from = min(Xm.p$soil_calcium),to = max(Xm.p$soil_calcium),length = 100))
Xmat <- model.matrix(~ soil_calcium,data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)
Mb5.p <- inla(outbreaks ~ soil_calcium,lincomb = lcb,
              control.predictor = list(link = 1,compute = TRUE),
              family = "binomial",Ntrials = 1,data = Xm.p)
Pred.pm <- Mb5.p$summary.lincomb.derived[,c("mean","0.025quant", "0.975quant")]
print(head(Mb5.p$summary.lincomb.derived), digits = 2)
Pred.marg <- Mb5.p$marginals.lincomb.derived
Pred.marg[[1]] 
MyFun <- function(x) {exp(x) / (1 + exp(x))}
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(MyFun , Pred.marg[[1]]))
inla.emarginal(MyFun, Pred.marg[[1]])
MyData$mu <- unlist(lapply(Pred.marg,function(x) inla.emarginal(MyFun, x)))
MyData$selo <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.025),inla.tmarginal(MyFun, x))))
MyData$seup <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.975),inla.tmarginal(MyFun, x))))
# Plot v.1:182   
p <- ggplot()
p <- p + geom_point(data = Xm.p,aes(y = outbreaks, x = soil_calcium),shape = 1,size = 1)
p <- p + xlab("Soil Calcium (cmol/kg)") + ylab("Probability of occurrence")
p <- p + theme(text = element_text(size=18)) 
p <- p + geom_line(data = MyData,aes(x = soil_calcium,y = mu),colour = "navyblue",size = 1.0)
p <- p + geom_ribbon(data = MyData,aes(x = soil_calcium,ymax = seup,ymin = selo),alpha = 0.2)
p


## soil water
MyData <- data.frame(soilwater = seq(from = min(Xm.p$soilwater),to = max(Xm.p$soilwater),length = 100))
Xmat <- model.matrix(~ soilwater,data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)
Mb5.p <- inla(outbreaks ~ soilwater,lincomb = lcb,
              control.predictor = list(link = 1,compute = TRUE),
              family = "binomial",Ntrials = 1,data = Xm.p)
Pred.pm <- Mb5.p$summary.lincomb.derived[,c("mean","0.025quant", "0.975quant")]
print(head(Mb5.p$summary.lincomb.derived), digits = 2)
Pred.marg <- Mb5.p$marginals.lincomb.derived
Pred.marg[[1]] 
MyFun <- function(x) {exp(x) / (1 + exp(x))}
inla.qmarginal(c(0.025, 0.5, 0.975), 
               inla.tmarginal(MyFun , Pred.marg[[1]]))
inla.emarginal(MyFun, Pred.marg[[1]])
MyData$mu <- unlist(lapply(Pred.marg,function(x) inla.emarginal(MyFun, x)))
MyData$selo <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.025),inla.tmarginal(MyFun, x))))
MyData$seup <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.975),inla.tmarginal(MyFun, x))))
# Plot v.1:182   
p <- ggplot()
p <- p + geom_point(data = Xm.p,aes(y = outbreaks, x = soilwater),shape = 1,size = 1)
p <- p + xlab("Soil Water (V%)") + ylab("Probability of occurrence")
p <- p + theme(text = element_text(size=18)) 
p <- p + geom_line(data = MyData,aes(x = soilwater,y = mu),colour = "navyblue",size = 1.0)
p <- p + geom_ribbon(data = MyData,aes(x = soilwater,ymax = seup,ymin = selo),alpha = 0.2)
p

#################################################################3


# HERE: predict and check the fit again:
# Fit of the model.
Pi <- Mb8.c$summary.fitted.values[1:N,"mean"]
Fit01 <- ifelse(Pi >= 0.52, 1, 0) # try different thresholds (0.6 best)
table(Fit01,mydata$outbreaks )
table(Fit01)
table(mydata$outbreaks)
Z <- table(mydata$outbreaks, Fit01)
Z







###########################################################################
### MODEL prediction ***********************

###########################################################################
# First, predict across Uganda then check the omission rate for wildlife data
# Calculate DIC and WAIC

#Set the working directory and load the data
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")

# Load prediction data
pred.data <- read.csv(file = "pred_data_1km_plusvars.csv",header = TRUE,na.strings = "NA",stringsAsFactors = TRUE,dec = ".")
pred.data$outbreaks <- rep(NA, nrow(pred.data))


# remove NAs from prediction data (save NAs as a separate file which you add to the final raster output)
summary(pred.data)
pred.nas1 <- pred.data[is.na(pred.data$soil_calcium),]
pred.nas2 <- pred.data[is.na(pred.data$soilwater),]
pred.nas <- rbind (pred.nas1, pred.nas2)

pred.data <- pred.data[!is.na(pred.data$soil_calcium),]
pred.data <- pred.data[!is.na(pred.data$soilwater),]
summary(pred.data)
dim(pred.data) # 210494

# Grid coords that were removed due to NAs in covariates
colnames(pred.nas)
pred.nas <- pred.nas[,c(3,4)]
head(pred.nas)
nas <- rep(NA,nrow(pred.nas))
pred.nas <- cbind(pred.nas, nas)
colnames(pred.nas) <- c("longitude", "latitude", "mean")
head(pred.nas)

list.pred <- list()
list.pred[[1]] <- pred.data[1:10000,]
list.pred[[2]] <- pred.data[10001:20000,]
list.pred[[3]] <- pred.data[20001:30000,]
list.pred[[4]] <- pred.data[30001:40000,]
list.pred[[5]] <- pred.data[40001:50000,]
list.pred[[6]] <- pred.data[50001:60000,]
list.pred[[7]] <- pred.data[60001:70000,]
list.pred[[8]] <- pred.data[70001:80000,]
list.pred[[9]] <- pred.data[80001:90000,]
list.pred[[10]] <- pred.data[90001:100000,]
list.pred[[11]] <- pred.data[100001:110000,]
list.pred[[12]] <- pred.data[110001:120000,]
list.pred[[13]] <- pred.data[120001:130000,]
list.pred[[14]] <- pred.data[130001:140000,]
list.pred[[15]] <- pred.data[140001:150000,]
list.pred[[16]] <- pred.data[150001:160000,]
list.pred[[17]] <- pred.data[160001:170000,]
list.pred[[18]] <- pred.data[170001:180000,]
list.pred[[19]] <- pred.data[180001:190000,]
list.pred[[20]] <- pred.data[190001:200000,]
list.pred[[21]] <- pred.data[200001:210494,]


list.pred.result <- list()

for(i in 1:21){
  pred.data1 <- list.pred[[i]]
  
  # Create columns for UTM in Km
  pred.data1$Xkm   <- pred.data1$xcoord #/ 1000 
  pred.data1$Ykm   <- pred.data1$ycoord #/ 1000
  
  # Standardize
  pred.data1$bio4c               <- MyStd(pred.data1$bio4)
  pred.data1$distance_to_waterc  <- MyStd(pred.data1$distance_to_water)
  pred.data1$elevationc          <- MyStd(pred.data1$elevation)
  pred.data1$soil_calciumc       <- MyStd(pred.data1$soil_calcium)
  pred.data1$soilwaterc          <- MyStd(pred.data1$soilwater)
  
  
  # Get the lego pieces for the  smoothers.
  Legoelevation.pred <- PredictMat(Legoelevation, pred.data1)
  
  
  # Give the basis functions unique names.
  colnames(Legoelevation.pred)    <- paste("elevation", 1:ncol_b, sep = "")
 
  
  
  Nx=nrow(pred.data1)
  Xm.pred <- data.frame(
    Intercept           = rep(1, Nx),
    distance_to_waterc  = pred.data1[,"distance_to_waterc"],
    elevationc  = pred.data1[,"elevationc"],
    soil_calciumc       = pred.data1[,"soil_calciumc"],
    soilwaterc          = pred.data1[,"soilwaterc"],
    elevation1 = Legoelevation.pred[,"elevation1"],
    elevation2 = Legoelevation.pred[,"elevation2"],
    elevation3 = Legoelevation.pred[,"elevation3"]
    )
  
  
  
  # Here is the stack for the observed data
  # Here is the stack
  N <- nrow(mydata)
  Stack <- inla.stack(
    tag = "Fit",
    data = list(outbreaks = mydata$outbreaks01),  
    A = list(1, 1, A),                 
    effects = list(  
      Intercept = rep(1, N),
      Xm        = Xm[,-1],
      w         = w.index))
  
  
  
  # The locations for which we want to do a prediction  
  # 
  
  LocPred <- pred.data1[,c("Xkm","Ykm")]
  LocPred <- as.matrix(LocPred)
  
  # Get the A matrix for this point:
  A.pred16 <- inla.spde.make.A(mesh = mesh, loc = LocPred)
  dim(A.pred16)
  # Instead of predicting for only 1 spatial point (which can be anywhere in
  # the mesh), you can also select a series of points for which you want to
  # predict. Just add them to LocPred. And these do not have to be actual 
  # sampling locations, but can be anywhere in the mesh!
  # Section 6.8.1 in Blangiardo and Cameletti does this for a large grid of 
  # spatial values, which are then plotted.
  
  # Here is the stack for the prediction Note that NA for NCalls.
  StackPred <- inla.stack(
    tag = "Predict",
    data = list(outbreaks = NA),  
    A = list(1, 1, A.pred16),               
    effects = list(  
      Intercept   = rep(1, Nx),
      Xm.pred     = Xm.pred[,-1],    
      w           = w.index))
  
  
  # We can combine the two stacks.         
  All.stacks <- inla.stack(Stack, StackPred)	              
  
  
  # ZT GAM with spatial correlation 
  # Check ommision rates 
  I1.Pred <- inla(outbreaks ~ -1 + Intercept + 
                    distance_to_waterc + soil_calciumc + soilwaterc +
                    elevation1 + elevation2 + elevation3,
                  family = "binomial",
                  data = inla.stack.data(All.stacks),
                  control.predictor = list(link = 1,
                                           A = inla.stack.A(All.stacks)))
  
  # Extract
  index.Fit  <- inla.stack.index(All.stacks, tag = "Fit")$data
  index.Pred <- inla.stack.index(All.stacks, tag = "Predict")$data
  # And we can extract the correct rows     
  Fit  <- I1.Pred$summary.fitted.values[index.Fit, c(1,2,3,5)]   #fitted values
  Pred <- I1.Pred$summary.fitted.values[index.Pred, c(1,2,3,5)]  #predicted values
  
  list.pred.result[[i]] <- Pred
}

Pred <- rbind(list.pred.result[[1]],
              list.pred.result[[2]],
              list.pred.result[[3]],
              list.pred.result[[4]],
              list.pred.result[[5]],
              list.pred.result[[6]],
              list.pred.result[[7]],
              list.pred.result[[8]],
              list.pred.result[[9]],
              list.pred.result[[10]],
              list.pred.result[[11]],
              list.pred.result[[12]],
              list.pred.result[[13]],
              list.pred.result[[14]],
              list.pred.result[[15]],
              list.pred.result[[16]],
              list.pred.result[[17]],
              list.pred.result[[18]],
              list.pred.result[[19]],
              list.pred.result[[20]],
              list.pred.result[[21]]
              
              )
pred.result <- cbind(pred.data, Pred) # combine original prediction data with the prediction results


#########################
### Convert to raster
head(pred.result)
summary(pred.result)
colnames(pred.result)

pred.result$longitude <- pred.result$xcoord
pred.result$latitude <- pred.result$ycoord
colnames(pred.result)

xyz.data.mean <- pred.result[,c(37,38,33)]
xyz.data.sd <- pred.result[,c(37,38,34)]
xyz.data.low <- pred.result[,c(37,38,35)]
xyz.data.up <- pred.result[,c(37,38,36)]


# Rbind results with locations that were removed due to NAs
xyz.data.mean <- rbind(xyz.data.mean, pred.nas)

sd.nas <- pred.nas
colnames(sd.nas) <- c("longitude", "latitude", "sd")
xyz.data.sd <- rbind(xyz.data.sd, sd.nas)

low.nas <- pred.nas
colnames(low.nas) <- c("longitude", "latitude", "0.025quant")
xyz.data.low <- rbind(xyz.data.low, low.nas)

up.nas <- pred.nas
colnames(up.nas) <- c("longitude", "latitude", "0.975quant")
xyz.data.up <- rbind(xyz.data.up, up.nas)

############################ 
## Irregular grid raster**
# set up an 'empty' raster, here via an extent object derived from your data
xyz.data <- as.matrix(xyz.data.mean)
colnames(xyz.data) <- c('X', 'Y', 'Z')
e <- extent(xyz.data[,1:2])
r <- raster(e, res = 1)
# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(xyz.data[, 1:2], r, xyz.data[,3], fun=mean)
crs(x) = "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"
plot(x, main = "Full model_mean")
### Save raster
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")
writeRaster(x, filename="Prediction_mean1km.tif", format="GTiff", overwrite=TRUE)

## sd
# set up an 'empty' raster, here via an extent object derived from your data
xyz.data <- as.matrix(xyz.data.sd)
colnames(xyz.data) <- c('X', 'Y', 'Z')
e <- extent(xyz.data[,1:2])
r <- raster(e, res = 1)
# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(xyz.data[, 1:2], r, xyz.data[,3], fun=mean)
crs(x) = "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"
plot(x, main = "Full model_sd")
### Save raster
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")
writeRaster(x, filename="Prediction_sd1km.tif", format="GTiff", overwrite=TRUE)


# up
xyz.data <- as.matrix(xyz.data.up)
colnames(xyz.data) <- c('X', 'Y', 'Z')
e <- extent(xyz.data[,1:2])
r <- raster(e, res = 1)
# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(xyz.data[, 1:2], r, xyz.data[,3], fun=mean)
crs(x) = "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"
plot(x, main ="Full model_up")
### Save raster
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")
writeRaster(x, filename="Prediction_up1km.tif", format="GTiff", overwrite=TRUE)

# low
xyz.data <- as.matrix(xyz.data.low)
colnames(xyz.data) <- c('X', 'Y', 'Z')
e <- extent(xyz.data[,1:2])
r <- raster(e, res = 1)
# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(xyz.data[, 1:2], r, xyz.data[,3], fun=mean)
crs(x) = "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"
plot(x, main = "Full model_low")
### Save raster
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")
writeRaster(x, filename="Prediction_low1km.tif", format="GTiff", overwrite=TRUE)


###########################################################################
# Extract the test data values from the predicted rasters.

## load test data (omission rate calculation)
#Set the working directory and load the data
setwd("~/CAMBRIDGE/6. Data_analysis/INLA Uganda/species")

test.data <- readOGR("NorthEastpoints_UTM36N.shp")
crs(test.data) <- "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"


# set up an 'empty' raster, here via an extent object derived from your data
xyz.data <- as.matrix(xyz.data.mean)
colnames(xyz.data) <- c('X', 'Y', 'Z')
e <- extent(xyz.data[,1:2])
r <- raster(e, res = 1)
# you need to provide a function 'fun' for when there are multiple points per cell
x <- rasterize(xyz.data[, 1:2], r, xyz.data[,3], fun=mean)
crs(x) = "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"
plot(x, main = "Full model_mean")

# WGS84 points
crs(x) = "+proj=utm +zone=36N +datum=WGS84 +units=km +ellps=WGS84 +towgs84=0,0,0"
points.d <- mydata
coordinates(points.d)= ~ longitude + latitude
plot(x, main = "Full model")
plot(test.data, col= "blue", add=T)
points(points.d)


# Calculate Omission
# Extract raster values for original data
## load test data
rasValue = raster::extract(x, test.data)
rasValue = as.data.frame(rasValue)
rasValue
omm.I1 <- length(rasValue[(rasValue$rasValue<0.52),])/nrow(rasValue)
omm.I1 # [1] [1] 0.1039604

#####################




