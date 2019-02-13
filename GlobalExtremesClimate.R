# Workspace setup
# Install packages if not already installed
required.packages <- c("ggplot2", "raster", "sp", "rgdal", "plyr", "ncdf4")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(required.packages, new.packages)

library(plyr)
library(rgdal)
library(maptools)
library(raster)
library(ncdf4)
options(scipen = 999) # turn off scientific notation
par(mar=c(1,1,1,1))

soil.template <- raster("E:/Powell_GeoEvoEco/Extreme_environments/Global/SoilGrids1km/Depth_0cm/BLDFIE_M_sl1_1km_ll.tif")
continents <-readOGR("E:/Powell_GeoEvoEco/Extreme_environments/Global/Continents", "Continents")

# Read in master list of variables used to define and delineate extreme environments
All.vars <- read.csv("E:/Powell_GeoEvoEco/Extreme_environments/Extreme_final_for_R.csv")

# Keep only WorldClim variables that the group decided to use:
WC.rows <- All.vars[which(All.vars$Source=="WorldClim" & All.vars$Use.=="Y"), ]
WC.vars <- as.character(levels(WC.rows$Variable))[WC.rows$Variable]
WC.expl <- as.character(levels(WC.rows$Explanation))[WC.rows$Explanation]
WC.use.hi <- as.character(levels(WC.rows$High_values_are_extreme))[WC.rows$High_values_are_extreme]
WC.use.low <- as.character(levels(WC.rows$Low_values_are_extreme))[WC.rows$Low_values_are_extreme]


# WorldClim variables to use, taken from Candidate Extreme Variables v2.xlsx
WC.vars <- substr(WC.vars, 4, nchar(WC.vars))
WC.vars <- as.character(sprintf("%02d", as.integer(WC.vars)))


for (V in 1:length(WC.vars)){ # begin loop across WorldClim variables
  # read in world clim raster
  WC.grid <- raster(paste("E:/Powell_GeoEvoEco/Extreme_environments/Global/WorldClim_v2_bio_30s/wc2.0_bio_30s_", WC.vars[V], ".tif", sep=""))
  
  # reproject so it aligns with the SoilGrids1km grids
  print(paste("processing WorldClim var:", WC.vars[V]))
  WC.grid <- projectRaster(WC.grid, soil.template, method="bilinear")
  
  # loop through continents
  
  WC.extremes <- raster()
  store.thresholds.this.WC.var <- data.frame()
  
  for (C in 1:nrow(continents)){ # begin loop across continents
    cont.name <- continents[C,]$CONTINENT
    cont.name <- as.character(levels(cont.name))[cont.name]
    continent <- continents[C, "CONTINENT"]
    print(paste(".... processing", cont.name))
    WC.continent <- mask(crop(WC.grid, continent), continent)
  
    # Identify extremes for this continent
    
        # only high values are extreme, top 10%
        if(WC.use.hi[V]=="Y" & WC.use.low[V]=="N"){ 
          threshold.high <- unname(quantile(WC.continent, probs = 0.90, na.rm = TRUE))
          threshold.low <- NA
          WC.continent.extreme <- reclassify(WC.continent, rcl=c(-Inf, threshold.high, 0, threshold.high, Inf, 1))}
        
        # only low values are extreme, bottom 10%
        if(WC.use.hi[V]=="N" & WC.use.low[V]=="Y"){
          threshold.high <- NA
          threshold.low <- unname(quantile(WC.continent, probs = 0.10, na.rm = TRUE))
          WC.continent.extreme <- reclassify(WC.continent, rcl=c(-Inf, threshold.low, -1, threshold.low, Inf, 0))}
        
        # both directions are extreme, top 5% and bottom 5%
        if(WC.use.hi[V]=="Y" & WC.use.low[V]=="Y"){ 
          threshold.high <- unname(quantile(WC.continent, probs = 0.95, na.rm = TRUE))
          WC.continent.top.05 <- reclassify(WC.continent, rcl=c(-Inf, threshold.high, 0, threshold.high, Inf, 1))
          
          threshold.low <- unname(quantile(WC.continent, probs = 0.05, na.rm = TRUE))
          WC.continent.bottom.05 <- reclassify(WC.continent, rcl=c(-Inf, threshold.low, -1, threshold.low, Inf, 0))
          WC.continent.extreme <- WC.continent.top.05 + WC.continent.bottom.05}
    
        # store threshold values
        row <- data.frame(WC.var.number=WC.vars[V], WC.variable=WC.expl[V], cont.name, threshold.high, threshold.low)
        store.thresholds.this.WC.var <- rbind(store.thresholds.this.WC.var, row)
    
        # plot a map of climate extremes for this continent
        setwd("E:/Powell_GeoEvoEco/Extreme_environments/Global/Climate_outputs/PDF_maps")
        pdf(paste("BIO", WC.vars[V], cont.name, "extremes.pdf", sep="_"), width=11, height=8.5)
        plot(WC.continent.extreme, breaks=c(-5,-0.5, 0.5, 5), col=c("blue","lightgray","red"), 
             main=paste("WorldClim BIO-", WC.vars[V], ":", WC.expl[V], sep=""), legend=FALSE)
        mtext(paste("Low threshold = ", threshold.low, ";  High threshold = ", threshold.high, sep=""))
        dev.off()
        
        # Create histogram of climate variable for this continent
        setwd("E:/Powell_GeoEvoEco/Extreme_environments/Global/Climate_outputs/histograms")
        pdf(paste("BIO", WC.vars[V], cont.name, "histogram.pdf", sep="_"), width=11, height=8.5)
        hist(values(WC.continent), breaks=1000,
             main=cont.name, xlab=WC.expl[V])
        abline(v=threshold.low, col="blue", lty=2)
        abline(v=threshold.high, col="red", lty=2)
        dev.off()
        
    
    # Mosaic continents together
        if (C==1){
          WC.extremes <- WC.continent.extreme
        } else {
          WC.extremes <- mosaic(WC.extremes, WC.continent.extreme, fun="max")
        }
        
    # Clean up
      rm(WC.continent, WC.continent.extreme, WC.continent.top.05, WC.continent.extreme)
      gc()
  } # end loop across continents
  

# Plot outputs
setwd("E:/Powell_GeoEvoEco/Extreme_environments/Global/Climate_outputs")
pdf(paste("BIO", WC.vars[V], "All_Continents_extremes.pdf", sep="_"), width=11, height=8.5)
plot(WC.extremes, breaks=c(-5,-0.5, 0.5, 5), col=c("blue","lightgray","red"),
     main=paste("WorldClim BIO-", WC.vars[V], ":", WC.expl[V], sep=""), legend=FALSE)
plot(continents, add=TRUE)
dev.off()
    
# Save outputs
setwd("E:/Powell_GeoEvoEco/Extreme_environments/Global/Climate_outputs/rasters")
writeRaster(WC.extremes, paste("Climate_extreme_BIO_", WC.vars[V], ".tif", sep=""))

setwd("E:/Powell_GeoEvoEco/Extreme_environments/Global/Climate_outputs/thresholds")
write.csv(store.thresholds.this.WC.var, paste("thresholds_WCvar_", WC.vars[V], ".csv", sep=""))

# Clean up
rm(WC.extremes, WC.grid)
    
} # end loop across WorldClim variables
  

