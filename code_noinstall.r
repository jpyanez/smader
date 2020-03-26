.libPaths("/home/jpyanez/mader/mb_lightweight/rpackages")
#install.packages("MBHdesign", repos='http://cran.us.r-project.org')
#install.packages("spsurvey", repos='http://cran.us.r-project.org')
#install.packages("maptools", repos='http://cran.us.r-project.org')
#install.packages("rgeos", repos='http://cran.us.r-project.org')
#install.packages("raster", repos='http://cran.us.r-project.org')
#install.packages("car", repos='http://cran.us.r-project.org')
#install.packages("sf", repos='http://cran.us.r-project.org')
#install.packages("tidyverse", repos='http://cran.us.r-project.org')
#install.packages("sp", repos='http://cran.us.r-project.org')
#install.packages("rgdal", repos='http://cran.us.r-project.org')
#install.packages("lakemorpho", repos='http://cran.us.r-project.org')
#install.packages("dplyr", repos='http://cran.us.r-project.org')
#install.packages("car", repos='http://cran.us.r-project.org')
#install.packages("doBy", repos='http://cran.us.r-project.org')
#install.packages("gdata", repos='http://cran.us.r-project.org')
#install.packages("ggplot2", repos='http://cran.us.r-project.org')
#install.packages("nngeo", repos='http://cran.us.r-project.org')
#install.packages("reshape2", repos='http://cran.us.r-project.org')
#install.packages("viridis", repos='http://cran.us.r-project.org')
#install.packages("scales", repos='http://cran.us.r-project.org')
#install.packages("parallel", repos='http://cran.us.r-project.org')
#install.packages("doSNOW", repos='http://cran.us.r-project.org')
#install.packages("foreach", repos='http://cran.us.r-project.org')
#install.packages("vegan", repos='http://cran.us.r-project.org')

library(parallel) 
Ncores <- 12 #(detectCores())-2 #Number of available CPU cores. will need this number later on.
print("Number of cores")
print(Ncores)
#will need this number later on, but the parallel package will interfere with the SNOW package
#but the parallel package will interfere with some other packages, so we will remove it now.
detach("package:parallel", unload = TRUE)


library(MBHdesign)
library(spsurvey)
library(maptools)
library(rgeos)
library(raster)
library(car)
library(sf)
library(tidyverse)
library(sp)
library(rgdal)
library(lakemorpho)
library(dplyr)
library(car)
library(doBy)
library(gdata)
library(ggplot2)
library(nngeo)
library(reshape2)
library(viridis)
library(scales)
library(stringi)




setwd("/home/jpyanez/mader/mb_lightweight")

### Setting outname ###
set.seed(sample(1:1000000,1))
outname=paste0("data","_",stri_rand_strings(1, 10, pattern = "[A-Za-z0-9]"),".Rdata")
print("OUTNAME")
print(outname)

#load("C:\\Users\\vwilgs\\Documents\\BOREAL POSITION\\Projects\\Boreal Monitoring Strategy\\BMS design test\\Test_Design_largeN.RData")
#setwd("C:/Users/vwilgs/Documents/BOREAL POSITION/Projects/Saskatchewan Breeding Bird Atlas/Boreal Site Selection/")



###################################################################################################################################
##########################################          IMPORT THE SAMPLING FRAME            ##########################################
###################################################################################################################################


# Albers equal-area projection parameters based on above

laea.proj <- "+proj=laea +lat_0=0 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" # use this because it is lccc10 projection
NAD.proj <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs " #Manitoba hexagons are in this projection


#samplesize<-c(57,20,34,67,12,11,1,11)  # target sample sizes for priority squares
sample_allocation<-read.csv("Sample_Size_Allocations.csv")
sample_allocation<-subset(sample_allocation, POL_DIV=="Manitoba")
samplesize<-as.vector(sample_allocation$sample.size)  #  #### 


#
REGION_NAM<-as.vector(sample_allocation$ECOREGION_ID)
sample.size<-data.frame(samplesize, REGION_NAM)
sample.size

eco.list <-REGION_NAM




###### NOTE: there are previous study locations  (3 hexagons) in Riding Mountain National Park that meet our definition of "legacy sites", so we will adjust sample sizes accordingly
sample.size$samplesize[sample.size$REGION_NAM == "Mid-Boreal Uplands"] <- 15- 12  #### original sample size was 15, but here I have reduced this to account for existing legacy sites
sample.size$samplesize[sample.size$REGION_NAM == "Mid-Boreal Lowland"] <- 56- 3  #### 
sample.size$samplesize[sample.size$REGION_NAM == "Selwyn Lake Upland"] <- 77  ####  Don't adjust sample size as these are have centroid in SK, just want to modify inclusion probabilities near them
sample.size$samplesize[sample.size$REGION_NAM == "Interlake Plain"] <- 36 -1  ####  


#sub<-subset(hexagons, LegacySite==1)
#table(sub$REGION_NAM)
#sample.size$Legacy<- c(50,11,26,42,13,7,1,15)
#colnames(sample.size)<-c("Target","REGION_NAM","Legacy")
#sample.size$samplesize<-sample.size$Target - sample.size$Legacy
#sample.size$samplesize<-ifelse(sample.size$samplesize< 1, 1,sample.size$samplesize)



#############################################################################################################################
###                                      
###                                      
###                                      		Select the REGION_NAM to draw from
###                                      
###                                      
#############################################################################################################################




#############################################################################################################################
###                                      
###                                      
###                                      TEST THE COMPETING DESIGNS by drawing repeated random samples
###                                      
###                                      
#############################################################################################################################


# set.seed(3911548)  # Set seed 
# create mdcaty.  keeps variable it is based on in output in addition to mdcaty.  useful sometimes
# also want sum of inclusion probabilities to sum to expected sample size.
#attframe <- read.dbf("Sample_Frame")
attframe <- read_sf("sample_frame.shp") 




###  SPECIFY THE STRATIFIED SAMPLING DESIGN
oversample.size = 0.2 ### 2x oversample


################### RESTART HERE
# DESIGN 3: FULL BMS DESIGN -----------------------------------------------------
print("Got to D3")
#attframe$mdcaty <- N * attframe$p/sum(attframe$p)
attframe$stratum <- attframe$REGION_NAM

#### STRATIFIED DESIGN WITH THE USE OF INCLUSION PROBABILITIES
# NO OVERSAMPLE because just testing the design with Panel one
Stratdsgn <- list("Aspen Parkland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Aspen Parkland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Aspen Parkland")$samplesize)*oversample.size),0, seltype="Continuous"),
                  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Boreal Transition")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Boreal Transition")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Coastal Hudson Bay Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Coastal Hudson Bay Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Coastal Hudson Bay Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Churchill River Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Churchill River Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Aspen Parkland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Hudson Bay Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Hudson Bay Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Hudson Bay Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Hayes River Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Hayes River Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Hayes River Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Interlake Plain"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Interlake Plain")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Interlake Plain")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Kazan River Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Kazan River Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Kazan River Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Lake Manitoba Plain"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Lake Manitoba Plain")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Lake Manitoba Plain")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Lac Seul Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Lac Seul Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Lac Seul Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Lake of the Woods"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Lake of the Woods")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Lake of the Woods")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Mid-Boreal Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Mid-Boreal Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Mid-Boreal Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Mid-Boreal Uplands")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Mid-Boreal Uplands")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Maguse River Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Maguse River Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Maguse River Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Selwyn Lake Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Selwyn Lake Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Selwyn Lake Upland")$samplesize)*oversample.size),0,seltype="Continuous")
)
Stratdsgn



n.simulations=1 #### run 100 random draws
full.design.samp<-list() ###

########################################################################################
####     Added Feb 2020: Initiate parallel processing, so draws don't take 3 days   ####
########################################################################################

library(doSNOW) #must use doSNOW instead of parallel because windows doesn't support multiple cores with parallel for some stupid reason
library(foreach)
#print("Made it after foreach")
cl <- makeCluster(Ncores) #number of CPU cores to use in cluster
registerDoSNOW(cl) #registers cluster on parallel backend
clusterEvalQ(cl, library("spsurvey"))
variables <- c("Stratdsgn", "attframe","sample.size","oversample.size")
clusterExport(cl, variables)

start.time <- Sys.time()
#### this differs because it saves every iteration, other designs ONLY saving the last iteration
pb <-txtProgressBar(max=500, style = 3)
progress <-function(n) setTxtProgressBar(pb,n)
opts<-list(progress =progress)
full.design.samp <- foreach (i = 1:n.simulations, .options.snow =opts)  %dopar% {
  full.design.samp[i] <- grts(design=Stratdsgn,
                              DesignID="ET_ID",
                              type.frame="finite",
                              src.frame="att.frame",
                              in.shape="Sample_Frame",
                              att.frame=attframe,
                              xcoord = "X",
                              ycoord = "Y",
                              stratum="REGION_NAM",
                              mdcaty="mdcaty",
                              shapefile=FALSE,
                              out.shape=paste0("Full_BMS_design", i)) 
}

end.time <-Sys.time()
draw.time<-end.time-start.time
draw.time
stopCluster(cl)
registerDoSEQ()

#set.seed(sample(1:1000000,1))
#outname=paste0("data","_",stri_rand_strings(1, 10, pattern = "[A-Za-z0-9]"),".Rdata")
saveRDS(full.design.samp,
        file= paste0("/home/jpyanez/mader/mb_lightweight/output_data/",outname))
#full.design.samp