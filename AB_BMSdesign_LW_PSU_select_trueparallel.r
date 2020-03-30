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




setwd("/home/jpyanez/mader/ab_lightweight")

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
aea.proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-96+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
laea.proj <- "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs" # use this because it is lccc10 projection
tmerc.proj <- "+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +datum=NAD83 +units=m +no_defs" #Convert to AB hexagons projection
NAD.proj <- "+proj=utm +zone=14 +ellps=GRS80 +datum=NAD83 +units=m +no_defs " #Alberta hexagons are in this projection

sample_allocation<-read.csv("Sample_Size_Allocations.csv")
sample_allocation<-subset(sample_allocation, POL_DIV=="Alberta")
samplesize<-as.vector(sample_allocation$sample.size)  #  #### 


#
REGION_NAM<-as.vector(sample_allocation$ECOREGION_ID)
sample.size<-data.frame(samplesize, REGION_NAM)
sample.size

eco.list <-REGION_NAM




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
                  "Athabasca Plain"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Athabasca Plain")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Athabasca Plain")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Boreal Transition"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Boreal Transition")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Boreal Transition")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Clear Hills Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Clear Hills Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Clear Hills Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Eastern Continental Ranges"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Eastern Continental Ranges")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Eastern Continental Ranges")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Hay River Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Hay River Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Hay River Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Mid-Boreal Uplands"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Mid-Boreal Uplands")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Mid-Boreal Uplands")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Northern Alberta Uplands"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Northern Alberta Uplands")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Northern Alberta Uplands")$samplesize)*oversample.size),0,seltype="Continuous"),
                  #"Northern Continental Divide"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Northern Continental Divide")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Northern Continental Divide")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Peace Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Peace Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Peace Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Slave River Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Slave River Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Slave River Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Tazin Lake Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Tazin Lake Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Tazin Lake Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Wabasca Lowland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Wabasca Lowland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Wabasca Lowland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Western Alberta Upland"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Western Alberta Upland")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Western Alberta Upland")$samplesize)*oversample.size),0,seltype="Continuous"),
                  "Western Boreal"=list(panel=c(PanelOne=subset(sample.size, REGION_NAM=="Western Boreal")$samplesize), over=ceiling((subset(sample.size, REGION_NAM=="Western Boreal")$samplesize)*oversample.size),0,seltype="Continuous")
                  
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
        file= paste0("/home/jpyanez/mader/ab_lightweight/output_data/",outname))
#full.design.samp