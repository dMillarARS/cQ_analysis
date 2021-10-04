
# Main script to separate baseflow, delineate storm events, 
# and run c-Q hysteresis index analyses

# Clear memory
rm(list=ls())

#################
# LOAD PACKAGES #
#################

library(dataRetrieval)
library(tidyverse)
library(viridis)

###################
# SET DIRECTORIES #
###################

# Define the input directory
input_dir <- "C:\\...."
output_dir <- "C:\\...."

#####################
# READ IN FUNCTIONS #
#####################

source(file.path(input_dir,"cQ_functions.R"))



# Download data from USGS...
#-------------------------------------------------------------
# Difficult Run
siteNumber <- "01646000"

siteInfo <- readNWISsite(siteNumber)

parametersOfInterest <- c("00060","99133")

# Read in the data
dataset <- readNWISuv(siteNumber,parametersOfInterest,
                      "2020-01-01","2021-01-01", tz="EST")


allInputData15Min <- dataset %>% 
  select(dateTime,
         X_00060_00000,
         X_99133_00000) %>%
  rename(datetime = dateTime,
         q_cfs = X_00060_00000,
         conc = X_99133_00000) %>%
  mutate(q_cms = q_cfs * 0.0283168) # convert ft^3 s^-1 to m^3 s^-1
#-------------------------------------------------------------

# Subset the data to the longest time period without missing concentration data
#---------------------------------

allInputData15Min$section_id <- numeric(length=nrow(allInputData15Min))

counter = 0

for (dt in 2:nrow(allInputData15Min)) {
  
  if (is.na(allInputData15Min$conc[dt]) == FALSE & is.na(allInputData15Min$conc[dt-1]) == TRUE) {
    
    counter = counter + 1 
    
    allInputData15Min$section_id[dt] = counter 
    
  }
  
  else if (is.na(allInputData15Min$conc[dt]) == FALSE & is.na(allInputData15Min$conc[dt-1]) == FALSE){
    
    allInputData15Min$section_id[dt] = allInputData15Min$section_id[dt-1]
    
  }
  
  else {allInputData15Min$section_id[dt] = 0}
  
  
}

LongestStretchID <-  allInputData15Min[complete.cases(allInputData15Min),] %>% 
                        count(section_id) %>% 
                            filter(n == max(n)) %>% 
                                select(section_id)

allInputData15Min <- allInputData15Min %>% filter(section_id == as.numeric(LongestStretchID))

#---------------------------------

# Name for associated output data
dataSetName <- "Difficult_Run"

# Chose constituent for plot axes labels (NO3, TOC, or turbidity)
constit <- "NO3"

# Interpolate missing flow data (only a few cases where flow was 0 for one contiguous time step)...
allInputData15Min$q_cms[allInputData15Min$q_cms == 0] <- NA
interpFlow <- approx(allInputData15Min$datetime,
                     allInputData15Min$q_cms,
                     xout=allInputData15Min$datetime)

allInputData15Min$q_cms <- interpFlow$y

allInputData15Min <- allInputData15Min %>% 
                        mutate(rescaled_conc = ((conc-min(conc))/(max(conc)-min(conc))*max(q_cms)))

# Vector containing candidate baseflow separation filter values
candidateFilterPara <- c(0.96,0.98,0.99)

# Vector containing candidate stormflow threshold values
candidateSfThresh <- c(0.05,0.1,0.5)

# Vector with interpolation intervals used for calculating HI
interp <- seq(0,1,0.01)


##########################################
# RUN ANALYSIS TO GET HYSTERESIS INDICES #
##########################################

batchRun1 <- batchRunBfAndEvSepForCQ(qInputs = allInputData15Min,
                                        bfSepPasses = 3,
                                        filterParam = candidateFilterPara,
                                        sfSmoothPasses = 4,
                                        sfThresh = candidateSfThresh,
                                        cInputs = allInputData15Min,
                                        timeStep = 15,
                                        minDuration = 12,
                                        maxDuration = 8000)

eventsDataAll1 <- getAllStormEvents(batchRun = batchRun1,
                                    timestep_min = 15)

batchRunFlowsLF1 <- batchRunflowCompare(qData = allInputData15Min,
                                         bfSepPasses = 3,
                                         filterPara = candidateFilterPara,
                                         sfSmoothPasses = 4)

eventsData1 <- stormEventCalcs(batchRun = batchRun1,
                               timestep_min = 15)


stormCounts1 <- stormCounts(batchRun1)


hysteresisData1 <- getHysteresisIndices(batchRun = batchRun1,
                                        xForInterp = interp,
                                        eventsData = eventsData1)

#---------------------------------

######################
# EXPORT OUTPUT DATA #
######################

write.csv(eventsData1,file = file.path(output_dir,paste(dataSetName,"_StormEventSummaryData.csv",sep="")))
write.csv(batchRunFlowsLF1,file = file.path(output_dir,paste(dataSetName,"_DischargeData.csv",sep="")))
write.csv(hysteresisData1,file = file.path(output_dir,paste(dataSetName,"_HysteresisData.csv",sep="")))
write.csv(eventsDataAll1,file = file.path(output_dir,paste(dataSetName,"_AllCQData.csv",sep="")))
write.csv(stormCounts1,file = file.path(output_dir,paste(dataSetName,"_StormCounts.csv",sep="")))


#########################################
# PLOT AND SAVE DATA - EVENT SEPARATION #
#########################################

# Make subfolder in output directory to save hydrograph plots
dir.create(file.path(output_dir, "Hydrographs"), showWarnings = FALSE)

# 1) Plot and save the hydrograph with input data

initialHydrograph <- ggplot(allInputData15Min,aes(x=datetime, y=q_cms)) +
                            geom_line(size=0.5, color="black") +
                            xlab(NULL) +
                            ylab(expression(paste("Total discharge (",m^3," ",s^-1,")"))) +
                            theme_bw() +
                            theme(text=element_text(size=18))

ggsave(file=file.path(output_dir,"Hydrographs",paste(dataSetName,"_TotalDischarge.jpeg")),
       initialHydrograph,
       width = 12, 
       height = 4, 
       units = "in",
       dpi=600)


# 2) Plot total discharge with baseflow

baseflowHydrograph <- ggplot() + 
                            geom_line(data=batchRunFlowsLF1, aes(x=datetime, y=total_flow), size=0.5, color="black") +
                            geom_line(data=batchRunFlowsLF1, aes(x=datetime, y=base_flow,color=filter_para), size=0.75) +
                            scale_color_brewer(palette = "Set1") +
                            xlab(NULL) +
                            ylab(expression(paste("Discharge (",m^3," ",s^-1,")"))) +
                            theme_bw() +
                            theme(legend.title = element_blank(),
                                  text=element_text(size=18))

ggsave(file=file.path(output_dir,"Hydrographs",paste(dataSetName,"_Baseflows.jpeg")),
       baseflowHydrograph,
       width = 14, 
       height = 4, 
       units = "in",
       dpi=600)


# 3) Plot smoothed storm flows

stormflowHydrograph <- ggplot() + 
  geom_line(data=batchRunFlowsLF1, aes(x=datetime, y=storm_flow,color=filter_para), size=0.75) +
  scale_color_brewer(palette = "Set1") +
  xlab(NULL) +
  ylab(expression(paste("Storm flow (",m^3," ",s^-1,")"))) +
  theme_bw() +
  theme(legend.title = element_blank(),
        text=element_text(size=18))

ggsave(file=file.path(output_dir,"Hydrographs",paste(dataSetName,"_StormflowsOnly.jpeg")),
       stormflowHydrograph,
       width = 14, 
       height = 4, 
       units = "in",
       dpi=600)


# 3a) Plot smoothed storm flows with storm flow thresholds

stormflowThreshHydrograph <- ggplot() + 
  geom_line(data=batchRunFlowsLF1, aes(x=datetime, y=storm_flow,color=filter_para), size=0.75) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = candidateSfThresh, linetype = "dashed", color = "black",alpha=0.5) +
  xlab(NULL) +
  ylab(expression(paste("Storm flow (",m^3," ",s^-1,")"))) +
  theme_bw() +
  theme(legend.title = element_blank(),
        text=element_text(size=18))

ggsave(file=file.path(output_dir,"Hydrographs",paste(dataSetName,"_StormflowsOnlyWithThresholds.jpeg")),
       stormflowThreshHydrograph,
       width = 14, 
       height = 4, 
       units = "in",
       dpi=600)


# 4) Plot batch run event separation hydrographs
eventsDataShaded1 <- eventsData1 %>% mutate(start = as.POSIXct(start,
                                                                   format("%Y-%m-%d %H:%M:%S"),tz="EST"),
                                            end = as.POSIXct(end,
                                                                   format("%Y-%m-%d %H:%M:%S"),tz="EST"),
                                            tops = max(allInputData15Min$q_cms),
                                            bottoms = 0)

batchEventSepPlot <- ggplot() + 
  geom_rect(data=eventsDataShaded1, mapping=aes(xmin=start, 
                                                xmax=end, 
                                                ymin=bottoms, 
                                                ymax=tops), fill="green", color="red", alpha=0.2) +
  
  geom_line(data=allInputData15Min, aes(x=datetime, y=q_cms), size=0.8, color="blue") +
  geom_line(data=allInputData15Min, aes(x=datetime, y=rescaled_conc), size=0.5, color="black",linetype="dashed") +
  facet_wrap(~ run_id, ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  xlab(NULL) +
  ylab(expression(paste("Discharge (",m^3," ",s^-1,")"))) +
  theme_bw() +
  theme(legend.title = element_blank(),
        text=element_text(size=18))


ggsave(file=file.path(output_dir,"Hydrographs",paste(dataSetName,"_BatchEventSeparationPlot.jpeg")),
       batchEventSepPlot,
       width = 14, 
       height = 10, 
       units = "in",
       dpi=600)

####################################
# PLOT AND SAVE DATA - c-Q RESULTS #
####################################

if (constit == "NO3") {
  
  makeCQPlotsNO3(batchRun1)
  makeHystFlushPlotsNO3(hysteresisData1)
  
} else if (constit == "TOC") {
  
  makeCQPlotsTOC(batchRun1)
  makeHystFlushPlotsTOC(hysteresisData1)
  
} else if (constit == "turbity") {
  
  makeCQPlotsTurb(batchRun1) 
  makeHystFlushPlotsTurb(hysteresisData1)
  
}
