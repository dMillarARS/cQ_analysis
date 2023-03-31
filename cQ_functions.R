
#############
# FUNCTIONS #
#############

#=======================#
# (1) Separate Baseflow #
#=======================#

baseflowSeparation <- function(datetime,
                               streamflow, 
                               filterPara=0.99, 
                               passes=3){
  
  #---------#
  # Inputs: # 
  #---------#
  
  # dateTime = datetime data set for input streamflow measurements
  # streamflow = vector data set of stream flow measurements matching dateTime
  # filterPara = recursive filter coefficient
  # passes = number of passes for the filter
  
  #----------#
  # Outputs: # 
  #----------#
  
  # outputs = dataframe containing datetime, baseflow, stormflow, and original total discharge (stream flow)

  
  # Initialize previous pass's baseflow vector, starting with stream flow
  bfPreviousData <- streamflow
  
  # Initialize vector to store most current iteration of baseflow
  bfData <- bfPreviousData
  
  # Iterate over the number of filter passes
  for (i in 1:passes) {
    
    # Create flag for forward (flag = 1) or backward (flag = 0) pass
    f_or_b_flag <- i%%2
    
    if (f_or_b_flag == 1) {
      
      for (j in 2:length(bfPreviousData)) {
        
        # Equation 5 from Neira et al. 2020 (doi: 10.1016/j.jhydrol.2020.125559)
        updatedBF <- filterPara*bfData[j-1] + (1-filterPara)/2*(bfPreviousData[j] + bfPreviousData[j-1])
        bfData[j] <- ifelse(updatedBF < bfPreviousData[j], updatedBF, bfPreviousData[j]) 
        
      }
    }
    
    if (f_or_b_flag == 0) {
      
      for (j in (length(bfPreviousData)-1):1) {
        
        # Equation 5 from Neira et al. 2020 (doi: 10.1016/j.jhydrol.2020.125559)
        updatedBF <- filterPara*bfData[j+1] + (1-filterPara)/2*(bfPreviousData[j] + bfPreviousData[j+1])
        bfData[j] <- ifelse(updatedBF < bfPreviousData[j], updatedBF, bfPreviousData[j]) 
        
      }
      
    }  
    
    bfPreviousData = bfData
  }
  
  base_flow <- bfData
  storm_flow <- streamflow - base_flow
  total_flow <- streamflow
  
  # consolidate all hydrograph data into a single dataframe
  outputs <- data.frame(datetime, 
                        base_flow,
                        storm_flow,
                        total_flow)
  
  return(outputs)
  
}
  

#======================#
# (2) Smooth Stormflow #
#======================#

smoothStormFlow <- function(dateTime,
                            stormFlow,
                            totFlow,
                            passes=4) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # dateTime = datetime data set for input storm flow measurements
  # stormFlow = vector data set of storm flow measurements (time series data, datetime not included here)
  # passes = number of passes for the filter
  
  #----------#
  # Outputs: # 
  #----------#
  
  # output = dataframe containing datetime and smoothed stormflow
  
  
  smoothStormFlow <- stormFlow
  
  for(pass in 1:passes){
    
    for (q in 2:(length(smoothStormFlow)-1)) {
      
      # Equation 3 from Tang and Carey 2017 (doi: 10.1002/hyp.11185)
      smoothStormFlow[q] <- (smoothStormFlow[q-1] + smoothStormFlow[q]*2 + smoothStormFlow[q+1])/4
      
    }
  }
  
  
  output <- cbind.data.frame(datetime = dateTime,
                                smooth_st_flow = smoothStormFlow,
                                total_flow = totFlow)
  
  return(output)
  
}


#======================================================#
# (3) Identify storm events using storm flow threshold #
#======================================================#

baseFlowEventSep <- function(smoothStFlow,
                             sfThresh) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # smoothStFlow = smoothed storm flow dataframe generated with SmoothStormFlow function
  # sfThresh = discharge threshold for storm flow, above which constitutes a storm event
  
  #----------#
  # Outputs: # 
  #----------#
  
  # flowData = dataframe containing extra columns defining when storm events occured (storm_y_n)
  #               and numbering each event (storm_id)
  
  flowData <- smoothStFlow %>% mutate(storm_y_n = ifelse((smooth_st_flow > sfThresh), "yes", "no"))
  
  flowData$storm_id <- numeric(length=nrow(flowData))
  
  counter = 0
  
  for (dt in 2:nrow(flowData)) {
    
    if (flowData$storm_y_n[dt] == "yes" & flowData$storm_y_n[dt-1] == "no") {
      
      counter = counter + 1 
      
      flowData$storm_id[dt] = counter 
      
    }
    
    else if (flowData$storm_y_n[dt] == "yes" & flowData$storm_y_n[dt-1] == "yes"){
      
      flowData$storm_id[dt] = flowData$storm_id[dt-1]
      
    }
    
    else {flowData$storm_id[dt] = 0}
    
    
  }
  
  return(flowData)
  
}


#==========================================================================#
# (4) Process individual storm events: normalize discharge and constituent #
# concentration and break out into rising and falling limbs                #
#==========================================================================#

processStormEventsWithConc <- function(stormIds,
                                       conc,
                                       timestep_min,
                                       minDuration_hrs,
                                       maxDuration_hrs,
                                       filterParam,
                                       sfThresh) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # stormIds = dataframe with the full time series data set of stormflow
  #               generated using the 'BaseFlowEventSep' function
  
  # conc = dataframe with full time series data set of solute concentrations
  
  #----------#
  # Outputs: # 
  #----------#
  
  # output = a list containing 3 lists as follows:
  #
  #             - fullStorms - a list of individual storm events (one df per event) 
  #                               with observed and normalized discharge
  #
  #             - risingLimbs - a list of individual storm rising limbs-only (one df per event) 
  #                                with observed and normalized discharge
  #
  #             - fallingStorms - a list of individual storm falling limbs-only (one df per event) 
  #                               with observed and normalized discharge
  
  # merge the concentration data into the stormflow dataframe
  stormIds <- left_join(stormIds, conc, by = "datetime")
  
  # find the peak discharge measurement for each storm event
  peakFlowPerEvent <- stormIds %>% 
                        filter(storm_id != 0) %>% 
                          group_by(storm_id) %>% 
                            slice(which.max(smooth_st_flow))
  
  # find the start and end dates for each storm event
  startAndEndDTs <- stormIds %>%
                      filter(storm_id != 0) %>%
                        group_by(storm_id) %>%
                          summarize(start_dt=min(datetime),
                                    end_dt=max(datetime))
  
  # combine these event data into one dataframe
   allEventDTs <- left_join(peakFlowPerEvent,
                            startAndEndDTs,
                            by="storm_id") %>%
                                select(!c(storm_y_n)) %>%
                                  rename(peak_dt = datetime,
                                         peak_sm_flow = smooth_st_flow)
  
  # create empty lists to store event data
  fullStorms <- list()
  risingLimbs <- list()
  fallingLimbs <- list()
  
  # loop through the events and fill the respective dataframes with the appropriate data
  for (storm in 1:nrow(allEventDTs)) {
    
    stormID <- paste("storm_",as.character(storm),"", sep = "")
    
    start_date <- allEventDTs$start_dt[allEventDTs$storm_id == storm]
    
    peak_date <- allEventDTs$peak_dt[allEventDTs$storm_id == storm]
    
    end_date <- allEventDTs$end_dt[allEventDTs$storm_id == storm]
    
    currentStorm <- stormIds %>% filter(storm_id == storm)
    
    strmLngth_hrs <- timestep_min*nrow(currentStorm)/60
    
    if (strmLngth_hrs >= minDuration_hrs & strmLngth_hrs <= maxDuration_hrs) {
      
    stormFlow <- stormIds %>% 
      filter(storm_id == storm) 
    
    stormFlow$time_step = seq(1,nrow(stormFlow),1)
    
    stormFlow <- stormFlow %>%
                    select(-storm_y_n) %>%
                      mutate(norm_tot_q = (total_flow-min(total_flow))/(max(total_flow)-min(total_flow)),
                             norm_c = (conc-min(conc))/(max(conc)-min(conc)),
                             norm_ts = (time_step-min(time_step))/(max(time_step)-min(time_step)))
    
    stormFlow$filter_para <- as.character(filterParam)
    
    stormFlow$sf_thresh <- as.character(sfThresh)
  
    # populate fullStorms with data frames for all storm flow data
    fullStorms[[stormID]] = stormFlow
    
    # populate risingLimbs with dataframes for all rising limbs
    risingLimb <- stormFlow %>% 
                    filter(datetime <= peak_date)
    
    risingLimbs[[stormID]] <- risingLimb
    
    # populate risingLimbs with dataframes for all rising limbs
    fallingLimb <- stormFlow %>% 
                    filter(datetime > peak_date)
    
    fallingLimbs[[stormID]] <- fallingLimb
    
    }
  }
  
  # combine these 3 lists into a final list as output that can be
  # unpacked later
  
  output <- list()
  
  output[["fullStorms"]] <- fullStorms
  output[["risingLimbs"]] <- risingLimbs
  output[["fallingLimbs"]] <- fallingLimbs
  
  return(output)
  
  
}


#===========================================================================#
# (5) Batch run baseflow separation and event separation using              #
# different combinations of filter parameters and stormflow thresholds      #
# and get normalized discharge and nutrient concentrations for c-Q analyses #
#===========================================================================#

batchRunBfAndEvSepForCQ <- function(qInputs,
                                    bfSepPasses,
                                    filterParam,
                                    sfSmoothPasses,
                                    sfThresh,
                                    cInputs,
                                    timeStep,
                                    minDuration,
                                    maxDuration) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # qInputs = dataframe with the original time series discharge data
  
  # bfSepPasses = number of passes to use with recursive filter for bf separation
  
  # filterParam = vector containing filter coefficient values to use for bf separation
  
  # sfSmoothPasses = number of passes to use for the stormflow smoothing filter
  
  # sfThresh = vector containing stormflow threshold values (same flow units as qData)
  
  # cInputs = time series solute concentration data at same time steps as qData
  
  # timeStep = time step of the input data (minutes)
  
  # minDuration = the minimum storm duration to be stored in output (hours)
  
  # maxDuration = the maximum storm duration to be stored in output (hours)
  
  #----------#
  # Outputs: # 
  #----------#
  
  # batchOutputs = an list of length equal to: length(filterParam)*length(sfThresh),
  #                    where each element is its own list containing 3 lists as follows:
  #
  #             - fullStorms - a list of individual storm events (one df per event) 
  #                               with observed and normalized discharge and concentration
  #
  #             - risingLimbs - a list of individual storm rising limbs-only (one df per event) 
  #                                with observed and normalized discharge and concentration
  #
  #             - fallingStorms - a list of individual storm falling limbs-only (one df per event) 
  #                               with observed and normalized discharge and concentration
  
  
  # Create list to store output lists
  batchOutputs <- list()
  
  # Create dataframe with column for filter parameter values
  # and one for stormflow threshold values
  #-------------------------------------------------------------------
  filterCoefCol <- rep(filterParam,length(sfThresh))
  
  sfThreshCol <- numeric(0)
  
  for (i in 1:length(sfThresh)) {
    
    sfThreshCol <- c(sfThreshCol,rep(sfThresh[i],length(filterParam)))
    
  }
  
  inputParas <- bind_cols(filter_para = filterCoefCol,
                          sf_threshhold = sfThreshCol)
  #-------------------------------------------------------------------
  
  # Loop through the different parameter combinations
  for (i in 1:nrow(inputParas)){
    
    filterParaValue <- inputParas$filter_para[i]
    sfThreshValue <- inputParas$sf_threshhold[i]
    
    paraCombo <- paste("FC = ",as.character(filterParaValue),", SFT = ",as.character(sfThreshValue), sep = "")
    
    baseFlow <- baseflowSeparation(datetime = qInputs$datetime,
                                   streamflow = qInputs$q_cms,
                                   filterPara = filterParaValue,
                                   passes = bfSepPasses)
    
    smoothStorm <- smoothStormFlow(dateTime = baseFlow$datetime,
                                   stormFlow = baseFlow$storm_flow,
                                   totFlow = baseFlow$total_flow,
                                   passes = sfSmoothPasses)
    
    stormFlowWithIds <- baseFlowEventSep(smoothStFlow = smoothStorm,
                                         sfThresh = sfThreshValue)
    
    eventOutput <- processStormEventsWithConc(stormIds = stormFlowWithIds,
                                              conc = cInputs,
                                              timestep_min = timeStep,
                                              minDuration_hrs = minDuration,
                                              maxDuration_hrs = maxDuration,
                                              filterParam = filterParaValue,
                                              sfThresh = sfThreshValue)
    
    batchOutputs[[paraCombo]] <- eventOutput
    
  }
  
  return(batchOutputs)
}


#========================#
# (6) Get all storm data #
#========================#

getAllStormEvents <- function(batchRun,
                              timestep_min) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # batchRun = list created by running the 'batchRunBfAndEvSepForCQ' function
  
  #----------#
  # Outputs: # 
  #----------#
  
  # eventsData = dataframe containing all storm event data
  
  eventsData <- data.frame()
  
  for (i in 1:length(batchRun)) {
    
    if (length(batchRun[[i]]$fullStorms) > 0) {
      
      for (j in 1:length(batchRun[[i]][["fullStorms"]])) {
        
        event <- batchRun[[i]][["fullStorms"]][[j]]
        
        event$run_id <- names(batchRun1[i])
        event$storm_id <- names(batchRun1[[1]][["fullStorms"]][j])
        
        eventsData <- bind_rows(eventsData,event)
        
      }
      
    }
  }
  return(eventsData)
}


#==========================================================================#
# (7) Batch run baseflow separation using different filter parameters      #
#==========================================================================#


batchRunflowCompare <- function(qData,
                                bfSepPasses,
                                filterParam,
                                sfSmoothPasses) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # qData = dataframe with the original time series discharge data
  
  # bfSepPasses = number of passes to use with recursive filter for bf separation
  
  # filterParam = vector containing filter parameter values to use for bf separation
  
  # sfSmoothPasses = number of passes to use for the stormflow smoothing filter
  
  #----------#
  # Outputs: # 
  #----------#
  
  # batchOutputLF = long format dataframe with with storm flow for each filter parameter
  
  # Create list to store output lists
  batchOutputs <- list()
  
  # Loop through the different parameter combinations
  for (i in 1:length(filterParam)){
    
    filterParaValue <- filterParam[i]
    
    filterParaChr <- paste("filterPara_",as.character(filterParaValue),"_flows", sep = "")
    
    filterParaColLab <- paste("filterPara_",as.character(filterParaValue), sep = "")
    
    baseFlow <- baseflowSeparation(datetime = qData$date,
                                   streamflow = qData$q_cms,
                                   filterPara = filterParaValue,
                                   passes = bfSepPasses)
    
    smoothStorm <- smoothStormFlow(dateTime = baseFlow$datetime,
                                   stormFlow = baseFlow$storm_flow,
                                   totFlow = baseFlow$total_flow,
                                   passes = sfSmoothPasses)

    smoothStormQ <- smoothStorm %>% 
                      select(smooth_st_flow)
    
    allFlows <- bind_cols(baseFlow,smoothStormQ)
    
    allFlows$para_label <- filterParaColLab
    
    allFlows$filter_para <- as.character(filterParam[i])
    
    batchOutputs[[filterParaChr]] <- allFlows
    
  }
  
  batchOutputLF <- bind_rows(batchOutputs)
  
  return(batchOutputLF)
}


#===================================================================#
# (8) Function to count storm events for each parameter combination #
#===================================================================#

stormCounts <- function(batchRun) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # batchRun = list created by running the 'batchRunBfAndEvSepForCQ' function
  
  #----------#
  # Outputs: # 
  #----------#
  
  # stormCounts = dataframe containing storm counts based on each parameter combo
  
  totalRuns <- length(batchRun1)
  
  runNames <- character(totalRuns)
  stormCounts <- integer(totalRuns)
  
  stormCounts <- bind_cols(run_id = runNames,
                           storm_counts = stormCounts)
  
  for (i in 1:length(batchRun1)) {
    
    stormCounts$run_id[i] <- names(batchRun1[i])
    stormCounts$storm_counts[i] <- length(batchRun1[[i]][["fullStorms"]])
    
  }
  
  return(stormCounts)
  
}


#===================================================================#
# (9) Calculate total discharge, duration, and discharge intensity  #
#     for every storm event                                         #
#===================================================================#

stormEventCalcs <- function(batchRun,
                            timestep_min) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # batchRun = list created by running the 'batchRunBfAndEvSepForCQ' function
  
  #----------#
  # Outputs: # 
  #----------#
  
  # eventsData = dataframe containing storm event duration, total discharge, and
  #                 discharge intensity
  
  eventsData <- cbind.data.frame(run_id = character(0),
                                 start = character(0),
                                 end = character(0),
                                 tot_q_m3 = numeric(0),
                                 duration_hrs = numeric(0),
                                 intensity_m3_hr = numeric(0),
                                 filter_para = character(0),
                                 sf_thresh = character(0))
  
  for (i in 1:length(batchRun)) {
    
    if (length(batchRun[[i]]$fullStorms) > 0) {
      
      for (j in 1:length(batchRun[[i]][["fullStorms"]])) {
        
        run_id <- names(batchRun[i])
        storm_id <- paste("storm_",as.character(j),sep="")
        start <- as.character(min(batchRun[[i]][["fullStorms"]][[j]]$datetime))
        end <- as.character(max(batchRun[[i]][["fullStorms"]][[j]]$datetime))
        tot_q_m3 <- sum((batchRun[[i]][["fullStorms"]][[j]]$smooth_st_flow*60*timestep_min))
        duration_hrs <- timestep_min*nrow(batchRun[[i]][["fullStorms"]][[j]])/60
        intensity_m3_hr <- tot_q_m3/duration_hrs
        filter_para <- batchRun1[[i]][["fullStorms"]][[j]]$filter_para[1]
        sf_thresh <- batchRun[[i]][["fullStorms"]][[j]]$sf_thresh[1]
        
        eventRow <- cbind.data.frame(run_id,
                                     storm_id,
                                     start, 
                                     end, 
                                     tot_q_m3,
                                     duration_hrs,
                                     intensity_m3_hr,
                                     filter_para,
                                     sf_thresh)
        
        eventsData <- bind_rows(eventsData,eventRow)
        
      }
      
    }
  }
  return(eventsData)
}


#=======================================================================#
# (10) Get hysteresis indices (HI) by interpolating normalized c-Q data #
#=======================================================================#

getHysteresisIndices <- function(batchRun,
                                 xForInterp,
                                 eventsData) {
  
  #---------#
  # Inputs: # 
  #---------#
  
  # batchRun = list created by running the 'batchRunBfAndEvSepForCQ' function
  
  # xForInterp = vector with interpolation intervals used for calculating HI
  
  # eventsData = dataframe created by running stormEventCalcs
  
  #----------#
  # Outputs: # 
  #----------#
  
  # hysteresis = long format dataframe hysteresis index output for all runs
  
  
  # start an empty data frame to store hysteresis data
  hysteresisData <- cbind.data.frame(q_quant = numeric(0),
                                     risingerp = numeric(0),
                                     fallingerp = numeric(0),
                                     hyst_index = numeric(0),
                                     flsh_index = numeric(0),
                                     run_id = character(0),
                                     storm_id = character(0))
  
  for (i in 1:length(batchRun)){
    
    if (length(batchRun[[i]]$fullStorms) > 0) {
      
      for (j in 1:length(batchRun[[i]]$fullStorms)) {
        
        # sort out rising limb data
        risingLimbData <- batchRun[[i]]$risingLimbs[[j]]
        qNormRising <- risingLimbData$norm_tot_q
        cNormRising <- risingLimbData$norm_c
        
        # sort out falling limb data
        fallingLimbData <- batchRun[[i]]$fallingLimbs[[j]]
        qNormFalling <- fallingLimbData$norm_tot_q
        cNormFalling <- fallingLimbData$norm_c
        
        # make sure there are at least 2 points to interpolate per event
        lenRising <- length(cNormRising)
        lenFalling <- length(cNormFalling)
        if ((lenRising > 1) & 
            (lenFalling > 1) &
            (length(unique(qNormRising)) > 1) &
            (length(unique(qNormFalling)) > 1)) {
          
          interpRising <- approx(qNormRising,
                                 cNormRising,
                                 xout = xForInterp)
          
          interpFalling <- approx(qNormFalling,
                                  cNormFalling,
                                  xout = xForInterp)
          
          # combine data and calculate c-Q indices
          cQInterp <- cbind.data.frame(q_quant = xForInterp,
                                       risingerp = interpRising[[2]],
                                       fallingerp = interpFalling[[2]])
          
          #calculate hysteresis index 
          cQInterp$hyst_index <- cQInterp$risingerp - cQInterp$fallingerp
          
          # calculate flushing index
          cQInterp$flsh_index <- batchRun[[i]]$risingLimbs[[j]]$norm_c[length(batchRun[[i]]$risingLimbs[[j]]$norm_c)] -
            batchRun[[i]]$fullStorms[[j]]$norm_c[1]
          
          cQInterp$run_id <- names(batchRun[i])  
          
          cQInterp$storm_id <- paste("storm_",j,sep="")
          
          hysteresisData <- bind_rows(hysteresisData,cQInterp)
          
        }
      }
    }
  }
  
  hysteresisData <- left_join(hysteresisData,eventsData, by = c("run_id", "storm_id"))
  
  return(hysteresisData)
}


#====================================================#
# (11) Save hysteresis index vs flushing index plots #
#====================================================#

makeHystFlushPlotsNO3 <- function(hysteresisData) {
  
  # Make subfolder in output directory to save plots
  dir.create(file.path(output_dir, "HystFlushPlots"), showWarnings = FALSE)
  
  
  hysteresisData <- hysteresisData[complete.cases(hysteresisData),] 
  
  allRunIds <- unique(hysteresisData$run_id)
  
  for (i in 1:length(allRunIds)) {
    
    indexData <- hysteresisData %>% 
                    filter(run_id == allRunIds[i]) %>%
                      select(run_id,hyst_index,flsh_index,duration_hrs,intensity_m3_hr)
    
    hiFiPlot <- ggplot(indexData,aes(x=flsh_index,y=hyst_index,color=intensity_m3_hr,size=duration_hrs)) +
      geom_point(alpha=0.3) +
      scale_size_continuous(range = c(4,10)) +
      scale_color_viridis(option = "plasma")+
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      xlim(c(-1,1)) +
      ylim(c(-1,1)) +
      xlab("Flushing Index") +
      ylab("Hysteresis Index") +
      ggtitle(expression(paste("N",O[3]^{phantom("A")~"-"}))) +
      theme_bw() +
      theme(text=element_text(size=30))
    
    ggsave(file=file.path(output_dir, "HystFlushPlots",paste(allRunIds[i],"_",dataSetName,"_HiFi.jpeg",sep="")),
           hiFiPlot,
           width = 12, 
           height = 8, 
           units = "in",
           dpi=600)
  }
}

makeHystFlushPlotsTOC <- function(hysteresisData) {
  
  # Make subfolder in output directory to save plots
  dir.create(file.path(output_dir, "HystFlushPlots"), showWarnings = FALSE)
  
  
  hysteresisData <- hysteresisData[complete.cases(hysteresisData),] 
  
  allRunIds <- unique(hysteresisData$run_id)
  
  for (i in 1:length(allRunIds)) {
    
    indexData <- hysteresisData %>% 
                    filter(run_id == allRunIds[i]) %>%
                       select(run_id,hyst_index,flsh_index,duration_hrs,intensity_m3_hr)
    
    hiFiPlot <- ggplot(indexData,aes(x=flsh_index,y=hyst_index,color=intensity_m3_hr,size=duration_hrs)) +
      geom_point(alpha=0.3) +
      scale_size_continuous(range = c(4,10)) +
      scale_color_viridis(option = "plasma")+
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      xlim(c(-1,1)) +
      ylim(c(-1,1)) +
      xlab("Flushing Index") +
      ylab("Hysteresis Index") +
      ggtitle("TOC") +
      theme_bw() +
      theme(text=element_text(size=30))
    
    ggsave(file=file.path(output_dir, "HystFlushPlots",paste(allRunIds[i],"_",dataSetName,"_HiFi.jpeg",sep="")),
           hiFiPlot,
           width = 12, 
           height = 8, 
           units = "in",
           dpi=600)
  }
}

makeHystFlushPlotsTurb <- function(hysteresisData) {
  
  # Make subfolder in output directory to save plots
  dir.create(file.path(output_dir, "HystFlushPlots"), showWarnings = FALSE)
  
  
  hysteresisData <- hysteresisData[complete.cases(hysteresisData),] 
  
  allRunIds <- unique(hysteresisData$run_id)
  
  for (i in 1:length(allRunIds)) {
    
    indexData <- hysteresisData %>% 
                    filter(run_id == allRunIds[i]) %>%
                        select(run_id,hyst_index,flsh_index,duration_hrs,intensity_m3_hr)
    
    hiFiPlot <- ggplot(indexData,aes(x=flsh_index,y=hyst_index,color=intensity_m3_hr,size=duration_hrs)) +
      geom_point(alpha=0.3) +
      scale_size_continuous(range = c(4,10)) +
      scale_color_viridis(option = "plasma")+
      geom_hline(yintercept=0) +
      geom_vline(xintercept=0) +
      xlim(c(-1,1)) +
      ylim(c(-1,1)) +
      xlab("Flushing Index") +
      ylab("Hysteresis Index") +
      ggtitle("Turbidity") +
      theme_bw() +
      theme(text=element_text(size=30))
    
    ggsave(file=file.path(output_dir, "HystFlushPlots",paste(allRunIds[i],"_",dataSetName,"_HiFi.jpeg",sep="")),
           hiFiPlot,
           width = 12, 
           height = 8, 
           units = "in",
           dpi=600)
  }
}


#=====================#
# (12) Save c-Q plots #
#=====================#

makeCQPlotsNO3 <- function(batchRun) {
  
  # Make subfolder in output directory to save c-Q plots
  dir.create(file.path(output_dir, "CQPlots"), showWarnings = FALSE)
  
  runIds <- names(batchRun)
  
  for (i in 1:length(runIds)) {
    
    oneRun <- batchRun[[runIds[i]]]$fullStorms
    
    oneRun <- do.call(rbind, oneRun)
    
    cQPlot <- ggplot(oneRun,aes(x=norm_tot_q,y=norm_c,color=norm_ts,group=storm_id)) +
      geom_point(size=3) + 
      facet_wrap(~ storm_id, ncol = 4) +
      scale_color_gradientn(colours = c("darkblue","green","yellow","red"),
                            breaks=c(0,1),
                            labels=c("Event Start", "Event End")) +
      geom_path(color="black") +
      xlab("Normalized discharge") +
      scale_x_continuous(breaks=seq(0, 1, 0.5)) +
      ylab(expression(paste("Normalized ",NO[3]^{phantom("A")~"-"}," concentration"))) +
      theme_bw() +
      theme(legend.title = element_blank(),
            text=element_text(size=30),
            panel.spacing = unit(1, "lines"))
    
    ggsave(file=file.path(output_dir, "CQPlots",paste(runIds[i],"_",dataSetName,"CQ.jpeg",sep="")),
           cQPlot,
           width = 18, 
           height = 12, 
           units = "in",
           dpi=600)
  }
  
}

makeCQPlotsTOC <- function(batchRun) {
  
  # Make subfolder in output directory to save c-Q plots
  dir.create(file.path(output_dir, "CQPlots"), showWarnings = FALSE)
  
  runIds <- names(batchRun)
  
  for (i in 1:length(runIds)) {
    
    oneRun <- batchRun[[runIds[i]]]$fullStorms
    
    oneRun <- do.call(rbind, oneRun)
    
    cQPlot <- ggplot(oneRun,aes(x=norm_tot_q,y=norm_c,color=norm_ts,group=storm_id)) +
      geom_point(size=3) + 
      facet_wrap(~ storm_id, ncol = 4) +
      scale_color_gradientn(colours = c("darkblue","green","yellow","red"),
                            breaks=c(0,1),
                            labels=c("Event Start", "Event End")) +
      geom_path(color="black") +
      xlab("Normalized discharge") +
      scale_x_continuous(breaks=seq(0, 1, 0.5)) +
      ylab("Normalized total organic carbon") +
      theme_bw() +
      theme(legend.title = element_blank(),
            text=element_text(size=30),
            panel.spacing = unit(1, "lines"))
    
    ggsave(file=file.path(output_dir, "CQPlots",paste(runIds[i],"_",dataSetName,"CQ.jpeg",sep="")),
           cQPlot,
           width = 18, 
           height = 12, 
           units = "in",
           dpi=600)
  }
  
}

makeCQPlotsTurb <- function(batchRun) {
  
  # Make subfolder in output directory to save c-Q plots
  dir.create(file.path(output_dir, "CQPlots"), showWarnings = FALSE)
  
  runIds <- names(batchRun)
  
  for (i in 1:length(runIds)) {
    
    oneRun <- batchRun[[runIds[i]]]$fullStorms
    
    oneRun <- do.call(rbind, oneRun)
    
    cQPlot <- ggplot(oneRun,aes(x=norm_tot_q,y=norm_c,color=norm_ts,group=storm_id)) +
      geom_point(size=3) + 
      facet_wrap(~ storm_id, ncol = 4) +
      scale_color_gradientn(colours = c("darkblue","green","yellow","red"),
                            breaks=c(0,1),
                            labels=c("Event Start", "Event End")) +
      geom_path(color="black") +
      xlab("Normalized discharge") +
      scale_x_continuous(breaks=seq(0, 1, 0.5)) +
      ylab("Normalized turbidity") +
      theme_bw() +
      theme(legend.title = element_blank(),
            text=element_text(size=30),
            panel.spacing = unit(1, "lines"))
    
    ggsave(file=file.path(output_dir, "CQPlots",paste(runIds[i],"_",dataSetName,"CQ.jpeg",sep="")),
           cQPlot,
           width = 18, 
           height = 12, 
           units = "in",
           dpi=600)
  }
  
}
