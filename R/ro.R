ro <- function(data,h=10,origins=10,call,value=NULL,
               ci=FALSE,co=FALSE,silent=FALSE,parallel=FALSE){
# Function makes Rolling Origin for the data using the call
#    Copyright (C) 2016  Yves Sagaert & Ivan Svetunkov
# Names of variables ivan41 and yves14 are given in order not to mess with the possible inner loops of "for(i in 1:n)" type.
  l.value <- length(value);
  if(!is.null(value)){
    for(ivan41 in 1:l.value){
      if(substring(value[ivan41],1,1)!="$"){
        value[ivan41] <- paste0("$",value[ivan41])
      }
    }
  }
  else{
    value <- "";
    l.value <- 1;
  }
  
# Check if the data is vector or ts
  if(!is.numeric(data) & !is.ts(data)){
    stop("The provided data is not a vector or ts object! Can't work with it!", call.=FALSE);
  }
  
  if(parallel==TRUE){
    if(!requireNamespace("foreach", quietly = TRUE)){
      stop("In order to run the function in parallel, 'foreach' package must be installed.", call. = FALSE)
    }
    if(!requireNamespace("parallel", quietly = TRUE)){
      stop("In order to run the function in parallel, 'parallel' package must be installed.", call. = FALSE)
    }
# Detect number of cores for parallel calculations
    crs <- min(parallel::detectCores() - 1, origins);

# Check the system and choose the package to use
    if(Sys.info()['sysname']=="Linux"){
      if(requireNamespace("doMC", quietly = TRUE)){
        doMC::registerDoMC(crs);
        cl <- NULL;
      }
      else if(requireNamespace("doParallel", quietly = TRUE)){
        cat(paste0("Setting up ", crs, " clusters using 'doParallel'..."));
        cat("\n");
        cl <- parallel::makeCluster(crs);
        doParallel::registerDoParallel(cl);
      }
      else{
        stop("Sorry, but in order to run the function in parallel, you need either 'doMC' (prefered) or 'doParallel' packages.",
             call. = FALSE)
      }
    }
    else{
# if(Sys.info()['sysname']=="Windows")
      if(requireNamespace("doParallel", quietly = TRUE)){
        cat(paste0("Setting up ", crs, " clusters using 'doParallel'..."));
        cat("\n");
        cl <- parallel::makeCluster(crs);
        doParallel::registerDoParallel(cl);
      }
      else{
        stop("Sorry, but in order to run the function in parallel, you need 'doParallel' package.",
             call. = FALSE)
      }
    }
  }

# Write down the start and frequency
  data.start <- start(data)
  data.freq <- frequency(data)
  
# Write down the entire dataset y and the original horizon
  y <- data
  hh <- h
  obs <- length(y)
  in.sample <- obs - origins

  actuals <- matrix(NA,nrow=h,ncol=origins)
  colnames(actuals) <- paste0("origin",c(1:origins))
  rownames(actuals) <- paste0("h",c(1:h))

  forecasts <- list(NA)
  if(silent==FALSE & parallel==FALSE){
    cat(paste0("Origins done:  "))
  }
  else if(silent==FALSE & parallel==TRUE){
    cat(paste0("Working..."))
  }

##### Start the main loop #####
  if(parallel==FALSE){
    if(co==FALSE){
      for(ivan41 in 1:origins){
# Adjust forecasting horizon to not exeed the sample size
        h <- min(hh,obs - (in.sample+ivan41-1))
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq)
        }
        else{
          data <- ts(y[ivan41:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq)
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
        for(yves14 in 1:l.value){
          forecasts[[(ivan41-1)*l.value+yves14]] <- eval(parse(text=paste0("o.m",value[yves14])))
        }
        actuals[1:h,ivan41] <- y[(in.sample+ivan41):(in.sample+ivan41+h-1)]
        if(silent==FALSE){
          cat(paste(rep("\b",nchar(ivan41)),collapse=""))
          cat(ivan41)
        }
      }
    }
    else{
      for (ivan41 in 1:origins){
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq)
        }
        else{
          data <- ts(y[ivan41:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq)
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
        for(yves14 in 1:l.value){
          forecasts[[(ivan41-1)*l.value+yves14]] <- eval(parse(text=paste0("o.m",value[yves14])))
        }
        actuals[,ivan41] <- y[(in.sample+ivan41-h+1):(in.sample+ivan41)]
        if(silent==FALSE){
          cat(paste(rep("\b",nchar(ivan41)),collapse=""))
          cat(ivan41)
        }
      }
    }
  }
  else{
##### Use foreach for the loop #####
# But first make the list of the needed packages to pass to doParallel
    callenvir <- globalenv();
    callpackages <- search();
    callpackages <- callpackages[c(-1,-length(callpackages))];
    callpackages <- callpackages[substring(callpackages,1,7)=="package"];
    callpackages <- substring(callpackages,9,nchar(callpackages));
    callpackages <- callpackages[callpackages!="timeDate"];
    callpackages <- callpackages[callpackages!="zoo"];
    callpackages <- callpackages[callpackages!="stats"];
    callpackages <- callpackages[callpackages!="graphics"];
    callpackages <- callpackages[callpackages!="grDevices"];
    callpackages <- callpackages[callpackages!="utils"];
    callpackages <- callpackages[callpackages!="datasets"];
    callpackages <- callpackages[callpackages!="methods"];

    if(co==FALSE){
      forecasts <- foreach::`%dopar%`(foreach::foreach(ivan41=1:origins, .packages=callpackages, .export=ls(envir=callenvir)),{
# Adjust forecasting horizon to not exeed the sample size
        h <- min(hh,obs - (in.sample+ivan41-1));
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq);
        }
        else{
          data <- ts(y[ivan41:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq);
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call));
# Save the forecast and the corresponding actuals in matrices
        for(yves14 in 1:l.value){
          forecasts[[yves14]] <- eval(parse(text=paste0("o.m",value[yves14])));
        }

        return(forecasts);
      })
    }
    else{
      forecasts <- foreach::`%dopar%`(foreach::foreach(ivan41=1:origins, .packages=callpackages, .export=ls(envir=callenvir)),{
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq)
        }
        else{
          data <- ts(y[ivan41:(in.sample-h+ivan41-1)],start=data.start,frequency=data.freq)
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
        for(yves14 in 1:l.value){
          forecasts[[yves14]] <- eval(parse(text=paste0("o.m",value[yves14])))
        }

        return(forecasts);
      })
    }

# Make the needed list if there were several values
    forecasts <- unlist(forecasts,recursive=FALSE);

# Check if the clusters have been made
    if(!is.null(cl)){
      parallel::stopCluster(cl);
    }

# Form matrix of actuals in a different loop...
    if(co==FALSE){
      for (ivan41 in 1:origins){
        actuals[1:h,ivan41] <- y[(in.sample+ivan41):(in.sample+ivan41+h-1)]
      }
    }
    else{
      for (ivan41 in 1:origins){
        actuals[,ivan41] <- y[(in.sample+ivan41-h+1):(in.sample+ivan41)]
      }
    }
  }

  if(silent==FALSE){
    cat("\n")
  }

  returned.list <- list(actuals)

# Create several matrices from the list
  if(all(value=="")){
    value <- names(forecasts[[1]]);
    if(length(names(forecasts[[1]][[1]]))>1){
      value.start <- 2;
    }
    else{
      value.start <- 1;
    }
    l.value <- length(value);
    forecasts <- as.list(unlist(forecasts,recursive=FALSE));
  }
  else{
    value <- substring(value,2,nchar(value));
    value.start <- 1;
  }

  for(ivan41 in value.start:l.value){
    stuff.max.length <- max(length(forecasts[[ivan41]]),length(forecasts[[l.value*(origins-1)+ivan41]]))
    stuff <- matrix(NA,nrow=stuff.max.length,ncol=origins)
    colnames(stuff) <- colnames(actuals)
    for(yves14 in 1:origins){
      stuff.length <- length(forecasts[[(yves14-1)*l.value+ivan41]])
      stuff[1:stuff.length,yves14] <- forecasts[[(yves14-1)*l.value+ivan41]]
    }
    returned.list[[ivan41+1]] <- stuff
  }

  names(returned.list) <- c("actuals",value)
  return(returned.list)
}