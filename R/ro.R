ro <- function(data,h=1,origins=1,call,value=NULL,
               ci=FALSE,co=FALSE,silent=FALSE,parallel=FALSE){
# Function makes Rolling Origin for the data using the call
#    Copyright (C) 2015  Yves Sagaert & Ivan Svetunkov
  l.value <- length(value);
  if(!is.null(value)){
    for(i in 1:l.value){
      if(substring(value[i],1,1)!="$"){
        value[i] <- paste0("$",value[i])
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
      for(i in 1:origins){
# Adjust forecasting horizon to not exeed the sample size
        h <- min(hh,obs - (in.sample+i-1))
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
        }
        else{
          data <- ts(y[i:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
        for(j in 1:l.value){
          forecasts[[(i-1)*l.value+j]] <- eval(parse(text=paste0("o.m",value[j])))
        }
        actuals[1:h,i] <- y[(in.sample+i):(in.sample+i+h-1)]
        if(silent==FALSE){
          cat(paste(rep("\b",nchar(i)),collapse=""))
          cat(i)
        }
      }
    }
    else{
      for (i in 1:origins){
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
        }
        else{
          data <- ts(y[i:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
        for(j in 1:l.value){
          forecasts[[(i-1)*l.value+j]] <- eval(parse(text=paste0("o.m",value[j])))
        }
        actuals[,i] <- y[(in.sample+i-h+1):(in.sample+i)]
        if(silent==FALSE){
          cat(paste(rep("\b",nchar(i)),collapse=""))
          cat(i)
        }
      }
    }
  }
  else{
##### Use foreach for the loop #####
# But first find out the environment of the "call" function
    callenvir <- environment(get(substring(call,1,which(strsplit(call, "")[[1]]=="(")[1]-1)));
#    callenvir <- globalenv();

    if(co==FALSE){
      forecasts <- foreach::`%dopar%`(foreach::foreach(i=1:origins, .export=ls(envir=callenvir)),{
# Adjust forecasting horizon to not exeed the sample size
        h <- min(hh,obs - (in.sample+i-1));
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+i-1)],start=data.start,frequency=data.freq);
        }
        else{
          data <- ts(y[i:(in.sample-h+i-1)],start=data.start,frequency=data.freq);
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call));
# Save the forecast and the corresponding actuals in matrices
        for(j in 1:l.value){
          forecasts[[j]] <- eval(parse(text=paste0("o.m",value[j])));
        }

        return(forecasts);
      })
    }
    else{
      forecasts <- foreach::`%dopar%`(foreach::foreach(i=1:origins, .export=ls(envir=callenvir)),{
# Make the in-sample
        if(ci==FALSE){
          data <- ts(y[1:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
        }
        else{
          data <- ts(y[i:(in.sample-h+i-1)],start=data.start,frequency=data.freq)
        }
# Evaluate the call string and save to object o.m.
        o.m <- eval(parse(text=call))
# Save the forecast and the corresponding actuals in matrices
        for(j in 1:l.value){
          forecasts[[j]] <- eval(parse(text=paste0("o.m",value[j])))
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
      for (i in 1:origins){
        actuals[1:h,i] <- y[(in.sample+i):(in.sample+i+h-1)]
      }
    }
    else{
      for (i in 1:origins){
        actuals[,i] <- y[(in.sample+i-h+1):(in.sample+i)]
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

  for(i in value.start:l.value){
    stuff.max.length <- max(length(forecasts[[i]]),length(forecasts[[l.value*(origins-1)+i]]))
    stuff <- matrix(NA,nrow=stuff.max.length,ncol=origins)
    colnames(stuff) <- colnames(actuals)
    for(j in 1:origins){
      stuff.length <- length(forecasts[[(j-1)*l.value+i]])
      stuff[1:stuff.length,j] <- forecasts[[(j-1)*l.value+i]]
    }
    returned.list[[i+1]] <- stuff
  }

  names(returned.list) <- c("actuals",value)
  return(returned.list)
}