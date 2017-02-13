plot.net <- function(x, r=1, ...){
# Plot MLP or ELM object

    neuron.col <- "lightgrey"
    xreg.col <- "lightblue"
    season.col <- "lightpink"

    ttl <- class(x)
    if (!any(sapply(ttl,function(x){x == c("elm","mlp","elm.fast")}))){
        stop("Model must be the output of either mlp or elm functions.")
    }

    # If elm.fast change title to ELM
    is.elm.fast <- any(ttl=="elm.fast")
    if (is.elm.fast && any(ttl == "elm")){
        ttl <- setdiff(ttl,"elm.fast")
    }
        
    # Check requested repetition
    if (is.elm.fast){
        reps <- length(x$b)
    } else {
        reps <- length(x$net$weights)    
    }
    if (r>reps){
        stop(paste0("Requested training repetition ",r," with only ", reps, " available."))
    }
    
    # Get network information
    if (is.elm.fast){
        layer.n <- 1 + 1 # +1 for output layer
        varnames <- rownames(x$W.in[[r]])[2:dim(x$W.in[[r]])[1]] # First is Bias
        if (is.null(varnames)){
            varnames <- paste0("X",1:(dim(x$W.in[[r]])[1]-1)) # -1 for Bias
        }
    } else {
        net <- x$net
        layer.n <- length(net$weights[[r]])
        varnames <- net$model.list$variables
    }
    layer.size <- vector("numeric",layer.n+1)
    layer.xx <- c(0,seq(0.1,0.9,length.out=layer.n+1),1)
    layer.yy <- vector("list",layer.n+1)
    layer.size[1] <- length(varnames)
    inputs.col <- rep(neuron.col,layer.size[1])
    inputs.col[grepl("Xreg.",varnames,fixed=TRUE)] <- xreg.col
    inputs.col[grepl("D",varnames,fixed=TRUE)] <- season.col
    layer.yy[[1]] <- seq(0,1,length.out=layer.size[1]+2)
    for (i in 1:layer.n){
        if (is.elm.fast){
            if (i == 1){
                layer.size[i+1] <- x$hd[r] # Single hidden layer
            } else {
                layer.size[i+1] <- 1 # Output
            }
        } else {
            layer.size[i+1] <- dim(net$weights[[r]][[i]])[2]
        }
        layer.yy[[i+1]] <- seq(0,1,length.out=layer.size[i+1]+2)
    }
    # Size of neurons
    rd <- max(0.015,1/((max(layer.size)+2)*1.75))
    rd <- min(rd,0.06)

    # Start plotting
    plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",bty="n",main=toupper(ttl))
    # Draw weights
    for (k in 1:(layer.n-as.numeric(ttl=="elm"))){
        if (layer.size[k] > 0){
            for (i in 1:layer.size[k]){
                if (layer.size[k+1] > 0){
                    for (j in 1:layer.size[k+1]){
                        lines(c(layer.xx[k+1]+rd,layer.xx[k+2]-rd),c(layer.yy[[k]][i+1],layer.yy[[k+1]][j+1]))
                    }
                }
            }
        }
    }

    # Draw output layer for ELM
    cmp <- rep("black",layer.size[layer.n])
    if (ttl == "elm"){
        ltt <- rep(1,layer.size[layer.n])
        # For neuralnets based ELM grey-out connections of unused neurons
        if (!is.elm.fast){
            w <- x$W[[r]][2:(layer.size[layer.n]+1)] != 0
            cmp[!w] <- "grey"
            ltt[!w] <- 2
        }
        if (layer.size[layer.n] > 0){
            for (i in 1:layer.size[layer.n]){
                lines(c(layer.xx[layer.n+1]+rd,layer.xx[layer.n+2]-rd),c(layer.yy[[layer.n]][i+1],layer.yy[[layer.n+1]][2]),col=cmp[i],lty=ltt[i])
            }
        }
    }
    
    # Draw neurons
    for (k in 1:(layer.n+1)){
        if (layer.size[k]>0){
            for (i in 1:layer.size[k]){
                if (ttl == "elm" & k == (layer.n)){
                    draw.circle(layer.xx[k+1],layer.yy[[k]][i+1],rd,col=neuron.col,border=cmp[i])    
                } else {
                    if (k == 1){
                      draw.circle(layer.xx[k+1],layer.yy[[k]][i+1],rd,col=inputs.col[i])
                    } else {
                      draw.circle(layer.xx[k+1],layer.yy[[k]][i+1],rd,col=neuron.col)
                    }
                }
                
            }
        }
    }
    
    # Draw inputs
    for (i in 1:layer.size[1]){
        arrows(0,layer.yy[[1]][i+1],layer.xx[2]-rd,length=0.1,code=2)
    }
    
    # Draw outputs
    for (i in 1:layer.size[layer.n+1]){
        arrows(layer.xx[layer.n+2]+rd,layer.yy[[layer.n+1]][i+1],1,length=0.1,code=2)
    }
    
    # Draw direct
    if (ttl == "elm"){
        if (x$direct == TRUE){
            wd <- x$W.dct[[r]]
            wd.n <- length(wd)
            if (wd.n > 0){
                for (i in 1:wd.n){
                    if (wd[i] != 0){
                        lines(c(layer.xx[2]+rd,layer.xx[layer.n+2]-rd),c(layer.yy[[1]][i+1],layer.yy[[layer.n+1]][2]),col="blue",lty=1)   
                    }
                }
            }
        }
    }
        
    # Add x-axis
    if (layer.size[1]>1){inp<-paste0("Inputs\n(",layer.size[1],")")}else{inp<-"Input"}
    if (layer.n>2){lay<-paste0("Hidden ",1:(layer.n-1))}else{lay<-"Hidden"}
    for (i in 1:(layer.n-1)){
        lay[i] <- paste0(lay[i],"\n(",layer.size[1+i],")")
    }
    axis(1,at=layer.xx[2:(layer.n+2)],labels=c(inp,lay,"Output"),col=NA)

}