plot.net <- function(fit,r=1){
# Plot MLP or ELM object

    ttl <- class(fit)
    neuron.col <- "lightgrey"
    xreg.col <- "lightblue"
    season.col <- "lightpink"
    
    if (!any(c("elm","mlp")==ttl)){
        stop("Model must be the output of either mlp or elm functions.")
    }
    
    net <- fit$net

    if (r>length(net$weights)){
        stop(paste0("Training repetition ",r," requested, with only ", length(net$weights), " available."))
    }
    
    layer.n <- length(net$weights[[r]])
    layer.size <- vector("numeric",layer.n+1)
    layer.xx <- c(0,seq(0.1,0.9,length.out=layer.n+1),1)
    layer.yy <- vector("list",layer.n+1)
    layer.size[1] <- length(net$model.list$variables)
    inputs.col <- rep(neuron.col,layer.size[1])
    inputs.col[grepl("Xreg.",net$model.list$variables,fixed=TRUE)] <- xreg.col
    inputs.col[grepl("D",net$model.list$variables,fixed=TRUE)] <- season.col
    layer.yy[[1]] <- seq(0,1,length.out=layer.size[1]+2)
    for (i in 1:layer.n){
        layer.size[i+1] <- dim(net$weights[[r]][[i]])[2]
        layer.yy[[i+1]] <- seq(0,1,length.out=layer.size[i+1]+2)
    }
    # Size of neurons
    rd <- max(0.015,1/((max(layer.size)+2)*1.75))
    rd <- min(rd,0.06)
    
    plot(NA,NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",xaxt="n",yaxt="n",bty="n",main=toupper(ttl))
    
    # Draw weights
    for (k in 1:(layer.n-as.numeric(ttl=="elm"))){
        for (i in 1:layer.size[k]){
            for (j in 1:layer.size[k+1]){
                lines(c(layer.xx[k+1]+rd,layer.xx[k+2]-rd),c(layer.yy[[k]][i+1],layer.yy[[k+1]][j+1]))
            }
        }
    }
    
    cmp <- rep("black",layer.size[layer.n])
    if (ttl == "elm"){
        w <- fit$W[[r]][2:(layer.size[layer.n]+1)] != 0
        
        cmp[!w] <- "grey"
        ltt <- rep(1,layer.size[layer.n])
        ltt[!w] <- 2
        for (i in 1:layer.size[layer.n]){
            lines(c(layer.xx[layer.n+1]+rd,layer.xx[layer.n+2]-rd),c(layer.yy[[layer.n]][i+1],layer.yy[[layer.n+1]][2]),col=cmp[i],lty=ltt[i])
        }
    }
    
    # Draw neurons
    for (k in 1:(layer.n+1)){
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
        if (fit$direct == TRUE){
            wd <- fit$W[[r]][(layer.size[layer.n]+2):length(fit$W[[r]])] != 0
            wd.n <- length(wd)
            wd.cmp <- rep("blue",wd.n)
            wd.ltt <- rep(1,wd.n)
            wd.cmp[!wd] <- "lightblue"
            wd.ltt[!wd] <- 2
            for (i in 1:wd.n){
                lines(c(layer.xx[2]+rd,layer.xx[layer.n+2]-rd),c(layer.yy[[1]][i+1],layer.yy[[layer.n+1]][2]),col=wd.cmp[i],lty=wd.ltt[i])
            }
        }
    }
        
    # Add x-axis
    if (layer.size[1]>1){inp<-"Inputs"}else{inp<-"Input"}
    if (layer.n>2){lay<-paste0("Hidden ",1:(layer.n-1))}else{lay<-"Hidden"}
    axis(1,at=layer.xx[2:(layer.n+2)],labels=c(inp,lay,"Output"),col=NA)

}