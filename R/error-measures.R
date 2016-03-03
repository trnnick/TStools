MPE <- function(actual,forecast,round=3)
{
# This function calculates Mean / Median Percentage Error
# actual - actual values,
# forecast - forecasted or fitted values.
    result <- round(mean((actual-forecast)/actual,na.rm=TRUE),digits=round);
    return(result);
}

MAPE <- function(actual,forecast,round=3){
# This function calculates Mean Absolute Percentage Error
# actual - actual values,
# forecast - forecasted values.
    result <- round(mean(abs((actual-forecast)/actual),na.rm=TRUE),round);
    return(result);
}

SMAPE <- function(actual,forecast,round=3)
{
# This function calculates Symmetric Mean / Median Absolute Percentage Error with
# sum of absolute values in the denominator
# actual - actual values,
# forecast - forecasted or fitted values.
    result <- round(mean(2*abs(actual-forecast)/(abs(actual)+abs(forecast)),na.rm=TRUE),digits=round);
    return(result);
}

MASE <- function(actual,forecast,scale,round=3){
# This function calculates Mean Absolute Scaled Error as in Hyndman & Koehler, 2006
# actual - actual values,
# forecast - forecasted values.
# scale - the measure to scale errors with. Usually - MAE of in-sample.
    result <- round(mean(abs(actual-forecast),na.rm=TRUE)/scale,round)
    return(result)
}

GMRAE <-function(actual,forecast,etalon,round=3){
# This function calculates Geometric Mean Relative Absolute Error
# actual - actual values,
# forecast - forecasted or fitted values.
# etalon - forecasted or fitted values of etalon method.
    result <- round(exp(mean(log(abs(actual-forecast)/abs(actual-etalon)),na.rm=TRUE)),round);
    return(result);
}

hm <- function(x,C=mean(x),digits=5,...)
{
# This function calculates half moment

    x <- x[!is.na(x)];
    result <- round(mean(sqrt(as.complex(x-C)),...),digits=digits);
    return(result);
}