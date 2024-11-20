# Base models functions
extract_bm <- function(x, type = "base"){
  type <- match.arg(type, c("base", "fc", "fit", "model", "res"))
  if(!("base_model" %in% class(x))){
    stop("x has to be a base_model obj.")
  }else{
    x[[type]]
  }
}

base_model <- function(model, x, h){
  class(model) <- model
  UseMethod("base_model", model)
}

base_model.snaive <- function(model = "snaive", x, h){
  fc <- forecast::snaive(ts(x, frequency = h), h)
  base <- unname(fc$mean)
  res <- unname(as.vector(as.numeric(x) - fitted(fc)))
  return(structure(list(model = "snaive",
                        fit = fc,
                        fc = fc,
                        base = base,
                        res = res), class = "base_model"))
}

base_model.arima <- function(model = "arima", x, h){
  fit <- forecast::auto.arima(x)
  fc <- forecast::forecast(fit, h = h)
  base <- unname(fc$mean)
  res <- unname(as.vector(as.numeric(x) - fitted(fit)))
  return(structure(list(model = "arima",
                        fit = fit,
                        fc = fc,
                        base = base,
                        res = res), class = "base_model"))
}

base_model.tbats <- function(model = "tbats", x, h){
  fit <- forecast::tbats(x, biasadj = TRUE)
  fc <- forecast::forecast(fit, h = h)
  base <- unname(fc$mean)
  res <- unname(as.vector(as.numeric(x) - fitted(fit)))
  return(structure(list(model = "tbats",
                        fit = fit,
                        fc = fc,
                        base = base,
                        res = res), class = "base_model"))
}

base_model.stlf <- function(model = "stlf", x, h){
  if(frequency(x) <= 1){
    freq <- 7
  }else{
    freq <- frequency(x)
  }
  x_ts <- ts(x, frequency = freq)
  
  fit <- fc <- forecast::stlf(x_ts, method = "ets", h = h)
  
  base <- unname(fc$mean)
  res <- unname(as.vector(as.numeric(x) - fitted(fit)))
  return(structure(list(model = "stlf",
                        fit = fit,
                        fc = fc,
                        base = base,
                        res = res), class = "base_model"))
}
