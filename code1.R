#generated data
data = rgamma(365,5,5)
p99 = quantile(data, .99) 
p95 = quantile(data, .95)
#functions to get output of 10 variables

RX5day = function(data){
  return(max(sapply(1:361, function(i) mean(data[i:(i+4)]))))
}

R99p = function(data, p99){
  return(sum(data[data>p99]))
}

Rx1day = function(data){return(max(data))}

R95p = function(data, p95){
  return(sum(data[data>p95]))
}

R95pT = function(data, p95){
  return(sum(data>p95))
}

SDII = function(data){
  return(sum(data)/sum(data>1))
} 

CWD = function(data){
  run_length = rle(data>1)
  return(max(run_length$lengths[run_length$values == 1]))
}

R10mm = function(data){
  return(sum(data>10))
}

PRCPTOT = function(data){
  return(sum(data[data>1]))
}

R20mm = function(data){
  return(sum(data>20))
}



precip_vars = function(data, p99, p95){
  Output = list(RX5day = RX5day(data),
                R99p = R99p(data, p99),
                Rx1day = Rx1day(data),
                R95p = R95p(data, p95),
                R95pT = R95pT(data, p95),
                SDII = SDII(data),
                CWD = CWD(data),
                R10mm = R10mm(data),
                PRCPTOT = PRCPTOT(data),
                R20mm = R20mm(data))
  
  return(Output)
}

precip_vars(data, p99, p95)
