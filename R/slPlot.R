slplot = function(obj){

  if(class(obj)=="hist") slplot.hist(obj)
  if(class(obj)=="slSimData") slplot.data(obj)
  if(class(obj)=="slSimSam") slplot.simsam(obj)

}

