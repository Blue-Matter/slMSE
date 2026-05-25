

# fleetlist = list(Long_highF,Short_lowf)
slCombineFleets = function(fleetlist, stock){
  nstock = length(stock)
  stocknames = names(stock)
  nfleet = length(fleetlist)
  fleetnames = names(fleetlist)
  sfl = list(); class(sfl) = "StockFleetList"
  for(ss in 1:nstock){
    sfl[[stocknames[ss]]] = list()
    class(sfl[[stocknames[ss]]]) = "FleetList"
    for(ff in 1:nfleet){
      sfl[[stocknames[ss]]][[fleetnames[ff]]] = fleetlist[[ff]]
    }
  }
  sfl 
}
