
# demo 2 - fleet high contrast setup

#' @export
demo_fleets = function(nyear, pyear, Seasons, nSim, nAges, stock){
  Short_lowF = slFleet(Name = "Short_lowF", sel50 = 1, selSLP = 4, relE_min = 0.3, relE_max = 0.75, maxF = 1.3, nYear = nYear, pYear = pYear, Seasons = Seasons, nSim = nSim, nAges = nAges, plot=T)
  Long_highF = slFleet(Name = "Long_highF", sel50 = 3.0, selSLP = 5, relE_loc = 0.1, relE_freq = 1, maxF = 0.65, smu = 0.7, scv = 0.15, nYear = nYear, pYear = pYear, Seasons = Seasons, nSim = nSim, nAges = nAges, plot=T)
  fleet = slCombineFleets(fleetlist = list( Short_lowF = Short_lowF, Long_highF = Long_highF), stock)
  fleet
}

# demo 3 - stock om

#' @export
demo_stocks = function(nyear, pyear, Seasons, nSim, nAges, CurrentYear){

  small = slStock(Name = "Small Pacific JFS", Species = "Dosidicus gigas", CommonName = "Jumbo Flying Squid",
                  nYear = nYear, pYear = pYear, Seasons = Seasons,              # Historical years, projection years, subyears
                  CurrentYear = CurrentYear, nSim = nSim,                       # Last historical year, number of simulations
                  rec_age = 1, nages = nAges, PlusGroup = F,                    # Season / subyear to recruitment, max seasons, senescence after nages
                  spawndist = c(0, 0.4, 0.5, 0.1),                              # Spawning distribution among sub-years
                  Len_age =  c(25.3, 31.1, 34.6, 36.7, 38.0, 38.8, 39.3, 39.6), # Length at age (cm) # round(40 * (1-exp(-0.5*2:9)),1)
                  Len_CV = 0.2, lenclasses = seq(10,100,by=2),
                  Wt_age = c(1.6, 3.0, 4.1, 4.9, 5.5, 5.8, 6.1, 6.2),           # Weight at age (kg)  # round(1E-4* (40 * (1-exp(-0.5*2:9)))^3,1)
                  M = 0.45,                                                     # Instantaneous natural mortality rate (per season)
                  Mat_age = c(0.18, 0.60, 0.88, 0.99, 1.00, 1.00, 1.00, 1.00),  # Spawning fraction at age ('maturity at age') # avec = ((1 : 8)-2) * 1.8; round(exp(avec)/(1+exp(avec)),2)
                  SR_type = "Ricker", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                  nareas = nAreas,
                  Frac_area = matrix(c(0.90, 0.80, 0.75, 0.60, 0.45, 0.30, 0.25, 0.20,   # Area 1
                                       0.05, 0.15, 0.20, 0.25, 0.30, 0.40, 0.40, 0.40,   # Area 2
                                       0.05, 0.05, 0.05, 0.15, 0.25, 0.30, 0.35, 0.40),  # Area 3
                                     byrow = T, nrow=nAreas),
                  prob_stay = 0.8)                                                       # Tendency to remain in area



  medium = slStock(Name = "Medium Pacific JFS", Species = "Dosidicus gigas", CommonName = "Jumbo Flying Squid",
                   nYear = nYear, pYear = pYear, Seasons = Seasons,              # Historical years, projection years, subyears
                   CurrentYear = CurrentYear, nSim = nSim,                       # Last historical year, number of simulations
                   rec_age = 1, nages = nAges, PlusGroup = F,                    # Season / subyear to recruitment, max seasons, senescence after nages
                   spawndist = c(0, 0.2, 0.6, 0.2),                              # Spawning distribution among sub-years
                   Len_age =  c(41.9, 50.1, 54.6, 57.0, 58.4, 59.1, 59.5, 59.7), # Length at age (cm)  # round(60 * (1-exp(-0.6*2:9)),1)
                   Len_CV = 0.2, lenclasses = seq(10,100,by=2),
                   Wt_age = c(7.4, 12.6, 16.2, 18.5, 19.9, 20.6, 21.1, 21.3),    # Weight at age (kg) # round(1E-4* (60 * (1-exp(-0.6*2:9)))^3,1)
                   M = 0.4,                                                     # Instantaneous natural mortality rate (per season)
                   Mat_age = c(0.14, 0.50, 0.86, 0.97, 1.00, 1.00, 1.00, 1.00),  # Spawning fraction at age ('maturity at age') # avec = ((1 : 8)-2) * 1.8; round(exp(avec)/(1+exp(avec)),2)
                   SR_type = "Ricker", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                   nareas = nAreas, Frac_area = matrix(c(rep(0.2,nAges),rep(0.5,nAges),rep(0.3,nAges)), byrow = T, nrow=nAreas),
                   prob_stay = 0.8)                                              # Tendency to remain in area



  large = slStock(Name = "Large Pacific JFS", Species = "Dosidicus gigas", CommonName = "Jumbo Flying Squid",
                  nYear = nYear, pYear = pYear, Seasons = Seasons,              # Historical years, projection years, subyears
                  CurrentYear = CurrentYear, nSim = nSim,                       # Last historical year, number of simulations
                  rec_age = 1, nages = nAges, PlusGroup = F,                    # Season / subyear to recruitment, max seasons, senescence after nages
                  spawndist = c(0, 0.4, 0.5, 0.1),                              # Spawning distribution among sub-years
                  Len_age =  c(23.1, 38.5, 48.9, 55.9, 60.5, 63.6, 65.7, 67.1), # Length at age (cm) # round(70 * (1-exp(-0.4*1:8)),1)
                  Len_CV = 0.2, lenclasses = seq(10,100,by=2),
                  Wt_age = c(5.7, 11.7, 17.4, 22.2, 25.8, 28.4, 30.3, 31.6),    # Weight at age (kg) # round(1E-4* (70 * (1-exp(-0.4*2:9)))^3,1)
                  M = 0.35,                                                      # Instantaneous natural mortality rate (per season)
                  Mat_age = c(0.03, 0.14, 0.50, 0.86, 0.97, 1.0, 1.0, 1.0, 1.0),# Spawning fraction at age ('maturity at age')# avec = ((1 : 9)-3) * 1.8; round(exp(avec)/(1+exp(avec)),2)
                  SR_type = "Ricker", h = 0.9, sigmaR = 1.0, trunc_sigmaR = 2.0, R_AC = 0.5, R0 = 1E6,
                  nareas = nAreas,
                  Frac_area = matrix(c(0.05, 0.05, 0.05, 0.15, 0.25, 0.30, 0.35, 0.40,           # Area 1
                                       0.05, 0.15, 0.20, 0.25, 0.30, 0.40, 0.40, 0.40,           # Area 2
                                       0.90, 0.80, 0.75, 0.60, 0.45, 0.30, 0.25, 0.20),          # Area 3
                                     byrow = T, nrow=nAreas),
                  prob_stay = 0.8)                                                               # Tendency to remain in area


  stock = list(Small_Pacific_JFS = small, Medium_Pacific_JFS = medium, Large_Pacific_JFS = large)
  stock

}
