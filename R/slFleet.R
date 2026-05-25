
#' Make a Short-Lived Species Fleet Object
#'
#' A function that builds the fleet component of an operating model for a short-lived species. Requires from
#' specifications for growth, stock-recruitment, natural mortality rate and spatial distribution. This
#' is a wrapper function and produces a fully customizable MSEtool object of class 'stock'. You can edit
#' that object to include very complex stock dynamics if desired.
#'
#' @param Name Character string. The name of the stock object. E.g., "North Pacific Neon Flying Squid"
#' @param Species Character string. The genus and species. E.g., "Ommastrephes bartramii"
#' @param CommonName Character string. The common name. E.g., "Neon Flying Squid"
#' @param nYear Positive Integer. The number of historical years. E.g., 10, (2015 - 2024)
#' @param pYear  Positive Integer. The number of projection years. E.g., 10 (2025 - 2034)
#' @param Seasons Positive Integer. The number of sub-year time steps. E.g., 12 for a monthly model
#' @param CurrentYear Integer. The Calendar year of the last historical year. E.g. 2024.
#' @param nSim Positive Integer. The number of independent simulations. E.g., 4 for testing, 24 for preliminary results, 48 for representative results, 192+ for informing management.
#' @param rec_age Positive Integer. The age the fish is recruited (in seasonal time steps). E.g., given a monthly model (Seasons = 12) a value of 4 means that the stock recruitment function calculates recruiment into the time step 4 months later.
#' @param nAges Positive Integer. The number of seasonal time steps that population dynamics will be calculated for. E.g., 24 would be two years in a monthly model. Note that this only necessary if creatures live that long or dynamics are not suitably approximated with a plus group.
#' @param Effort A matrix of positive real numbers nSim x time steps (nYear x Seasons). Given a default q value of 1 this is the apical fishing mortality rate. Defaul is 'NA' and in this case effort pattern is simulated with a seasonal and temporal trend.
#' @param sel50 Positive real number. The age (in seasons) that 50% of individuals are selected, in a logistic selectivity model. E.g., 6 in a monthly model would be 50% selected after 6 months.
#' @param selSLP Positive real number. The slope of the logistic selectivity model. E.g., 2 in a monthly model.
#' @param plot Boolean. Should plots be presented on the creation of the fleet object.
#' @return An object of MSEtool class fleet
#' @author T. Carruthers
#' @examples
#' A_short_lived_fleet<- slFleet("Default fleet values for demo")
#' class(A_short_lived_fleet)
#' @seealso \link{slFleet} for making a short-lived stock object and \link{slOM} for specifying the entire operating model from fleet and stock objects.
#' @export
slFleet = function(Name = "A fleet", nYear = 10, pYear = 10, Seasons = 12, CurrentYear = 2026,
                   nSim = 4, rec_age = 1, nAges = 24, Effort = NA, sel50 = 6, selSLP = 2,
                   plot = F){

  # Default args for testing:
  # Name = "A fleet"; nYear = 10; pYear = 10; Seasons = 12; nArea = CurrentYear = 2026; nSim = 4; rec_age = 4; nages = 24; Effort = NA;
  # sel50 = 6; selSLP = 2

  fleet = Fleet()
  fleet@Name = Name
  fleet@nYear = nYear
  fleet@pYear = pYear
  fleet@CurrentYear = CurrentYear
  fleet@Years = CurrentYear+ (-(nYear-1):pYear)
  fleet@Seasons = Seasons
  fleet@nSim = nSim

  if(class(Effort)!="array"){
    Effort = Effort_sim(nSim, nYear, Seasons, nArea, ymin = 0.25, yfac = 0.5, ECV = 0.15, maxF = 0.1, plot=plot)
  }

  fleet@Effort = Effort(Effort=Effort)

  fleet@Catchability = Catchability(Efficiency = 1)

  svec = ((rec_age : nAges)-sel50)*selSLP
  fleet@Selectivity = Selectivity(MeanAtAge = exp(svec)/(1+exp(svec)))

  # Optional slots:
  # Distribution - only needed if you don't want openMSE to solve the distribution according to vulnerable biomass

  fleet

}

SL_fleet_check = function(fleet){
  # fleet spec checks

}



# === Fleet internal ==========================================================================================

# Invent Effort

Effort_sim = function(nSim, nYear, Seasons, nArea, ymin = 0.25, yfac = 0.5, ECV = 0.15, maxF = 0.1, plot=F){
  nt =  nYear*Seasons
  Effort = array(0,c(nSim,nt))
  ye = 1 + ymin + (sin(((9+1:(nYear))/3))*1)
  se = dnorm(1:Seasons,Seasons/2,Seasons/5)
  seind = t(array(1:Seasons,c(nt,nSim)))
  yind = t(array(rep(1:nYear,each=Seasons),c(nt,nSim)))
  Effort[] = ye[yind] * se[seind] * trlnorm(nt*nSim,1,ECV)
  Effort = Effort / apply(Effort,1,max) * maxF
  if(plot)matplot(t(Effort),type="l",lty=1)
  Effort
}



