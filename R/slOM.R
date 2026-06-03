
#' Construct Complete Short-Lived Species Operating Model
#'
#' A function that combines stock and fleet object in a complete operating model for a short-lived species. Can
#' include lists of fleets and stocks for more complex operating models
#'
#' @param Name Character string. The name of the operating model. E.g., "NP NFS Ref OM #3: high M, high h, low rec"
#' @param Agency Character string. The relevant fishery agency "North Pacific Fishery Council (NPFC)"
#' @param Author Character string. Name of the author(s) of the operating model E.g., "A. Person"
#' @param Email Character string. Email address of the corresponding author
#' @param Region Character string. The location of the stock / management area.
#' @param ComplexName Character string. The name of the multi-stock complex (if more than 1 stock is being simulated).
#' @param Latitude Real number. Degrees N. Negative is southern hemisphere
#' @param Longitude Real number. Degrees E. Negative is western hemisphere
#' @param Sponsor Character string. The name of the person(s) or organization(s) that paid for the research.
#' @param nSim Positive Integer. The number of independent simulations. E.g., 4 for testing, 24 for preliminary results, 48 for representative results, 192+ for informing management.
#' @param nYear Positive Integer. The number of historical years. E.g., 10, (2015 - 2024)
#' @param pYear  Positive Integer. The number of projection years. E.g., 10 (2025 - 2034)
#' @param Seasons Positive Integer. The number of sub-year time steps. E.g., 12 for a monthly model
#' @param CurrentYear Integer. The Calendar year of the last historical year. E.g. 2024.
#' @param Interval Positive integer. The management interval (in seasons) - how frequently is new advice provided. In a monthly model, a value of 6 would mean every 6 months after June and December.
#' @param Seed Real number. The seed for sampling of random variables.
#' @param stock A stock object or list of stock objects
#' @param fleet A fleet object or list of fleet objects
#' @return An object of MSEtool class om
#' @author T. Carruthers
#' @examples
#' A_short_lived_om <- slOM("Default fleet values for demo")
#' class(A_short_lived_om)
#' @seealso \link{slStock} for making a short-lived stock object and \link{slFleet} for making a short-lived fleet object.
#' @export
slOM = function(Name = "Short-lived simulation", Agency = "A fishery agency", Author = "A fishery analyst",
                 Email = "a.person@email.com", Region = "A fishery management area", ComplexName = "Short-Lived Complex", Latitude = NA, Longitude = NA,
                 Sponsor = "A generous funder", nSim = 4, nYear = 10, pYear = 10, Seasons = 12, CurrentYear = 2026,
                 Interval = 6, Seed = 1, stock = NA, fleet = NA, obs = NA){

  #  Name = "Short-lived simulation"; Agency = "A fishery agency"; Author = "A fishery analyst"; Email = "a.person@email.com"; Region = "A fishery management area"; ComplexName = "Short-Lived Complex"; Latitude = NA; Longitude = NA
  #  Sponsor = "A generous funder"; nSim = 4; nYear = 10; pYear = 10; Seasons = 12; CurrentYear = 2026; Interval = 6; Seed = 1; obs = NA; stock = NA; fleet = NA

 # om = OM()     # Operating model object
  om = new('om')
  # Dimensions and labels
  om@Name = Name
  om@Agency = Agency
  om@Author = Author
  om@Email = "An email address"
  om@Region = "A region"
  om@Latitude = 0
  om@Longitude = 0
  om@Sponsor = "A sponsor"
  om@nSim = nSim
  om@nYear = nYear
  om@pYear = pYear
  om@CurrentYear = CurrentYear
  om@Seasons = Seasons
  om@Interval = Interval
  om@Seed = Seed

  #data = new('data')
  #data@YearLH = CurrentYear
  #om@Data=data
  # imp

  if(!(class(obs)=="obs")) obs = slObs()
  if(!(class(stock)%in%c("list","stock"))) stock = slStock()
  if(!(class(fleet)=="StockFleetList")) fleet = slFleet()

  nstocks = length(stock)
  om@Stock = stock
  om@Fleet = fleet
  om@Obs = obs # you need to add biomass selectivity to survey - also catch at length is currently done post hoc

  om = Populate(om)

  if(nstocks > 1){                   # Complex specification
    om@Complexes = list(1:nstocks)
    names(om@Complexes) = ComplexName
  }

  om

}

slOM_check=function(){
  # om spec checks

}

