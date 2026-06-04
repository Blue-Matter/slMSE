
# ==== Fake data MP =============================================================================

fakeMP = function(Data, repyr = 2028, outdir = "C:/GitHub/slMSE/data", filenam = "Example_Data", Eff=0.5){
  ad=Advice(Effort=Eff)
  if(max(Data@Years)>repyr){
    assign("Example_Data", Data)
    do.call(save, list("Example_Data",file=paste0(outdir,"/Example_Data.rda")))
    cat(paste0("Example data file with data in year ",repyr, " outputted to: ",outfile, "\n"))
    stop()
  }
  ad
}
class(fakeMP) = "mp"


# ==== Generic MPs ==============================================================================

Eff = function(Data, rel_eff = 1){
  Advice(Effort = rel_eff)
}
class(Eff) = "mp"



# ==== Index based MPs ==========================================================================

#' A Generic Index Ratio Management Procedure (OpenMSE v2.0)
#'
#' A function for calculating TAC management advice
#'
#' @param Data An object of openMSE class 'Data'
#' @param HCR_ICP 2-position vector. The Index control points - index relative to last historical year. These are a pair of x-axis control points for a hockey-stick control rule. Phrased in terms of the levels of the generic survey index for the last historical year.
#' @param HCR_LCP  2-position vector. Management Lever control points - TAC relative to Landings in last historical year. These are a pair of y-axis control points for a hockey-stick control rule. Phrased in terms of the levels of catch over the last historical year.
#' @param HCR_sens Sensitivity of the response. A value of 1 is proportional.
#' @param HCR_up_max The maximum rate of advice increase. Positive imperfect fraction. E.g., 0.2 is a maximum management increase of +20 per cent.
#' @param HCR_down_max The maximum rate of advice decrease. Positive fraction. E.g. 0.8 is a maximum management decreate of 80 per cent.
#' @param HCR_smooth The effective number of parameters of a polynomial smoother applied to the input indices. Defaults to NA - raw data. A value of 0.1 means that there will be length(Index)*0.1 number of smoothing parameters. Hence, higher values are less smoothing. This parameterization keeps smoothing consistent as the length of the time series increases in projections.
#' @param plot Should plots be produced that explain MP calculations (e.g. HCRs)? Boolean.
#' @examples
#' GIR(Example_Data)
#' @author T. Carruthers
#' @export
GIR = function(Data, HCR_ICP = c(0.5, 1.5), HCR_LCP = c(0.5, 1.5), HCR_sens = 1, HCR_up_max = 0.2, HCR_down_max = 0.2, HCR_smooth = 0.2, plot = F){

  ad = Advice()
  # other things can go in here!
  ad = do_HCR(Data, HCR_ICP, HCR_LCP, HCR_sens, HCR_up_max, HCR_down_max,
              HCR_smooth, plot = plot, onlyHCR = onlyHCR, ad)

  ad

}
class(GIR) = "mp"

# myMSE = Project(hist, MPs = "fakeMP")
#  load("C:/GitHub/slMSE/data/Example_Data.rda"); Data = Example_Data; HCR_ICP = c(0.5, 1.5); HCR_LCP = c(0.5,1.5); HCR_sens = 1; HCR_smooth = 0.2
#  HCR_up_max = 0.2; HCR_down_max = 0.2; plot=T; onlyHCR = F

do_HCR = function(Data, HCR_ICP, HCR_LCP, HCR_sens, HCR_up_max, HCR_down_max, HCR_smooth, plot, onlyHCR, ad){

  # Notes:
  # Relative TAC level (relative to last historical year) goes in in Data@Misc$rel_lev_vals and out in ad@Misc$rel_lev_vals (for determining acceptible rate of change)

  d = get_dims(Data); nt = d$nt; nf = d$nf;
  season = get_season(Data)

  ad = Advice()
  TAC = rep(NA,d$nf)

  old_rel_lev_vals = Data@Misc$rel_lev_vals
  if(is.null(old_rel_lev_vals))old_rel_lev_vals = 1 # if first advice year
  new_rel_lev_vals = NA

  Index = Data@Survey@Value[,1] / mean(Data@Survey@Value[refhistts(Data, 1:4),1]) # calibrated to be 1 over last historical year

  if(season == 1){ # Do control calcs for season 1 - updates are annual so can skip this for seasons 2-4
    if(!is.na(HCR_smooth)) Index = Index_smooth(Index, HCR_smooth, plot = plot, plotname = "HCR Index Smooth #") # polynomial smoother parameterized as 'effective number of parameters' - a fraction of the total length of the time series
    xCP = HCR_ICP
    yCP = HCR_LCP
    curx = mean(Index[refrecentts(Data)])                      # last historical annual index lagged by specified year
    rel_lev_val = old_rel_lev_vals
    new_rel_lev_vals = do_HockStick(curx, xCP, yCP, reflev=1, rel_lev_val, HCR_up_max_hh = HCR_up_max, HCR_down_max_hh = HCR_down_max, HCR_sens_hh = HCR_sens, plot=plot, onlyHCR = onlyHCR, levlab="TAC (relative to catch in last historical year)", Data=Data)
  }else{           # same as season 1 adjustment
    new_rel_lev_vals = old_rel_lev_vals
  }

  ad@TAC = Data@Landings@Value[refhistts(Data, Season),] * new_rel_lev_vals
  ad@Misc$rel_lev_vals = new_rel_lev_vals
  ad

}

do_HockStick = function(curx, xCP, yCP, reflev, rel_lev_val, HCR_up_max_hh, HCR_down_max_hh, HCR_sens_hh, plot=F, onlyHCR=F, levlab="Management Lever",Data){

  if(curx<xCP[1]){
    val = yCP[1]
  } else if(curx>xCP[2]){
    val = yCP[2]
  } else {
    val = yCP[1] + (yCP[2]-yCP[1])* (curx-xCP[1])/(xCP[2]-xCP[1])
  }
  ref_update = val/reflev                          # compared to historical ref level
  prop_change = ref_update/rel_lev_val             # size of update relative to last change
  sens_change = exp(log(prop_change)*HCR_sens_hh) # apply sensitivity control parameter
  update = update2 = sens_change
  if((update-1)>HCR_up_max_hh) update2 = 1 + HCR_up_max_hh
  if((update-1)< (-HCR_down_max_hh)) update2 = 1 - HCR_down_max_hh
  delta = update2 * rel_lev_val

  if(plot|onlyHCR){

    plot(c(0,max(xCP,curx)*1.2),c(0,max(yCP,reflev)*1.1),col="white",xlab="Index relative to last historical year",ylab="");grid()
    mtext(levlab,2,line=2.8,cex=0.9)

    lines(c(0,xCP[1]),rep(yCP[1],2),lwd=4); lines(xCP,yCP,lwd=4); lines(c(xCP[2],xCP[2]*1000),rep(yCP[2],2),lwd=4)
    abline(v=curx,col="blue",lwd=2)
    abline(h=reflev,col="orange",lwd=2)
    abline(h=val,col="blue",lwd=2)
    abline(h=rel_lev_val*reflev,col="grey",lwd=2)
    abline(h=update*reflev*rel_lev_val,col="green",lty=2,lwd=1.5)
    abline(h=update2*reflev*rel_lev_val,col="red", lty=2, lwd=2)

    legend('topleft',legend=c(
      paste0("Reference (",round(reflev,2),")"),
      paste0("Previous (",round(rel_lev_val,3),")"),
      paste0("With control rule (",round(ref_update,3),")"),
      paste0("With sensitivity (",round(update*rel_lev_val,3),")"),
      paste0("With change constraints (",round(update2*rel_lev_val,3),")")),
      cex=0.7,text.col=c("orange","grey","blue","green","red"),bty="n")

    xs = max(xCP,curx)*1.2 * c(0.96,0.92,0.88)
    acond = function(x,y)(abs(1-x/y)>0.025)
    if(acond(rel_lev_val*reflev,val))  arrows(xs[1],rel_lev_val*reflev,xs[1],val,0.15,col="blue",lwd=2)
    if(val>1E5){if(acond(val, update*reflev*rel_lev_val)){ arrows(xs[2],val,xs[2], update*reflev*rel_lev_val,0.15,col="green",lwd=2)}}
    if(acond(update*reflev*rel_lev_val, update2*reflev*rel_lev_val)) arrows(xs[3],update*reflev*rel_lev_val,xs[3],update2*reflev*rel_lev_val,0.15,col="red",lwd=2)
    mtext("HCR",line=0.3,font=2,cex=0.9)
  }

  delta
}



get_dims = function(Data){
  nt = nrow(Data@Landings@Value)
  nf = ncol(Data@Landings@Value)
  data.frame(nt=nt,nf=nf)
}


get_season = function(Data){
  rep(1:4,1000)[max(Data@Years)-Data@YearLH+1] # This is the season following the data of this season - so one time step ahead and the correct index for TAC and Effort distribution
}

refhistts = function(Data, Season){ match(Data@YearLH, Data@Years) - (4-Season)}
refrecentts = function(Data){ length(Data@Years) - (3:0)}


Index_smooth<-function(xx, enp.mult, plot=F, plotname=""){
  tofill<-!is.na(xx)
  xx[xx==0]<-1E3
  predout<-rep(NA,length(xx))
  dat<-data.frame(x=1:length(xx),y=log(xx))
  enp.target<-sum(tofill)*enp.mult
  out<-loess(y~x,dat=dat,enp.target=enp.target)
  predout[tofill]<-exp(predict(out))
  if(plot){
    ts_ind = (1:length(xx))[!is.na(xx)]
    plot(xx,type="p",xlab="Year",ylab="Index",ylim=c(0,max(xx,na.rm=T)*1.05),main="",col="black");grid()
    legend('topright',legend=c("Index",paste0("Smoothed (",enp.mult,")")),text.col=c("black","red"),bty='n')
    lines((1:length(xx))[!is.na(xx)],predout[ts_ind],col="#ff000090",lwd=2)
    mtext(plotname,line=0.3,cex=0.9,font=2)
  }
  predout
}










