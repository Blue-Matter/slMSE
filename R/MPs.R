
# fake data MP

fakeMP = function(Data, repyr = 2028, outdir = "C:/temp", filenam = "Example_Data", Eff=0.5){

  ad=Advice(Effort=Eff)
  if(max(Data@Years)==repyr){
    Example_Data = Data
    outfile =  paste0(outdir,"/",filenam,".rds")
    saveRDS(Example_Data, outfile)
    cat(paste0("Example data file with data in year ",repyr, " outputted to: ",outfile, "\n"))
    stop()
  }
  ad

}
class(fakeMP) = "mp"




# Generic MPs

Eff = function(Data, rel_eff = 1){
  Advice(Effort = rel_eff)
}
class(Eff) = "mp"



# Index based MPs

