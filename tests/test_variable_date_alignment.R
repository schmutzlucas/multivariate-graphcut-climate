source("R/load.R");load_mvgc();require_packages("ncdf4")
root<-tempfile("variable_dates_");dir.create(root)
model<-"FIXTURE";variables<-c("pr","tas","psl")
write_variable<-function(v,start_year) {
 folder<-file.path(root,model,v);dir.create(folder,recursive=TRUE)
 path<-file.path(folder,paste0(v,"_",model,"_",start_year,"0101-19531231.nc"))
 dims<-list(ncdf4::ncdim_def("lon","degrees_east",0:1),ncdf4::ncdim_def("lat","degrees_north",0:1),
  ncdf4::ncdim_def("time",paste0("days since ",start_year,"-01-01 00:00:00"),
   seq_len((1954-start_year)*365L)-.5,calendar="365_day"))
 n<-ncdf4::nc_create(path,ncdf4::ncvar_def(v,switch(v,pr="mm day-1",tas="K",psl="Pa"),dims,missval=-999,prec="double"))
 # Include negative/zero/positive/missing precipitation to exercise every backend.
 count<-(1954-start_year)*365L
 values<-if(v=="pr")array(rep(rep(c(-1e-6,0,.5,1,2,NA_real_),length.out=count),each=4),c(2,2,count))else array(.5,c(2,2,count))
 ncdf4::ncvar_put(n,v,values);ncdf4::nc_close(n)
}
write_variable("pr",1950);write_variable("tas",1949);write_variable("psl",1950)
bounds<-array(c(rep(0,12),rep(1,12)),c(2,2,3,2))
common<-list(model=model,variables=variables,data_dir=root,year_present=1951L,year_future=1952L,
 lon=0:1,lat=0:1,range_var=bounds,season_start_md=c(1,1),season_end_md=c(12,31),
 pdf2_names=c("pr_tas","pr_psl","tas_psl"),nbins3d=8L,nbins2d=16L,nbins1d=32L,
 strict_time_checks=TRUE,season_selection_mode="auto",time_decode_mode="auto",requested_year_policy="warn",
 diagnostics=FALSE,verbose=FALSE,compute_hist_nd_fun=compute_histND)
a<-do.call(compute_model_pdfs_fractional_legacy,common)
b<-do.call(compute_model_pdfs_fractional_optimized,c(common,list(spatial_block_size=1L)))
c<-do.call(compute_model_pdfs_fractional_optimized,c(common,list(spatial_block_size=NULL)))
stopifnot(identical(a$pdf3_present,b$pdf3_present),identical(a$pdf3_future,b$pdf3_future),
 all(apply(b$pdf3_present,c(1,2),sum)==1))
for(component in c("pdf3_present","pdf3_future","pdf2_present","pdf2_future","pdf1_present","pdf1_future","mean_present","mean_future"))
 stopifnot(identical(a[[component]],b[[component]]),identical(a[[component]],c[[component]]))
raw<-rep(c(-1e-6,0,.5,1,2,NA_real_),length.out=1460)[366:730]
expected_mean<-mean(precipitation_qc(raw),na.rm=TRUE)
stopifnot(max(abs(a$mean_present$pr-expected_mean))<1e-15)
plans<-prepare_k_sensitivity_plans(model,variables,root,1951:1952,0:1,0:1,verbose=FALSE)
stopifnot(length(plans[[model]]$selected_year)==730L)
block<-read_k_sensitivity_block(plans[[model]],variables,1L)
stopifnot(all(block$pr>=0,na.rm=TRUE),any(block$pr==0,na.rm=TRUE),anyNA(block$pr))
cat("PASS legacy/full/blocked PMFs, original-scale QC means and k-sensitivity agree across unequal file start years\n")
