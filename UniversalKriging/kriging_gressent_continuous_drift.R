# Set directory
setwd("/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/03_Mapping/") # !!! MUST BE UPDATED BY THE USER !!!

# Import libraries 
library(RGeostats) # !!! MUST BE UPDATED BY THE USER !!!
library(maptools);library(RColorBrewer)
library(fields);library(raster)
library(rgdal);library(gstat)
library(XML);library(akima)
library(rgeos);library(sp)
library(spacetime);library(chron)
library(plotrix); library(Rcpp)
library(RcppArmadillo); library(MASS)
library(foreach)
library(doParallel)
library(sf)
#################################################################################
#               AQ MAPPING USING SENSOR FIXED AND MOBILE DATA                   #
#           Kriging with an external drift, use RGeostats packages              #               
#                      Created 23/10/2018 updated 18/05/2020                    #
#            Author: Alicia Gressent (INERIS) alicia.gressent@ineris.fr         #
#################################################################################

#####################################
#            FUNCTIONS              #          
#####################################

# Function to calculate the Variance of Measurement Errors (VME) from sensor data
vme2df <- function(obs_data,U) {
  unique_data <- aggregate(.~lat+lon+id_sensor,data=obs_data,unique,na.action=na.pass)
  vme_pol_all <- c()
  for (l in 1:length(unique_data[,1])){
    loni=unique_data[l,2]
    lati=unique_data[l,1]
    pol_conc = obs_data[which(obs_data$lat == lati & obs_data$lon == loni),4]# 4 -< represents no2 measurements
    N = length(pol_conc)
    if (length(pol_conc) > 1){
      var1 <- (sd(pol_conc)/sqrt(N))**2                              # dispersion
      var2 <- (U**2/N)*sum(pol_conc[2:length(pol_conc)]**2)          # instrument uncertainties
      vme_pol_sens <- var1 + var2                                    # sum of variances
    }else{
      vme_pol_sens <- -9999
    }          
    vme_pol_all[l] = vme_pol_sens
  }
  vme_pol_sens = vme_pol_all
  return(vme_pol_sens)
}

# Function for residuals calculation of the regression between the estimation and the drift (cross validation)
reg_VC<-function(formula,df) {
  res=NULL
  for (i in 1:nrow(df)) {
    reg=lm(formula,data=df[-i,])
    restmp=predict(reg,newdata=df[i,])-df$pol[i]
    res=c(res,restmp)
  }
  return(res)
}

#####################################
#           INITIALIZATION          # !!! MUST BE UPDATED BY THE USER !!!          
#####################################

# Define projections
CRS_RDnew <- st_crs(28992)$proj4string # RD new
CRS_WGS84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Situation
city="Amsterdam"                   # Estimation location
pol="NSL_LUR_NO2_predictions"      # Pollutant
model="NSL-LUR"                    # Dispersion model
drift_type="2019-annual-mean"   # Definition of the drift

# Estimation time
estim_YYYY <- "2019" # Extract year


# Directory paths
indir_drift <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Data/Predictions/"
indir_airview <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Data/GoogleAir/"
file_drift <- "fishnet_AMS_NSL_500m_SummarizeWithin_NSL_LUR_predictions_.shp"                                                            # Drift file
file_airview <-"googleLabsAirQualityAmsterdam_latlong.csv"
outdir <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Scripts/KED/output/" # Directory for output figures
#outdir2 <-"/ccc/work/cont004/ineris/gressena/microcapteurs_carto/SESAM/OUTPUTS/"        # Directory for output files

# Read drift file
drift <- st_read(paste0(indir_drift,file_drift))
ggplot() + 
  geom_sf(data = drift, size = 3, color = "black", fill = "cyan1") + 
  ggtitle("drift") + 
  coord_sf()
drift_sp <- as(drift, Class="Spatial")
drift

drift_points <- read_csv2(file=paste0(indir_drift, "fishnet_AMS_NSL_500m_SummarizeWithin_NSL_LUR_predictions_spatialPoints.csv"))
colnames(drift_points)
drift_points$

drift_points_df <- as.data.frame(drift_points)
drift_points_df <- subset(drift_points_df, select=-c(OID_,mean_nsl_lur_no2_predictions,Point_Count,ORIG_FID))
colnames(drift_points_df)
drift_points_df$x <- as.numeric(drift_points_df$POINT_X)
drift_points_df$y <- as.numeric(drift_points_df$POINT_Y)
coordinates(drift_points_df)=~x+y
proj4string(drift_points_df)=CRS_RDnew 
drift_points_df <- spTransform(drift_points_df,CRS_RDnew)
drift_points_df_wgs84 <- spTransform(drift_points_df,CRS_WGS84)

data <- subset(data, select=-c(lat,lon))
tmp = data@coords; dlon=tmp[,1]; dlat=tmp[,2]
data$lon <- dlon; data$lat <- dlat
class(data); summary(data)
#data <- data[data$lon >= X1 & data$lon <= X2 & data$lat >= Y1 & data$lat <= Y2,] # subset data to the domain limits
data_wgs84 <- spTransform(data,CRS_WGS84)

drift_points_matrix <- data.matrix(drift_points_df)

# Parameters for grid definition
#X1=345000; X2=365000; Y1=6680000; Y2=6701000    # Estimation domain limits
xmin=110763.1;
ymin= 476842.9;
xmax= 133763.1;
ymax = 493842.9;
res= 500                                         # Grid Resolution in meters

ncols= (xmax-xmin)/res                                 # Grid dimensions
nrows= (ymax-ymin)/res                                  # Grid dimensions

# Other variables
U_ms <- 0.75        # Uncertainty for mobile sensors (75%)
U_fs <- 0.50        # Uncertainty for fixed sensors  (50%)

# Constants
R <- 6371e3         # Medium earth radius in meters

# Storage
nbr_pts_all <-c()   # Number of data points for kriging

# Color palette for mapping
d=rbind(c(38,48,132),c(61,99,174),c(114,201,195),c(220,225,30),c(240,78,34),c(133,22,24))
palette=colorRampPalette(rgb(d[,1],d[,2],d[,3],maxColorValue=255))(128)
mycol=palette

#####################################
#               GRID                #         
#####################################

# Define interpolation grid from the dispersion model grid
print("CREATE GRID")

X<-c() # x vector
for (i in 1:ncols){if (i==1){X[i]=xmin}else{X[i]=X[i-1]+res}}
Y<-c() # y vector
for (i in 1:nrows){if (i==1){Y[i]=ymin}else{Y[i]=Y[i-1]+res}}
grid1=expand.grid(X,Y) # create grid from x and y columns
data1 <- rep(0,length(grid1)) # first attribute pseudo-data
grid1$data1 <- data1
grid1 = subset(grid1, Var1 >= xmin & Var1 <= xmax & Var2 >=ymin & Var2 <= ymax) # subset full grid to the domain geographical limits if necessary
coordinates(grid1)=~Var1+Var2 # convert grid to a spatial object
proj4string(grid1)=CRS_RDnew # define grid projection
grid_coords = grid1@coords # extract coordinates

# Define a temporary grid (grid_tmp) as the grid to use for the interpolation
grid_tmp <- data.frame(lon=seq(1:length(grid_coords[,1])),
                       lat=seq(1:length(grid_coords[,1])),
                       Drift=seq(1:length(grid_coords[,1])))
grid_tmp$lon <- grid_coords[,1]
grid_tmp$lat <- grid_coords[,2]
grid_tmp$Drift = grid_tmp$Drift*1 # add pseudo drift to the grid (no matter the values)
grid_tmp$Long <- grid_tmp$lon; grid_tmp$Lat <- grid_tmp$lat
coordinates(grid_tmp)=~Long+Lat
proj4string(grid_tmp)=CRS_RDnew # define projection

#####################################
#             READ DRIFT            #         
#####################################

print("READ DRIFT")

# names(drift)[names(drift) == "NSL_LUR_x"] <- "lon"
# names(drift)[names(drift) == "NSL_LUR_y"] <- "lat"
# #drift_sub = subset(drift, lon >= X1 & lon <= X2 & lat >= Y1 & lat <= Y2) # subset full grid to the data geographical limits

ras_drift=raster(list(x=sort(unique(coordinates(grid_tmp)[,1])),y=sort(unique(coordinates(grid_tmp)[,2])),z=matrix(drift$mean_nsl_l,nrow=length(sort(unique(coordinates(grid_tmp)[,1]))),byrow=F))) # create a raster
plot(ras_drift)

# # Save drift as a spatial dataframe
# spmodel <- drift_sub
# spmodel$Long <- spmodel$lon; spmodel$Lat <- spmodel$lat
# coordinates(spmodel)=~Long+Lat
# proj4string(spmodel)=CRS_L93
# spmodel <- spTransform(spmodel,CRS_WGS84)
# spmodel <- subset(spmodel, select=-c(lat,lon))
# tmp = spmodel@coords; dlon=tmp[,1]; dlat=tmp[,2]
# spmodel$lon <- dlon
# spmodel$lat <- dlat


class(drift_sp); summary(drift_sp)
names(drift_sp)[names(drift_sp) == "mean_nsl_l"] <- "pol" 



#####################################
#         READ SENSOR DATA          #         
#####################################

# Data from mobile sensors
print("READ MOBILE AIRVIEW DATA")
airview_df <- read.csv2(paste0(indir_airview, file_airview), sep=",")

# remove 0 (-999999) values
airview_df <- airview_df[airview_df$mixed_no2 > 0,]

# make all variables numeric
#airview_df <- sapply(airview_df, as.numeric)
airview_df <- as.data.frame(airview_df)

# load(paste(indir2,file_ms,sep=","))
# data_ms_preproc <- data_preproc
# names(data_ms_preproc)[names(data_ms_preproc) == "pol"] <- "pol_ori"
# data_ms_preproc = subset(data_ms_preproc, selec=-c(run,pm25)) # !!! MUST BE UPDATED BY THE USER !!! 
# data_ms_preproc = na.omit(data_ms_preproc)

#####################################
#               MAPPING             #         
#####################################

print("START MAPPING LOOP OVER HOURS")

  
# Set new grid for interpolation
grid <- grid_tmp

### Select sensor data and calculate the Variance of Measurement Errors (VME)

# 1) Data from mobile sensors
colnames(airview_df)
names(airview_df)[names(airview_df) == "long"] <- "lon"
names(airview_df)[names(airview_df) == "ogc_fid"] <- "id_sensor"
vme_pol_msens_all=vme2df(airview_df,U_ms)

vme_pol_msens_all = replace(vme_pol_msens_all,vme_pol_msens_all==0,1e-30) # replace null values with very low values => avoid NAN in the colorscale
uncert_max = max(vme_pol_msens_all)
vme_pol_msens_all <- replace(vme_pol_msens_all, vme_pol_msens_all==-9999, uncert_max*2)#Replace -9999 values with maximum uncertainty
colnames(airview_df)
airview_df = subset(airview_df, selec=-c(mixed_ufp, mixed_bc, latlong, X, X_id))
colnames(airview_df)
airview_df$mixed_no2 <- as.numeric(airview_df$mixed_no2)
#data_ms <- aggregate(.~lat+lon,data=airview_df,mean,na.action=na.pass) # aggregate on location
#data_ms = subset(data_ms, selec=-c(id_sensor))
#names(data_ms)[names(data_ms) == "pol_bgcorr"] <- "pol"
#names(data_ms)[names(data_ms) == "datetime"] <- "date"
data_ms <- airview_df
data_ms$vme <- vme_pol_msens_all


### Define data for kriging and concatenate into a final dataframe
data <- data_ms; kriging_data <- 'MS'

### Aggregate sensor data at each measurement position and transform the dataframe to a spatial object
#data <- aggregate(.~lat+lon,data=data_tmp,mean,na.action=na.pass) # mean over position to avoid duplicate data on a similar location ==> error for kriging
data$Long <- as.numeric(data$lon); data$Lat <- as.numeric(data$lat)
coordinates(data)=~Long+Lat
proj4string(data)=CRS_RDnew 
#data <- spTransform(data,CRS_RDnew)
data <- subset(data, select=-c(lat,lon))
tmp = data@coords; dlon=tmp[,1]; dlat=tmp[,2]
data$lon <- dlon; data$lat <- dlat
class(data); summary(data)
#data <- data[data$lon >= X1 & data$lon <= X2 & data$lat >= Y1 & data$lat <= Y2,] # subset data to the domain limits
data_wgs84 <- spTransform(data,CRS_WGS84)


# Look at the number of data points for kriging
nbr_points = length(data)
print(paste0("NBR POINTS FOR KRIGING = ",nbr_points))

### Plot drift and data 

# Color and scale for drift
zl <- c(10,30)
drift_sp[which(drift_sp@data$pol < zl[1])] <- zl[1]
drift_sp[which(drift_sp@data$pol > zl[2])] <- zl[2]
#ncol=100 
ncol = ncols
brks <- seq(zl[1], zl[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol <- colorRampPalette(mycol)(ncol)

# Color and scale for data
data_plt <- data
data_plt@data$mixed_no2[which(data_plt@data$mixed_no2 < zl[1])] <- zl[1]+0.1
data_plt@data$mixed_no2[which(data_plt@data$mixed_no2 > zl[2])] <- zl[2]-0.1
cut_pol <- cut(data_plt@data$mixed_no2,brks)
col_pol <- mycol_ncol[cut_pol]

# Color and scale for VME
data_plt$vme = data_plt$vme
vme_lim=c(0,1500)
ncol2=9
brks2 <- seq(vme_lim[1], vme_lim[2], length.out = ncol2+1)
brkslab2 <- format(brks2, scientific=FALSE, digits=2)
indbrks2 <-  seq(1,length(brks2), by = 15)
mycol2 <- brewer.pal(9,"YlOrRd")
mycol_ncol2 <- colorRampPalette(mycol2)(ncol2)
data_plt@data$vme[which(data_plt@data$vme < vme_lim[1])] <- vme_lim[1]
data_plt@data$vme[which(data_plt@data$vme > vme_lim[2])] <- vme_lim[2]
cut_vme <- cut(data_plt@data$vme,brks2)
col_vme <- mycol_ncol2[cut_vme]

# Do plot
png(paste(outdir, city,"_drift_",pol,"_",estim_YYYY,"_",kriging_data,"_CONC.png",sep=""),width=900, height=900)
plot(drift_sp, col=mycol_ncol, zlim=zl, breaks=brks, interpolate=FALSE, maxpixels=500000000,
     main="test", cex.main=2,
     xlab="Longitude",ylab="Latitude", cex.lab=1.5, cex.axis=1.5,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.5,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.5))
plot(data_plt,pch=21,col='black',bg=col_pol,add=T,cex=2)
dev.off()


png(paste(outdir, city,"_drift_",pol,"_",estim_YYYY,"_",kriging_data,"_VME.png",sep=""),width=900, height=900)
plot(drift_sp, col=mycol_ncol, zlim=zl, breaks=brks, interpolate=FALSE, maxpixels=500000000,
     main="Drift & VME", cex.main=2,
     xlab="Longitude",ylab="Latitude", cex.lab=1.5, cex.axis=1.5,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.5,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.5))
plot(data_plt,pch=21,col='black',bg=col_vme,add=T,cex=0.5)
image.plot(legend.only=TRUE, add=TRUE, horizontal = FALSE, zlim=vme_lim,breaks=brks2,
           col=mycol_ncol2, legend.shrink=.7,legend.width=1.,legend.mar=8,
           legend.args=list(text=expression(paste("VME (", mu, "g/", m^3, ")",)^2),cex=1,line=1,font=1,col='black'),
           axis.args=list(at=brks2,labels=format(brks2,digits=2),cex.axis=1,col.axis="black",col="black"))
dev.off()

### Calculate correlation between drift and data
### & Add compute the drift value for the corresponding airview points
data2 <- as.data.frame(data)
model_pol <- rep(0, length(data2[,1]))

for (ll in 1:length(data_wgs84)){
  data_wgs84_tmp = data_wgs84[ll,]
  dist_vector = spDistsN1(drift_points_df_wgs84,data_wgs84_tmp,longlat = TRUE)
  dist_vector = dist_vector*1e3 # km to m
  dist_min=min(dist_vector)
  idx=which(dist_vector <= dist_min)
  tmp=drift_sp[idx,]
  tmp2=as.data.frame(tmp)
  pol_tmp2=mean(tmp2$pol)
  model_pol[ll]=pol_tmp2
}


data$Drift <- model_pol # add drift to the dataframes
data2$MODEL_POL <- model_pol
driftd <- drift
#coordinates(driftd)=~lon+lat
#proj4string(driftd)=CRS_L93
#grid$Drift = driftd$drift_pol # add drift to the grid spdataframe
correlation = round(cor(data$mixed_no2,data$Drift),2)
correlation

# Do plot
png(filename=paste0(outdir,"Correlation_",estim_YYYY,".png"), width=600, height=600, type="cairo",bg = "white")
par(mar=c(7,9,2,9)) # margin bot left top right (need space for the lgd
plot(data$mixed_no2,data$Drift, pch=16, col="black",main="a)",cex.main=2, xlab=bquote(LCS ~ .(pol) ~ (mu*g/m^3)),ylab=bquote(Drift ~ .(pol) ~ (mu*g/m^3)),cex.lab=1.8, cex.axis=1.8)
reg=lm(mixed_no2~Drift,data)
abline(reg,col="red", lwd=4)
abline(0,1,lwd=4,lty=2,col="red")
legend("topleft",paste("R=",correlation,sep =""),box.col=0,cex=1.8)
dev.off()

### Start kriging

print("START KRIGING")


# 1) Create database for kriging
data.u.neigh=neigh.create(ndim=2,type=0) #ndim: spatial dimension and type=0: unique neighborhood
data.db=db.create(data[,c("lon","lat","mixed_no2","vme","Drift")],flag.grid=F,ndim=2,autoname=F)
grid.db=db.create(grid[,c("lon","lat","Drift")],flag.grid=F,autoname=F)
data.db=db.locate(data.db,'Drift',"f",1)
data.db=db.locate(data.db,'vme',"v",1)
grid.db=db.locate(grid.db,'Drift',"f",1)
unique.neigh=neigh.create(type=0,ndim=2)

# 2) Calculate the residuals between data and drift
data$resvcdrift=reg_VC(mixed_no2~Drift,data) # residuals calculation
data.u.neigh=neigh.create(ndim=2,type=0)
data2.db=db.create(data[,c("lon","lat","pol","vme","Drift","resvcdrift")],flag.grid=F,ndim=2,autoname=F) # create database
data2.db=db.locate(data2.db,"resvcdrift","z")

# 3) Calculate the variogram map based on the residuals
vm=vmap.calc(data2.db)
png(filename=paste0(outdir,"Map_variogram_",city,"_",pol,"_",estim_period,"_",kriging_data,".png"), width=900, height=600, type="cairo",bg = "white") # do plot
par(mar=c(6,9,2,9)) # margin bot left top right (need space for the lgd
plot(vm,npairdw=TRUE,npairpt=TRUE,cex.main=1.2,cex.lab=1.2, cex.axis=1.2)
dev.off()

# 4) Calculate the variogram based on the residuals
vario_res_drift=vario.calc(data2.db) # experimental variogram on residuals
vario_res_drift_model=model.auto(vario_res_drift,draw=F,c('Nugget Effect','Spherical','Exponential','Gaussian','Linear','Cubic'),wmode=2,auth.aniso=TRUE,xlab='distance',ylab="variogram") # modeled variogram
png(filename=paste0(outdir,"Variogram_",estim_date2,"_",estim_HH_end,"h.png"), width=700, height=600, type="cairo",bg = "white")
par(mar=c(7,9,2,9)) # margin bot left top right (need space for the lgd
plot(vario_res_drift, npairdw=TRUE,  npairpt=TRUE,main = "b)",ylim=c(0,35),lwd=3,cex.main=2,xlab="Distance h (m)",ylab="γ(h)",cex.lab=1.8, cex.axis=1.8)
plot(vario_res_drift_model,add=T,col='springgreen3',lwd=3)
dev.off()

# 5) Prepare database for kriging
data2.db=db.locate(data2.db,"pol","z")
data2.db=db.locate(data2.db,"Drift","f",1)
data2.db=db.locate(data2.db,"vme","v")
krige_ED=kriging(data2.db,grid.db,model=vario_res_drift_model,uc=c("1","f1"),neigh=unique.neigh) # do kriging

# 6) Allocate results to the grid and save
grid$Pred_EDK=db.extract(krige_ED,'Kriging.pol.estim') # estimation
grid$StDev_EDK=db.extract(krige_ED,'Kriging.pol.stdev') # standard deviation
Pred_EDK=db.locate(krige_ED,9,NA)
StDev_EDK=db.locate(krige_ED,8,NA)
write.table(grid,file=paste(outdir2,"EDK_grid_",city,"_",pol,"_",estim_period,"_",kriging_data,'.csv',sep=''),row.names=F,col.names=T,sep=',')

### Prediction and mapping error

# Extract prediction and errors, and convert into a matrix to create a raster
Pred_EDK2=matrix(grid$Pred_EDK,nrow=length(sort(unique(coordinates(grid)[,1]))))
StDev_EDK2=matrix(grid$StDev_EDK,nrow=length(sort(unique(coordinates(grid)[,1]))))

# Create a raster and define color scale for the prediction
ras_Pred=raster(list(x=sort(unique(coordinates(grid)[,1])),y=sort(unique(coordinates(grid)[,2])),z=Pred_EDK2))
if (pol=="PM10"){ zl <- c(10,30)}
ras_Pred[which(ras_Pred@data@values < zl[1])] <- zl[1]
ras_Pred[which(ras_Pred@data@values > zl[2])] <- zl[2]

# Do plot
png(filename=paste0(outdir,"Pred_",city,"_",pol,"_",estim_period,"_",kriging_data,".png"), width=900, height=900, type="cairo",bg = "white")
plot(ras_Pred,col=mycol_ncol,scale=1, zlim=zl,
     xlab='Longitude',  ylab='Latitude', cex.lab=1.2, cex.axis=1.2,
     main="Fused map", cex.main=1.2,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.2,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.2))
dev.off()

# Create a raster and define color scale for the errors
ras_StDev=raster(list(x=sort(unique(coordinates(grid)[,1])),y=sort(unique(coordinates(grid)[,2])),z=StDev_EDK2))
if (pol=="PM10"){ zl <- c(0,10)}
ras_StDev[which(ras_StDev@data@values < zl[1])] <- zl[1]
ras_StDev[which(ras_StDev@data@values > zl[2])] <- zl[2]
ncol=100
brks <- seq(zl[1], zl[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol2 <- colorRampPalette(mycol)(ncol)

# Do plot
png(filename=paste0(outdir,"Stdev_",city,"_",pol,"_",estim_period,"_",kriging_data,".png"), width=900, height=900, type="cairo",bg = "white")
plot(ras_StDev, scale=1, col=mycol_ncol2, zlim=zl,
     xlab='Longitude', ylab='Latitude', cex.lab=1.2, cex.axis=1.2,
     main="Error map", cex.main=1.2,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.2,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.2))
dev.off()


#####################################
#       MAPPING PERFORMANCE         #         
#####################################

### Kriging performance

print("KRIGING PERFORMANCE")

# Assess kriging performance by cross validation and save for postprocessing
data2.db=db.locate(data2.db,"pol","z")
data2.db=xvalid(data2.db,vario_res_drift_model,neigh.create(ndim=2,type=0),uc=c('1','f1'),radix='KDE_drift') # Perform cross-validation
data2.db=db.locate(data2.db,"KDE_drift.pol.stderr",NA)
data2$estim_KDE=data2$pol+db.extract(data2.db,'KDE_drift.pol.esterr')
save(data2,file=paste0(outdir2,"EDK_data4perf_",city,"_",pol,"_",estim_period,"_",kriging_data,".Rda"))

# Do plot
png(filename=paste0(outdir,"CV_",city,"_",pol,"_",estim_period,"_",kriging_data,"_VME.png"), width=1200, height=800, type="cairo",bg = "white")
par(mfrow=c(1,2), cex.lab=2, cex.main=2, cex.axis=2, mar=c(5,5,5,5))
plot(data2$pol, data2$estim_KDE, xlim=c(floor(min(data2$pol)),floor(max(data2$pol))+1), ylim=c(floor(min(data2$pol)),floor(max(data2$pol))+1),main =" Scatter plot observations - estimations", ylab =" Estimations (CV)",xlab =" Observations",pch =15)
lines(data$pol,data$pol)
legend("topleft",box.col=0,paste("cor=",round(cor(data2$pol,data2$estim_KDE),2)),cex =2)
hist(data2.db[,12] ,col ="slategray", border ="black", main ="Residuals from cross validation",xlab ="Residuals",ylab ="Density",prob =T)
dev.off()
  

