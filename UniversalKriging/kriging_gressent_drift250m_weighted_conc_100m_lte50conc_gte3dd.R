# Set directory
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
library(lattice)
library(dplyr)
library(tidyverse) # import csv

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
vme2df <- function(obs_data, U) {
  vme_pol_all <- c()
  # Iteration through all raster grid cells 
  for (l in 1:length(unique(obs_data$id_cell))) {
    cell <- unique(obs_data$id_cell)[l]
    pol_conc <- obs_data[which(obs_data$id_cell == cell), 7]$pol
    #N <- length(pol_conc)
    point_count <- as.integer(obs_data[which(obs_data$id_cell == cell), 6]$Point_Count)
    N <- point_count + as.integer(obs_data[which(obs_data$id_cell == cell), 3]$sum_drivedays)
    # VME calculated for grid cells with more than one obs 
    # if we would compute below when N > 1 instead of point_count > 1, we will get NA values, since the s.d. is NA then
    if (point_count > 1) {
      var1 = (obs_data[which(obs_data$id_cell == cell), 1]$std_no2_data / sqrt(N))**2 
      var2 = (U ** 2 / N) * obs_data[which(obs_data$id_cell == cell), 4]$sum_no2_squared
      vme_pol_sens = var1 + var2
    } else {
      vme_pol_sens = -9999 
    } 
    vme_pol_all[l] <- vme_pol_sens
  }
  vme_pol_sens <- vme_pol_all
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
  print(res)
  return(res)
}

#####################################
#           INITIALIZATION          #       
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
indir_drift <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Data/GridData/"
indir_airview <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Data/GoogleAir/Amsterdam/"#
file_drift <- "fishnet_250m_AirView_Utrecht_NO2_LUR_predictions_Drift_AMS_bufferLTE300/fishnet_250m_AirView_Utrecht_NO2_LUR_predictions_Drift_AMS_bufferLTE300.shp"                          # Drift file
file_drift_points <- "fishnet_250m_AirView_Utrecht_NO2_LUR_predictions_Drift_AMS_bufferLTE300.csv"
file_airview <- "AMS_fishnet_NSL_100m_SummarizeWithin_NO2_5Jul_MidPoints_ConcLTE50_gte3dd_points.csv"
outdir <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Scripts/KED/output/" # Directory for output figures
outdir2 <- "C:/Users/Jorane Rogier/Documents/Master/DSthesis/Data/Predictions/"# Directory for output predictions


# Read drift file (which is a .shp file)
# This is done here, to get the parameters from the grid definition
drift <- st_read(paste0(indir_drift,file_drift))
(names(drift)[names(drift) == "mean_lur_p"] <- "pol")
drift_sp <- as(drift, Class="Spatial")
(names(drift_sp)[names(drift_sp) == "mean_lur_p"] <- "pol")


# Parameters for grid definition
xmin=110763.1;
ymin= 476842.9;
xmax= 133513.1;
ymax = 493842.9;
res= 250                                         # Grid Resolution in meters
res_airview = 100 # to compute the VME. 

ncols= (xmax-xmin)/res                                 # Grid dimensions
nrows= (ymax-ymin)/res                                  # Grid dimensions

n_cols_airview = (xmax-xmin)/res_airview
nrows_airview = (ymax-ymin)/res_airview

# Other variables
U_ms <- 0.25        # Uncertainty for mobile sensors (75%)

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

# Define interpolation grid from the LUR model grid.
# This grid does not contain any actual values yet.
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

# Also read the drift points file
drift_points <- read_csv2(file=paste0(indir_drift, file_drift_points))
(names(drift_points)[names(drift_points) == "mean_lur_predictions_and_measured"] <- "pol")
drift_points_df <- as.data.frame(drift_points)
drift_points_df <- subset(drift_points_df, select=-c(OID_,ORIG_FID, Point_Count))
drift_points_df$x <- as.numeric(drift_points_df$POINT_X)
drift_points_df$y <- as.numeric(drift_points_df$POINT_Y)
coordinates(drift_points_df)=~x+y
proj4string(drift_points_df)=CRS_RDnew 
drift_points_df <- spTransform(drift_points_df,CRS_RDnew)
drift_points_df_wgs84 <- spTransform(drift_points_df,CRS_WGS84)

# create raster 
ras_drift=raster(list(x=sort(unique(coordinates(grid_tmp)[,1])),y=sort(unique(coordinates(grid_tmp)[,2])),z=matrix(drift$pol,nrow=length(sort(unique(coordinates(grid_tmp)[,1]))),byrow=F))) # create a raster
plot(ras_drift)

# Print info from drift spatial class
class(drift_sp); summary(drift_sp)


#####################################
#         READ SENSOR DATA          #         
#####################################
# Create grid with lower resolution that will be used for the VME assignment

ext <- raster::extent(xmin, xmax, ymin, ymax) 
grid_s <- st_bbox(ext) %>%
  st_make_grid(cellsize = res_airview, n=c(n_cols_airview, nrows_airview), what = 'polygons') %>%
  st_set_crs(28992)#

# Add cell identifier 
grid_s <- grid_s %>%
  st_sf() %>% 
  mutate(id_cell = seq_len(nrow(.)))


airview_ams <-read_csv2(file=paste0(indir_airview, file_airview))
(names(airview_ams)[names(airview_ams) == "mean_conc_weighted_dd"] <- "pol")

#### Need to randomly subsample 5000 points, otherwise there are too many points for kriging
## airview_ams <- airview_ams[sample(nrow(airview_ams),5000),]
## Not anymore: now subset of original data is taken, where a point has 10+ driving days
colnames(airview_ams)
to.remove <- c("OID_","ROAD_FID", "Shape_Length", "Shape_Area")
airview_ams <- airview_ams[,!names(airview_ams) %in% to.remove]

data_sub <- airview_ams[airview_ams$pol != 0,]

data_sf <- data_sub %>%
  st_as_sf(coords = c('POINT_X','POINT_Y'), crs = 28992, remove=FALSE) %>%
  st_transform(crs = 28992)

# Join data with grid_s to get the sensor observations per grid cell 
data_join <- st_join(data_sf, grid_s)
View(data_join)

data_join <- data_join %>%
  arrange(id_cell)

# Export the ids
unique <- unique(data_join$id_cell)

# Remove the geometry to convert to dataframe to use function
data_join <- st_set_geometry(data_join, NULL)

###############################################################################################################
######################################               MAPPING             ######################################         
###############################################################################################################

print("START MAPPING")

# Set new grid for interpolation
grid <- grid_tmp

### Select sensor data and calculate the Variance of Measurement Errors (VME)

# 1) Data from mobile sensors, compute VME
vme_pol_msens_all <- vme2df(data_join, U_ms)
vme_pol_msens_all = replace(vme_pol_msens_all,vme_pol_msens_all==0,1e-30) # replace null values with very low values => avoid NAN in the colorscale
(uncert_max = max(vme_pol_msens_all))
vme_pol_msens_all <- replace(vme_pol_msens_all, vme_pol_msens_all==-9999, uncert_max*2)#Replace -9999 values with maximum uncertainty
# aggregate
data_ms <- aggregate(.~POINT_X+POINT_Y,data=data_join,mean,na.action=na.pass) # aggregate on location

vme_data <- as.data.frame(vme_pol_msens_all)
# Add id_cell to merge and get the VME for each sensor obs per grid cell
vme_data$id_cell <- unique
data_ms <- merge(data_join, vme_data, by = 'id_cell')
data_ms$vme <- data_ms$vme_pol_msens_all


# Aggregegate sensor data at each measurement position and transform the data frame to a spatial object
# In gressent this is done twice as well, but would expect that this only needs to be done once
data_ms <- aggregate(.~POINT_X+POINT_Y,data=data_ms,mean,na.action=na.pass) # mean over position to avoid duplicate data on a similar location ==> error for kriging

# Convert to spatial object 
data_ms$Long <- data_ms$POINT_X
data_ms$Lat <- data_ms$POINT_Y
coordinates(data_ms) <- ~Long + Lat
proj4string(data_ms) <- CRS_RDnew
data_ms
data_ms <- subset(data_ms, select = -c(POINT_X, POINT_Y, vme_pol_msens_all, id_cell))
tmp <- data_ms@coords; dlon <- tmp[,1]; dlat <- tmp[,2]
data_ms$lon <- dlon; data_ms$lat <- dlat
class(data_ms); summary(data_ms)
data_wgs84 <- spTransform(data_ms, CRS_WGS84)
data_RDnew <- data_ms

kriging_data <- 'MS'


# Look at the number of data points for kriging
nbr_points = length(data_RDnew)
print(paste0("NBR POINTS FOR KRIGING = ",nbr_points))


### Plot drift and data 
# Color and scale for drift
zl_drift <- c(min(drift_points$pol),max(drift_points$pol)+1)
drift_sp[which(drift_sp@data$pol < zl_drift[1])] <- zl_drift[1]# reassign everything < 10 to 10
drift_sp[which(drift_sp@data$pol > zl_drift[2])] <- zl_drift[2]# reassign everything > 30 to 30
#ncol=100 
ncol = ncols# if this is changed, no continuous colors are plotted
brks <- seq(zl_drift[1], zl_drift[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol <- colorRampPalette(mycol)(ncol)

# Color and scale for airview data
zl_airview <- c(min(data_sub$pol),max(data_sub$pol))
data_plt <- data_RDnew
data_plt@data$pol[which(data_plt@data$pol < zl_airview[1])] <- zl_airview[1]+0.1
data_plt@data$pol[which(data_plt@data$pol > zl_airview[2])] <- zl_airview[2]-0.1
cut_pol <- cut(data_plt@data$pol,brks)
col_pol <- mycol_ncol[cut_pol]

# Color and scale for VME
data_plt$vme = data_plt$vme
vme_lim=c(0,uncert_max*2)
ncol2=6
brks2 <- seq(vme_lim[1], vme_lim[2], length.out = ncol2+1)
brkslab2 <- format(brks2, scientific=FALSE, digits=2)
indbrks2 <-  seq(1,length(brks2), by = 15)
mycol2 <- brewer.pal(9,"YlOrRd")
mycol_ncol2 <- colorRampPalette(mycol2)(ncol2)
data_plt@data$vme[which(data_plt@data$vme < vme_lim[1])] <- vme_lim[1]
data_plt@data$vme[which(data_plt@data$vme > vme_lim[2])] <- vme_lim[2]
cut_vme <- cut(data_plt@data$vme,brks2)
col_vme <- mycol_ncol2[cut_vme]



#####
##### Do plot
#####


#spplot(drift_sp, zcol ='pol', xlab='Longitude' ,ylab="Latitude", cex.lab=1.5, cex.axis=1.5)

# plot drift pol concentration, overlaid with target variable
#png(paste(outdir, city,"_drift_",pol,"_",estim_YYYY,"_",kriging_data,"_CONC.png",sep=""),width=900, height=900)
plot(ras_drift, col=mycol_ncol, zlim=zl_drift, breaks=brks, interpolate=FALSE, maxpixels=500000000,
     main="Drift", cex.main=2,
     xlab="Longitude",ylab="Latitude", cex.lab=1.5, cex.axis=1.5,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.5,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.5))
plot(data_plt,pch=21,col='black',bg=col_pol,add=T,cex=1)# plot target variable
dev.off()

# plot drift & VME
#png(paste(outdir, city,"_drift_",pol,"_",estim_YYYY,"_",kriging_data,"_VME.png",sep=""),width=900, height=900)
par(mar=c(7,9,2,9)) # margin bot left top right (need space for the lgd
plot(ras_drift, col=mycol_ncol, zlim=zl_drift, breaks=brks, interpolate=FALSE, maxpixels=500000000,
     main="Drift & VME", cex.main=2,
     xlab="Longitude",ylab="Latitude", cex.lab=1.5, cex.axis=1.5,
     legend.width = 1., legend.shrink=0.7, legend.mar=4,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.5,line=1,font=1, 
                      x = "right", bty="n", inset=c(-0.50,0), xpd = TRUE),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1))
plot(data_plt,pch=21,col='black',bg=col_vme,add=T,cex=0.5)
image.plot(legend.only=TRUE, add=TRUE, horizontal = FALSE, zlim=vme_lim,breaks=brks2,
           col=mycol_ncol2, legend.shrink=.7,legend.width=1.,legend.mar=12,
           legend.args=list(text=expression(paste("VME (", mu, "g/", m^3, ")",)^2),cex=1.5,line=1,font=1,col='black'),
           axis.args=list(at=brks2,labels=format(brks2,digits=2),cex.axis=1,col.axis="black",col="black"))



### Calculate correlation between drift and data
### & Add compute the drift value for the corresponding airview points
data2 <- as.data.frame(data_wgs84)
model_pol <- rep(0, length(data2[,1]))

# compute correlation based on wgs84 projection (for this, needed drift as points)
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


data_RDnew$Drift <- model_pol # add drift to the dataframes
data2$MODEL_POL <- model_pol
#driftd <- drift
#coordinates(driftd)=~lon+lat
#proj4string(driftd)=CRS_L93
#grid$Drift = driftd$drift_pol # add drift to the grid spdataframe
grid$Drift <- drift_sp$pol
(correlation = round(cor(data_RDnew$pol,data_RDnew$Drift),2))

# Plot the correlation
#png(filename=paste0(outdir,"Correlation_",estim_YYYY,".png"), width=600, height=600, type="cairo",bg = "white")
par(mar=c(7,9,2,9)) # margin bot left top right (need space for the lgd
plot(data_RDnew$pol,data_RDnew$Drift, pch=16, col="black",main="a)",cex.main=2, xlab=bquote(AirView ~ NO2 ~ (mu*g/m^3)),ylab=bquote(Drift ~ NO2 ~ (mu*g/m^3)),cex.lab=1.8, cex.axis=1.8)
reg=lm(Drift~pol,data_RDnew)
abline(reg,col="red", lwd=4)
abline(0,1,lwd=4,lty=2,col="red")
r2 = summary(reg)$r.squared
legend("topleft",paste("r=",correlation, "\n", "R2=", round(r2, 2), sep =""),box.col=0,cex=1.8)
dev.off()

### Start kriging

print("START KRIGING")


# 1) Create database for kriging
#data.u.neigh=neigh.create(ndim=2,type=2, nmini=2, nmaxi = 20) #ndim: spatial dimension and type=0: unique neighborhood -> originally 0
data.db=db.create(data_RDnew[,c("lon","lat","pol","vme","Drift")],flag.grid=F,ndim=2,autoname=F)
data.db=db.locate(data.db,'Drift',"f",1)
data.db=db.locate(data.db,'vme',"v",1)
grid.db=db.create(grid[,c("lon","lat","Drift")],flag.grid=F,autoname=F)# grid is also in RDnew
grid.db=db.locate(grid.db,'Drift',"f",1)

# 2) Calculate the residuals between data and drift
data_RDnew$resvcdrift=reg_VC(pol~Drift,data_RDnew) # residuals calculation # sometimes saved at data_RDnew@data$resvcdrift
#data.u.neigh=neigh.create(ndim=2,type=0)
data2.db=db.create(data_RDnew[,c("lon","lat","pol","vme","Drift","resvcdrift")],flag.grid=F,ndim=2,autoname=F) # create database
data2.db=db.locate(data2.db,"resvcdrift","z")

# 3) Calculate the variogram map based on the residuals
vm=vmap.calc(data2.db)
#png(filename=paste0(outdir,"Map_variogram_",city,"_",pol,"_",estim_YYYY,"_",kriging_data,".png"), width=900, height=600, type="cairo",bg = "white") # do plot
par(mar=c(6,9,2,9)) # margin bot left top right (need space for the lgd
plot(vm,npairdw=TRUE,npairpt=TRUE,cex.main=1.2,cex.lab=1.2, cex.axis=1.2)

# 4) Calculate the variogram based on the residuals
vario_res_drift=vario.calc(data2.db, lag=200, nlag=40) # experimental variogram on residuals
# wmode=2: The weight is proportional to the number of pairs and inverse proportional to the distance
vario_res_drift_model=model.auto(vario_res_drift,draw=F,c('Nugget Effect','Spherical','Exponential','Gaussian','Linear','Cubic'),wmode=2,auth.aniso=TRUE,xlab='distance',ylab="variogram") # modeled variogram
#png(filename=paste0(outdir,"Variogram_",estim_YYYY,"h.png"), width=700, height=600, type="cairo",bg = "white")
par(mar=c(7,9,2,9)) # margin bot left top right (need space for the lgd
plot(vario_res_drift, npairdw=TRUE,  npairpt=TRUE,main = "b)",ylim=c(0,200),lwd=3,cex.main=2,xlab="Distance h (m)",ylab="γ(h)",cex.lab=1.8, cex.axis=1.8)
plot(vario_res_drift_model,add=T,col='springgreen3',lwd=3)


# 5) Prepare database for kriging
# type = 0 -> unique neighborhood, too many points for this data
# type = 1 -> bench neighborhood, 'width' specifies the slice around the target point, within which all data points will be to the neigh
# type = 2 -> moving neighborhood, 'nmini' and 'nmaxi' define the min & max nr of points in the neigh
# choose 2, since otherwise there will sometime be target points within neighbors?
#unique.neigh=neigh.create(ndim=2,type=2, nmini=5, nmaxi = 50, radius=c(1000,1000))# TODO: check this -> orginally this:neigh.create(type=0,ndim=2)
unique.neigh=neigh.create(type=0,ndim=2)
#unique.neigh=neigh.create(type=1, width=500)
data2.db=db.locate(data2.db,"pol","z")
data2.db=db.locate(data2.db,"Drift","f",1)# f is default variable for calling function external drift
data2.db=db.locate(data2.db,"vme","v")# v is default variable for calling measurement error variance

krige_ED=kriging(data2.db, grid.db,model=vario_res_drift_model,uc=c("1","f1"),neigh=unique.neigh, calcul = "point") # do kriging

# 6) Allocate results to the grid and save
grid$Pred_EDK=db.extract(krige_ED,'Kriging.pol.estim') # estimation
grid$StDev_EDK=db.extract(krige_ED,'Kriging.pol.stdev') # standard deviation
Pred_EDK=db.locate(krige_ED,9,NA)
StDev_EDK=db.locate(krige_ED,8,NA)
write.table(grid,file=paste(outdir2,"EDK_grid_drift250mNSLLURmin3DrivingDays_lag200_uniqueneigh",city,"_",pol,"_",estim_YYYY,"_",kriging_data,'.csv',sep=''),row.names=F,col.names=T,sep=';')
#grid <- read_csv(file=paste0(outdir2,"EDK_grid_",city,"_",pol,"_",estim_YYYY,"_",kriging_data,'.csv'))
#View(krige_ED)

prediction_df <- data.frame("pred"=ras_Pred@data@values, "x"=grid.db@items$lon, "y"=grid.db@items$lat)
write.csv(prediction_df, paste0(outdir2, "predictions_gressent_drift250m_airview_mte3dd_200lag40neigh_allneighbours.csv"))


### Prediction and mapping error
# Extract prediction and errors, and convert into a matrix to create a raster
Pred_EDK2=matrix(grid$Pred_EDK,nrow=length(sort(unique(coordinates(grid)[,1]))))
StDev_EDK2=matrix(grid$StDev_EDK,nrow=length(sort(unique(coordinates(grid)[,1]))))
min(StDev_EDK2)
min(Pred_EDK2)
max(StDev_EDK2)
max(Pred_EDK2)

#Drift_matrix=matrix(grid$Drift,nrow=length(sort(unique(coordinates(grid)[,1]))))
#coords = which(is.na(Pred_EDK2), arr.ind=TRUE)
#Pred_EDK2[coords] <- Drift_matrix[coords]
#c_pred <- which(Pred_EDK2>100, arr.ind=TRUE)
#c_pred_neg <- which(Pred_EDK2<0, arr.ind=TRUE)
#Pred_EDK2[c_pred]
#Pred_EDK2[c_pred_neg]
#Pred_EDK2[c_pred] <- Drift_matrix[c_pred]
#Pred_EDK2[c_pred_neg] <- Drift_matrix[c_pred_neg]

# Create a raster and define color scale for the prediction
ras_Pred=raster(list(x=sort(unique(coordinates(grid)[,1])),y=sort(unique(coordinates(grid)[,2])),z=Pred_EDK2))
zl_fusion <- c(min(Pred_EDK2),max(Pred_EDK2)+1)
ras_Pred[which(ras_Pred@data@values < zl_fusion[1])] <- zl_fusion[1]
ras_Pred[which(ras_Pred@data@values > zl_fusion[2])] <- zl_fusion[2]
ncol = ncols# if this is changed, no continuous colors are plotted
brks <- seq(zl_fusion[1], zl_fusion[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol <- colorRampPalette(mycol)(ncol)

# Do plot
#png(filename=paste0(outdir,"Pred_",city,"_",pol,"_",estim_YYYY,"_",kriging_data,".png"), width=900, height=900, type="cairo",bg = "white")
par(mar=c(7,9,2,9))
plot(ras_Pred,col=mycol_ncol,scale=1, zlim=zl_fusion,
     xlab='Longitude',  ylab='Latitude', cex.lab=1.2, cex.axis=1.2,
     main="Fused map", cex.main=1.2,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.2,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.2))
dev.off()

# Create a raster and define color scale for the errors
ras_StDev=raster(list(x=sort(unique(coordinates(grid)[,1])),y=sort(unique(coordinates(grid)[,2])),z=StDev_EDK2))
zl_error <- c(min(StDev_EDK2),max(StDev_EDK2)+1)
ras_StDev[which(ras_StDev@data@values < zl_error[1])] <- zl_error[1]
ras_StDev[which(ras_StDev@data@values > zl_error[2])] <- zl_error[2]
ncol=100
brks <- seq(zl_error[1], zl_error[2], length.out = ncol+1)
brkslab <- format(brks, scientific=FALSE, digits=2)
indbrks <-  seq(1,length(brks), by = 15)
mycol_ncol2 <- colorRampPalette(mycol)(ncol)

# Do plot
#png(filename=paste0(outdir,"Stdev_",city,"_",pol,"_",estim_period,"_",kriging_data,".png"), width=900, height=900, type="cairo",bg = "white")
plot(ras_StDev, scale=1, col=mycol_ncol2, zlim=zl_error,
     xlab='Longitude', ylab='Latitude', cex.lab=1.2, cex.axis=1.2,
     main="Error map", cex.main=1.2,
     legend.width = 1.2, legend.shrink=0.75,
     legend.args=list(text=expression(paste("(", mu, "g/", m^3, ")")),cex=1.2,line=1,font=1),
     axis.args = list(at=brks[indbrks], labels=brkslab[indbrks], cex.axis=1.2))


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
save(data2,file=paste0(outdir2,"EDK_grid_drift250mNSLLURmin3DrivingDays_lag200_uniqueneigh_",city,"_",pol,"_",estim_YYYY,"_",kriging_data,".Rda"))

# Do plot
png(filename=paste0(outdir,"CV_EDK_grid_drift250mNSLLURmin3DrivingDays_lag200_uniqueneigh_",city,"_",pol,"_",estim_YYYY,"_",kriging_data,"_VME.png"), width=1200, height=800, type="cairo",bg = "white")
par(mfrow=c(1,2), cex.lab=2, cex.main=1.5, cex.axis=2, mar=c(5,5,5,5))
plot(data2$pol, data2$estim_KDE, xlim=c(floor(min(data2$pol)),floor(max(data2$pol))+1), ylim=c(floor(min(data2$pol)),floor(max(data2$pol))),main =" Scatter plot observations - estimations", ylab =" Estimations (CV)",xlab =" Observations",pch =15)
lines(data$pol,data$pol)
legend("topright",box.col=1,paste("cor=",round(cor(data2$pol,data2$estim_KDE),2)),cex =1)
hist(data2.db[,12] ,col ="slategray", border ="black", main ="Residuals from cross validation",xlab ="Residuals",ylab ="Density",prob =T)
dev.off()

final.cor <- round(cor(data2$pol,data2$estim_KDE),2)
final.R2 <- final.cor**2
n <- length(data2$pol)
k <- 1
adj.R2 <- 1 -((1-final.R2)*(n-1)/(n-k-1))

