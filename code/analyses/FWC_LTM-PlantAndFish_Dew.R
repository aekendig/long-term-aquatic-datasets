
# Amy,
# See below for my approach to combining LTM plant data into a single data frame
# I combined them into different data frames by year, but you can rbind all of those together
# by doing something like:
#   do.call(rbind, list(csv1, csv2, ....))
# The way I approach it is to key out the plant codes and replace them with the scientific name. If 
# there is no match in the key, it will keep the code or name that is there.  
# After mashing together all of the data frames, creating new columns if they didn't exist in the last 
# data frame (done with ldply), then I melt it all down into a long format data frame that has one row 
# per plant entry. That may not be what is useful for you, but it seems like the cleanest way to analyze 
# the data.  If you need to match up common sites for different plant species, you could maybe assign an id by 
# just going df$id<-1:nrow(df) right before you melt the dataframe.
# 
# If you want to meet up and run through any of this, just let me know.





library(plyr)
library(dplyr)
library(FNN)
library(rgeos)
library(raster)
library(reshape2)


#pointfiles<-'S:/HSC/IPM/Biobase/ShinyAppData'

#all files are in a folder according to the year sampled.  Get each of those years (2015,2016,2017)
years<-list.dirs(pointfiles, full.names=F, recursive=F)

#Make a list of all files in each year folder for each of the data types (point intercept, Biovolume, bathymetry, bottom hardness)
point_fi<-lapply(years, function(x){ list.files(path=paste(pointfiles,'/',x,sep=""),pattern='\\.csv', full.names=T)})
pointfiles<-point_fi[[3]]  #this is the latest year



##############
# we just start by putting all of one year's point intercept csv's in a single folder:
# pointfol<-"C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/Projects/BioBase/ShinyAppData/2020"
pointfol <- "original-data/fwri-point-intercept-data/2020"

## i used pattern= here to make sure I was getting just the point intercept csvs and not any other ones
pointfiles<-list.files(pointfol, full.names=T, pattern="PlantData")



## this ldply function runs everything inside and then mashes all of the results into a single data.frame
  ## it is useful because it will bind data that has the same column name, but if there is a new one, it will 
    # create a new column (unlike rbind() )

# AK: below doesn't work with the folders of data I have, can manually check for differences from my data processing

point<-ldply(pointfiles ,function(x){
  
  biopoint<-read.csv(x) #read in the point data

  # rename all headers that are associated with lat/longs as just X and Y
  colnames(biopoint)[tolower(names(biopoint)) %in% c("longitude","long","lng","x")]<-"X"
  colnames(biopoint)[tolower(names(biopoint)) %in% c("latitude","lat","y")]<-"Y"
  
  biopoint<-biopoint[!is.na(biopoint$X),] #get rid of NA lines (YOU MAY WANT TO SKIP THIS AND SEE WHAT IS NA, JUST TO CHECK WE'RE NOT MISSING ANYTHING IMPORTANT)
  
  #replace missing values with zeros.  These NAs actually represent 0's, so we will replace them as such.  If a site truly was missed, there 
    #are the No_Access columns in the data to say that the site was "Not Accessible" that year. Some of them say why they were not accessible like NA_ISLAND
  0->biopoint[is.na(biopoint)] 
  
  
  ### here are all of the codes that I have seen be written in the column name, instead of having the code written. We will replace them
  names(biopoint)[names(biopoint)=="DogFennel"]<-"DOFE"
  names(biopoint)[names(biopoint)=="Dodder"]<-"DODD"
  names(biopoint)[names(biopoint)=="Milkweed"]<-"MIWE"
  names(biopoint)[names(biopoint)=="Sphagnum.Moss"]<-"MOSS"
  names(biopoint)[names(biopoint)=="Shade.Mudflower"]<-"SHMU"
  names(biopoint)[names(biopoint)=="Submersed.Sagittaria"]<-"SAGI"
  names(biopoint)[names(biopoint)=="SSAG"]<-"SAGI"
  names(biopoint)[names(biopoint)=="VirginiaButtonweed"]<-"BUWE"
  names(biopoint)[names(biopoint) %in% c("RedMaple","Red.Maple", "Maple")]<-"ACRU"
  names(biopoint)[names(biopoint) %in% c("Oak", "OAK")]<-"QUER"
  names(biopoint)[names(biopoint) %in% c("Pigweed", "pigweed")]<-"PIGW"
  names(biopoint)[names(biopoint) %in% c("Lygodium")]<-"OWCF"
  names(biopoint)[names(biopoint) %in% c("Sapium.sebiferum")]<-"CHTA"
  names(biopoint)[names(biopoint) %in% c("Begonia_Cucullata")]<-"BECU"
  ####
  
  
  # I don't know if this will be useful to you or not, but here I calculate the average distance that the points are spaced out in
    # the grid in order to use them for mapping.  It doesn't hurt to leave this, but if point spacing isn't important for you, you 
      # can delete this section
  ##### Calculate the distance of the point samples, in order to get the size of the species points
  biopoin<-biopoint
  coordinates(biopoin)<-cbind(biopoin$X, biopoin$Y)
  
  suppressWarnings(
  proj4string(biopoin)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  suppressWarnings(
  biopoin<-spTransform(biopoin, CRS("+init=EPSG:3857")))
  
  pointdis<-mean(knn.dist(biopoin@coords, k=1))
  ###
  biopoint$pointdis<-rep(pointdis, nrow(biopoint))
  #####
  
  
  
  biopoint  # for this function, you need to call this at the end
})

# this is what you get after the above function. You'll notice that the point distance column is in the middle, because the first data 
  #frame did not have all those other columns.  We'll have to fill in the 0's again (you will see that a few lines down)...I should have just done it here and not in the function, but
    # it can be helpful to see what columns existed for which lakes when you're checking your data
head(point)   


## here i read in the key from FWRI (this is the one that is always on the other excel sheet with the data)
key<-read.csv("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/rtables/shinyapps/whoml/www/ltmplants/biobasespecieskey_2019.csv")


###mash the species names together for each point 
##Make a point layer that has all species found
speciess<-names(point)[names(point) %in% key$P_code] #all the species found in the key

#here are all the ones that didn't make it, we'll have to go back and insert those in the key unless they are "Lake","Date","Crew","Site","Y","X", etc
names(point)[!names(point) %in% key$P_code]  


pointall<-point[c("Lake","Date","X","Y","pointdis",speciess)]  ## make a new data frame that skips anything that we didn't key out # YOU MAY NOT WANT TO DO THIS STEP
0->pointall[is.na(pointall)]  # put in the other 0's

# here I actually drop 0's out of my dataset and then make two new columns that have all of the names combined of everything that was found at that spot
  # I THINK YOU CAN SKIP THESE THREE LINES. YOU PROBABLY WANT TO PRESERVE 0'S FOR THE ANALYSIS
    # if you are going to run this, make sure that ONLY your first 5 columns are non-species columns because this is looking at column 6 through the end
pointall<-pointall[rowSums(pointall[6:length(pointall)]==0) !=length(pointall[6:length(pointall)]) ,] #drop rows that are 0 for all species
pointall$spp<-apply(pointall[6:length(pointall)], 1, function(x) paste0(key$scientific_name[match(names(which(!x==0)),key$P_code)],collapse=", ") )
pointall$sppcommon<-apply(pointall[6:(length(pointall)-1)], 1, function(x) paste0(key$common_name[match(names(which(!x==0)),key$P_code)],collapse=", ") )


## here I melt the data frame into long format, keeping the first 5 columns and having two columns beside that, one for the species name and another for its density value
pointallmelt<-reshape2::melt(pointall, id=names(pointall)[1:5])
pointallmelt<-pointallmelt[pointallmelt$value>0,]

# then, finally i add the scientific and common names to the data set
pointallmelt$scientific_name<-key$scientific_name[match(pointallmelt$variable, key$P_code)]
pointallmelt$common_name<-key$common_name[match(pointallmelt$variable, key$P_code)]

head(pointallmelt)

# write.csv(pointallmelt, paste0(pointfol,"/pointall.csv"), row.names=F)
#write.fst(pointallmelt, paste0(pointfol,"/pointall.fst"))



#############
## from here you can probably add a year column to each of these and then combine them
  # as you have already experienced, the date column is not very reliable because they only fill that in if there wasn't all 0's found at the point
    # if they find all 0's then they just leave the whole thing blank. ...so you can't just run year(pointallmelt$Date) to get the year





############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
## fish data! 


# Just to recap on the fish data
# This is the same data that was used to create that summary you got from Hoyer.  You should proabably ask him
# for the FWRI LTM Fish sampling protocol manual.  If he doesn't have it, I can get it for you.
# 
# The data is mainly in two tables that are summary tables from the database.  One is the "fish data" 
# (transectdata db1 is what i named it) which is pretty much a row for every fish sampled (however, if 
# they catch a bunch of 1 inch fish, they sometimes batch weigh them, so there is a Count column that 
# shows how many fish were associated with that data. This column is usually 1, except for those times).
# 
# The other main table is the "event" data (I named it transectdata db2), and that is basically the data
# associated with driving that transect.  
# 
# The weird thing about these tables is that they restructed their database and had to create a new master id.  The id you should 
# use is the SchemaMasterID. 
# 
# I think for these analyses, only "Electrofishing" gear and "Standard" sampling type should be used. They do a lot of "hobby shocking"
# and that isn't part of the standard protocol.  I believe spring shocking (which only goes after bass) is 900 seconds and fall shocking
# (which is community shocking) is 600 seconds.  We may want to filter for only those efforts too.  You can see what the target was (LMB for 
# just bass or ALL for community) in the Target.Species column. 
# 
# If you want to do any spatial analysis with this data, don't use the beginlat, beginlong data, because it is riddled with errors. I have
# their transect plans for most of these, and those are more reliable. 
# 
# In additon to these two tables, I have a plant table that gives the domininant species per "event". You'll see the %submersed that is 
# the same value that is already in the 'event' summary table.
# Another summary table is the one that summarizes number of individual bass sampled at each transect 'event'. This is broken into the 
# different size classes.  This data is what is useful for examining if there is any relationship with having higher substock one year and 
# that leading to more fish in the larger size classes in 2 to 3 years (the 'strong year class' argument). 
# 
# I think since we have the data, it is good to go ahead and merge in LakeWatch data and confirm what many other studies have found
# in the past that more Phosphorus in a system leads to more and larger fish. 



# Reading Maceina's review of Hoyer and Canfield's 1996 paper, we may want to take everything we are analyzing in these data and 
# split them into 'large' and 'small' waterbodies.  ...something to think about after doing the analyses with all data




library(fst)
library(data.table)
library(plotly)

# AK note: I believe this is Allfish.csv, but only LMB data
# fish<-read.fst("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/R Scripts/transectdata1_db.fst") #the "fish" data
fish <- read.fst("original-data/dew-fish-analysis/transectdata1_db.fst")
head(fish)

fishh<-fish[fish$Species == "LMB",] # only keep rows that are Largemouth bass
fishh<-fishh[fishh$Type == "Standard",]  #only standard sampling
fishh<-fishh[fishh$Gear=="Electrofishing",]  #only sampling with Electrofishing
fishh$Total.Length<-as.numeric(fishh$Total.Length)  # turn these columns numeric
fishh$Total.Weight<-as.numeric(fishh$Total.Weight)
fishh$Total.Weight<-fishh$Total.Weight / as.numeric(fishh$Count)  ## turn this into a weight per fish

fishh$EffortSeconds<-as.numeric(fishh$EffortSeconds)

## remove rows where they only got a length of the fish and not get the weight
fishh<-fishh[!is.na(fishh$Total.Weight),]

## remove large typos (get rid of fish weight records over 20lbs)
fishh<-fishh[fishh$Total.Weight < 9000,]

# take out the entry where they didn't put in an effort amount
fishh<-fishh[which(fishh$EffortSeconds>0),] 

#############
# AK: this is metadata/environmental information for electrofishing event
## merge in the "event" (transect) data to the electrofishing dataset
# event<-read.fst("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/R Scripts/transectdata2_db.fst")
event <- read.fst("original-data/dew-fish-analysis/transectdata2_db.fst")
head(event)

## take avg submersed plant cover by lake and year 
event$X..Submersed<-as.numeric(event$X..Submersed) #there will be some 'NAs by coersion' because the SQL database uses "NULL", not "NA"
event$X..Submersed[is.na(event$X..Submersed)]<-0  #if submersed is NA, put 0
event$X..Total.Coverage<-as.numeric(event$X..Total.Coverage)
event$X..Total.Coverage[is.na(event$X..Total.Coverage)]<-0

# Make a summary dataframe, creating an average SAV percent by waterbody and year
avgsub<-as.data.frame(
  data.table(event)[, list(avgsubmersed=mean(X..Submersed, na.rm=T)), 
                                        by=list(WaterBodyID, Year)]
  )


###############################################
# #merge average submersed in here and save these for later
event<-merge(event, avgsub, by=c("WaterBodyID", "Year"), all.x=T) #merge the average back onto the main 'event' dataset
##

####################################################


#################################################################################################
#################################################################################################
## Look at avg biomass per lake vs mean % submersed 

# first take total mass of fish per transect , then average that by lake
  #these are two summary tables, the first one goes through all the fish rows per transect ('event') and the second one does it by
    # waterbody id and year
transectsum<-data.table(fishh)[, list( TotalWeight = sum(Total.Weight, na.rm = T) ),
                               
                               by=list(WaterBodyID, Water.Body, County, Year,SchemaMasterID, EffortSeconds)]

transectsum$TotalWeightPerEffort<-transectsum$TotalWeight/transectsum$EffortSeconds
fishyearlakesum<-as.data.frame(transectsum[, list(AvgBiomassPerTransect=mean(TotalWeight, na.rm=T), AvgWeightPerUnitEffort=mean(TotalWeightPerEffort, na.rm=T)), 
                                           by=list(WaterBodyID, Water.Body, County, Year)])


##
#merge in avg submersed per lake and year that we summarized already in the above section
yearlake<-merge(fishyearlakesum, avgsub, by=c("WaterBodyID", "Year"), all.x=T)



### plot lakewide avg submersed vs avg transect biomass

avgbiomassmod<-lm(AvgWeightPerUnitEffort~avgsubmersed, data=yearlake)  # check out a simple linear model
summary(avgbiomassmod)

plot_ly() %>% 
  add_trace(data = yearlake, x = yearlake$avgsubmersed, y = yearlake$AvgWeightPerUnitEffort, type="scatter", name="TotalWeightvsSubmersed", mode="markers",  # mode="lines+markers",
            text=paste0(yearlake$Water.Body, " - Lakewide Avg Subm. Cover: ", round(yearlake$avgsubmersed)),
            marker = list(size = 10,
                          color = 'rgba(193, 182, 255, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
  add_trace(x=yearlake$avgsubmersed, y=as.numeric(predict(avgbiomassmod)), name="Trend",mode='lines') %>% 
  
  layout(title = 'Lakewide Mean Submersed Veg(%) vs Lakewide Mean Biomass Per Unit Effort', showlegend=FALSE,
         yaxis = list(title="Weight (g)"),
         xaxis = list(title="Total Submersed"))




######################
# what about lakes that maintain more SAV over time?  Does looking at it that way instead of single years
  # yield a different result?

################
## take out single years and do a whole lake average
avgbiomassperlake<-as.data.frame(data.table(yearlake)[, list(AvgBiomassPerTransect=mean(AvgBiomassPerTransect, na.rm=T),AvgWeightPerUnitEffort=mean(AvgWeightPerUnitEffort, na.rm = T) ,avgsubmersed=mean(avgsubmersed, na.rm=T)), 
                                                      by=list(WaterBodyID, Water.Body, County)]) # took out year here

avgbiomassmod<-lm(AvgWeightPerUnitEffort~avgsubmersed, data=avgbiomassperlake)
summary(avgbiomassmod)
plot_ly() %>% 
  add_trace(data = avgbiomassperlake, x = avgbiomassperlake$avgsubmersed, y = avgbiomassperlake$AvgWeightPerUnitEffort, type="scatter", name="TotalWeightvsSubmersed", mode="markers",  # mode="lines+markers",
            text=paste0(avgbiomassperlake$Water.Body, " - Lakewide Avg Subm. Cover: ", round(avgbiomassperlake$avgsubmersed)),
            marker = list(size = 10,
                          color = 'rgba(193, 182, 255, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
  add_trace(x=avgbiomassperlake$avgsubmersed, y=as.numeric(predict(avgbiomassmod)), name="Trend",mode='lines') %>% 
  
  layout(title = 'Lakewide Mean SAV (%) vs Lakewide Mean LMB Biomass per Unit Effort (seconds)', showlegend=FALSE,
         yaxis = list(title="Weight (g)"),
         xaxis = list(title="Mean Total Submersed"))





###########################################################################################################################
## what about the maximum weight bass that has been sampled per waterbody over the study period? Do lakes that maintain higher
  # SAV % yield larger maximum weight bass?

## summarize by the maximum weight per waterbody and extract the event ID for where that fish was sampled
fishsum<-data.table(fishh)[, list(MaxWeight = max(Total.Weight, na.rm=T)  , YearofMaxWeight = Year[which(Total.Weight == max(Total.Weight, na.rm=T))]  ,  EventIDofMaxWeight= SchemaMasterID[which(Total.Weight == max(Total.Weight, na.rm=T))]   ), 
                           by=list(WaterBodyID, Water.Body, County)]

head(fishsum)
nrow(fishsum) # 230 waterbodies


## filter out columns in the event data
eventt<-event[c("SchemaMasterID","X..Total.Coverage", "X..Submersed", "Bottom.Type.1", "avgsubmersed")]
head(eventt)


fishmer<-merge(fishsum, eventt,by.x="EventIDofMaxWeight" ,by.y="SchemaMasterID", all.x=T)
head(fishmer)

summary(fishmer)


head(avgbiomassperlake)
# merge in record bass
fishmerrec<-merge(avgbiomassperlake, fishmer[,c("WaterBodyID", "MaxWeight")], by=c("WaterBodyID"), all.x=T)
head(fishmerrec)



avgmaxmod<-lm(MaxWeight~avgsubmersed, data=fishmerrec)  #simple linear model
summary(avgmaxmod)
plot_ly() %>% 
  add_trace(data = fishmerrec, x = fishmerrec$avgsubmersed, y = fishmerrec$MaxWeight, type="scatter", name="TotalWeightvsSubmersed", mode="markers",  # mode="lines+markers",
            text=paste0(fishmerrec$Water.Body, " - Lakewide Avg Subm. Cover: ", round(fishmerrec$avgsubmersed)),
            marker = list(size = 10,
                          color = 'rgba(193, 182, 255, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
  #add_trace(x=fishmerrec$avgsubmersed, y=as.numeric(predict(avgmaxmod)), name="Trend",mode='lines') %>% 
  
  layout(title = 'Lakewide Mean SAV (%) vs Lakewide Max Weight (all years)', showlegend=FALSE,
         yaxis = list(title="Weight (g)"),
         xaxis = list(title="Total Submersed"))




#############################################################################################



#################################################################################################
#################################################################################################
#################################################################################################
# lw<-read.csv("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/R Scripts/PvsBass/LakewatchData.csv")
lw <- read.csv("original-data/dew-fish-analysis/LakewatchData.csv")
head(lw)


lwsummary<-data.table(lw)[, list(tenyearavg = mean(AnnualMeanTP[Year %in% 2010:2020], na.rm=T),   fiveyearavg = mean(AnnualMeanTP[Year %in% 2015:2020])), 
                          by=list(GNIS_ID, Lake, County, Latitude, Longitude)]


## Let's merge the GNIS id and coordinates onto this dataset using the FWC_id as a common key
## this is a key from the data you gave me a while back.  You may have a more up to date list now
# gkey<-read.csv("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/R Scripts/PvsBass/GNIS_key.csv") 
gkey <- read.csv("original-data/dew-fish-analysis/GNIS_key.csv")

head(gkey)
gkey$GnisID<-as.numeric(gkey$GnisID)


## merge in the biggest bass to this dataset
fishmerk<-merge(fishmer, gkey, by.x="WaterBodyID", by.y="WaterBodyId", all.x=T)
head(fishmerk)
summary(fishmerk)



#lwfish<-merge(fishmerk, lwsummary, by.x="GnisID", by.y="GNIS_ID", all.x=T)

fishmerk$TPtenyearavg<-lwsummary$tenyearavg[match(fishmerk$GnisID, lwsummary$GNIS_ID)]
fishmerk$TPfiveyearavg<-lwsummary$fiveyearavg[match(fishmerk$GnisID, lwsummary$GNIS_ID)]

head(fishmerk)


#plot(fishmerk$TPtenyearavg,fishmerk$MaxWeight)
fishmerk<-fishmerk[!is.na(fishmerk$TPtenyearavg),]
fishmerkmod<-lm(MaxWeight~TPtenyearavg, data=fishmerk)  ## simple linear model of maximum weight bass vs mean P levels
summary(fishmerkmod)

plot_ly() %>% 
  add_trace(data = fishmerk, x = fishmerk$TPtenyearavg, y = fishmerk$MaxWeight, type="scatter", 
            text=paste0(fishmerk$Water.Body, " - P: ", round(fishmerk$TPtenyearavg)),
            marker = list(size = 10,
                          color = 'rgba(255, 182, 193, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
  add_lines(x=fishmerk$TPtenyearavg, y=predict(fishmerkmod), name="Trend") %>% 
  
  layout(title = 'Ten Year Mean TP vs Max Weight', showlegend=F,
         yaxis = list(title="Weight (g)"),
         xaxis = list(title="Total P"))


 ### WE PROBABLY WANT TO ALSO DO THIS ^ FOR BIOMASS PER UNIT EFFORT



##################################################################################################


# pl<-read.csv("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/R Scripts/plant_type.csv")
pl <- read.csv("original-data/dew-fish-analysis/plant_type.csv")
head(pl)

## merge the dominant sav to the transect biomass sum and compare hydrilla to other

head(transectsum)

transectsum$dominantSAV<-pl$Sub1Scientific[match(transectsum$SchemaMasterID, pl$SchemaMasterID)]
transectsum$percentsubmersed<-as.numeric(pl$PercentSubmersed[match(transectsum$SchemaMasterID, pl$SchemaMasterID)])
transectsumy<-transectsum[!transectsum$dominantSAV=="NULL",]  #get rid of transects that had no SAV
transectsumy$hy<-ifelse(transectsumy$dominantSAV=="Hydrilla verticillata", "Hydrilla","Other")



plot_ly(transectsumy, y=transectsumy$TotalWeightPerEffort, color=transectsumy$hy, type="box") %>% 
  layout(title = 'Total LMB Biomass by Dominant SAV Species', showlegend=F, yaxis = list(title="Biomass Per Unit Effort per Transect (g/sec)"))

t.test(transectsumy$TotalWeightPerEffort[transectsumy$hy=='Hydrilla'], transectsumy$TotalWeightPerEffort[transectsumy$hy=='Other'])
0.5274110 / 0.5730358  #hydrilla supports 92% the biomass compared to "other"

## this is very rough - we probably need to make sure all of the plants they put as 'dominant SAV' are actually SAV, and filter out
#the ones that are not, then we should filter out other exotics and run an exotic vs native, as well as hydrilla vs native
  ## we may want to just check if this result is just a function of biovolume. Hydrilla grows thicker than native SAV, so it 
    ## would be important to see if that is the reason that it supports less fish and/or if there are other reasons
 
## we may want to also look into if some plants support more bass than others (can do lake bottom types too)
plot_ly(transectsumy, y=transectsumy$TotalWeight, color=transectsumy$dominantSAV, type="box")



#############################################################################################################
#############################################################################################################
## Mean harvestable bass size vs SAV

# defining >= 250mm as adult like in hoyer and canfield

fish250<-fishh[fishh$Total.Length>=250,]
fish250<-fish250[!is.na(fish250$Total.Length),]
transectsum250<-data.table(fish250)[, list( meanWeight = mean(Total.Weight, na.rm=T) ), 
                                    by=list(WaterBodyID, Water.Body, County, Year,SchemaMasterID)]



transectsum250$dominantSAV<-pl$Sub1Scientific[match(transectsum250$SchemaMasterID, pl$SchemaMasterID)]
transectsum250$percentsubmersed<-as.numeric(pl$PercentSubmersed[match(transectsum250$SchemaMasterID, pl$SchemaMasterID)])
transectsum250$percentsubmersed[is.na(transectsum250$percentsubmersed)]<-0

head(transectsum250)


fishlakeyearsum250<-data.table(transectsum250)[, list( meanWeight = mean(meanWeight, na.rm=T), meanSAV=mean(percentsubmersed, na.rm=T) ), 
                                               by=list(WaterBodyID, Water.Body, County, Year)]


#### Mean weight of adult bass vs % SAV
meanweightmod<-lm(meanWeight~meanSAV, data=fishlakeyearsum250)  #simple linear model
summary(meanweightmod)
plot_ly() %>% 
  add_trace(x = fishlakeyearsum250$meanSAV, y = fishlakeyearsum250$meanWeight, type="scatter", name="TotalWeightvsSubmersed", mode="markers",  # mode="lines+markers",
            text=paste0(fishlakeyearsum250$Water.Body, " - Lakewide Avg Subm. Cover: ", round(fishlakeyearsum250$meanSAV)),
            marker = list(size = 10,
                          color = 'rgba(193, 182, 255, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
   add_trace(x=fishlakeyearsum250$meanSAV, y=as.numeric(predict(meanweightmod)), name="Trend",mode='lines') %>% 
  
  layout(title = 'Percent SAV vs Mean Adult Bass (>=250mm TL) Weight', showlegend=FALSE,
         yaxis = list(title="Mean Adult Weight (g)"),
         xaxis = list(title="Total Submersed"))



#############################################################################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################
######################################################################################################
## lmb size categories


# This is the summary table that shows number of individual bass in each size class for each transect event. Using this we 
# can investigate if having higher numbers of substock in year n-2 and n-3 leads to an increased number of fish in the larger
# size classes 
## we can summarize the same thing from the 'fish' and 'event' datasets, so if any important parameters are missing, we can do
## that instead

# cat<-read.csv("C:/users/alex.dew/downloads/lmb_sizecategories.csv")
cat <- read.csv("original-data/dew-fish-analysis/lmb_sizecategories.csv")
head(cat)
unique(cat$EffortMinutes)  #  HERE THEY SUMMARIZED EFFORT IN MINUTES, NOT SECONDS - should be 15 for standard bass sampling
cat$EffortMinutes<-as.numeric(cat$EffortMinutes)

## summary table that calculates fish per minute per waterbody per year per size class
catsum<-as.data.frame(data.table(cat)[, list(FishPerMin=sum(Total, na.rm = T)/sum(EffortMinutes, na.rm = T)), 
                                      by=list(WaterBody, WaterBodyID, Year, StockSizeGroup)])

## going from long format to wide, so that we can have a different column for each size class
catd<-reshape2::dcast(data=catsum, WaterBody+WaterBodyID+Year~StockSizeGroup, sum,value.var="FishPerMin")
head(catd)


##################################
## look at substock in year n-2 and n-3 compared to other size classes 

## create new columns to put the value for substock n-2 and n-3
catd$nminus2<-NA
catd$nminus3<-NA

## loop through and put the appropriate value in the columns
for (j in unique(catd$WaterBodyID)){
  
  for (i in min(catd$Year[catd$WaterBodyID==j]):max(catd$Year[catd$WaterBodyID==j])){
    
    subwb<-catd[catd$WaterBodyID == j,]
    
    if(length(subwb$Substock[which(subwb$Year == i-2)])==0){
      catd$nminus2[catd$Year==i & catd$WaterBodyID==j]<-NA
    } else {
      catd$nminus2[catd$Year==i & catd$WaterBodyID==j]<-subwb$Substock[which(subwb$Year == i-2)]
    }
    
    if(length(subwb$Substock[which(subwb$Year == i-3)])==0){
      catd$nminus3[catd$Year==i & catd$WaterBodyID==j]<-NA
    } else {
      catd$nminus3[catd$Year==i & catd$WaterBodyID==j]<-subwb$Substock[which(subwb$Year == i-3)]
    }
    
  }
}


head(catd)


#####
#n - 2 substock vs trophy
catdn2<-catd[!is.na(catd$nminus2),]  # the first two years in the data will be NA, so remove those

## trophy+ (24" and above) vs substock year n-2
modn2<-lm(data=catdn2, `>Trophy`~nminus2)
summary(modn2)
plot_ly() %>% 
  add_trace(x=catdn2$nminus2, y=catdn2$`>Trophy`, type="scatter", text=catdn2$WaterBody) %>% 
  #add_trace(x=catdn2$nminus2, y=as.numeric(predict(modn2)), name="Trend",mode='lines') %>% 
  layout(title = 'Substock year n-2 vs Trophy+ Bass',
         yaxis = list(title="Trophy CPUE"),
         xaxis = list(title="Substock CPUE year n-2"))


#n - 3 substock vs trophy
catdn3<-catd[!is.na(catd$nminus3),]

modn3<-lm(data=catdn3, `>Trophy`~nminus3)
summary(modn3)
plot_ly() %>% 
  add_trace(x=catdn3$nminus3, y=catdn3$`>Trophy`, type="scatter", text=catdn3$WaterBody) %>% 
  #add_trace(x=catdn3$nminus3, y=as.numeric(predict(modn3)), name="Trend",mode='lines') %>% 
  layout(title = 'Substock year n-3 vs Trophy+ Bass',
         yaxis = list(title="Trophy CPUE"),
         xaxis = list(title="Substock CPUE year n-3"))


###  ###
#n - 2 substock vs  mem-trophy
catdn2<-catd[!is.na(catd$nminus2),]

modn2<-lm(data=catdn2, `Mem-Trophy`~nminus2)
summary(modn2)
plot_ly() %>% 
  add_trace(x=catdn2$nminus2, y=catdn2$`Mem-Trophy`, type="scatter", text=catdn2$WaterBody) %>% 
  #add_trace(x=catdn2$nminus2, y=as.numeric(predict(modn2)), name="Trend",mode='lines') %>% 
  layout(title = 'Substock year n-2 vs Memorable-Trophy Bass',
         yaxis = list(title="Memorable-Trophy Bass CPUE"),
         xaxis = list(title="Substock CPUE year n-2"))


#n - 3 substock vs Mem-Trophy
catdn3<-catd[!is.na(catd$nminus3),]
#catdn3<-catdn3[catdn3$`Mem-Trophy`>0,]
modn3<-lm(data=catdn3, `Mem-Trophy`~nminus3)
summary(modn3)
plot_ly() %>% 
  add_trace(x=catdn3$nminus3, y=catdn3$`Mem-Trophy`, type="scatter", text=catdn3$WaterBody) %>% 
  #add_trace(x=catdn3$nminus3, y=as.numeric(predict(modn3)), name="Trend",mode='lines') %>% 
  layout(title = 'Substock year n-3 vs Memorable-Trophy Bass',
         yaxis = list(title="Memorable-Trophy Bass CPUE"),
         xaxis = list(title="Substock CPUE year n-3"))



###  ###
#n - 2 substock vs legal size (all fish >= 16")
catdn2<-catd[!is.na(catd$nminus2),]

modn2<-lm(data=catdn2, Legal~nminus2)
summary(modn2)
plot_ly() %>%
  add_trace(x=catdn2$nminus2, y=catdn2$Legal, type="scatter", text=catdn2$WaterBody) %>%
  #add_trace(x=catdn2$nminus2, y=as.numeric(predict(modn2)), name="Trend",mode='lines') %>%
  layout(title = 'Substock year n-2 vs >=16" Bass',
         yaxis = list(title=">=16 in. Bass CPUE"),
         xaxis = list(title="Substock CPUE year n-2"))


#n - 3 substock vs trophy and pref-mem
catdn3<-catd[!is.na(catd$nminus3),]

modn3<-lm(data=catdn3, Legal~nminus3)
summary(modn3)
plot_ly() %>%
  add_trace(x=catdn3$nminus3, y=catdn3$Legal, type="scatter", text=catdn3$WaterBody) %>%
  #add_trace(x=catdn3$nminus3, y=as.numeric(predict(modn3)), name="Trend",mode='lines') %>%
  layout(title = 'Substock year n-3 vs >=16" Bass',
         yaxis = list(title=">=16 in. Bass CPUE"),
         xaxis = list(title="Substock CPUE year n-3"))




##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
## we can also use this data to see if the number of individuals in each size class change with submersed veg %

# size classes vs % submersed


subsum<-as.data.frame(data.table(event)[, list( meanSub=mean(X..Submersed, na.rm=T)  ), 
                                        by=list( WaterBodyID, Year)])


catd<-merge(catd, subsum, by=c("WaterBodyID", "Year"), all.x=T)
head(catd)

modsub<-lm(data=catd, `>Trophy`~meanSub)
summary(modsub)
plot_ly() %>% 
  add_trace(x=catd$meanSub, y=catd$`>Trophy`, type="scatter", text=catd$WaterBody) %>% 
  add_trace(x=catd$meanSub, y=as.numeric(predict(modsub)), name="Trend",mode='lines') %>% 
  layout(title = 'Mean Submersed % vs Number of Trophy+ (24"+) Bass',
         yaxis = list(title="Trophy CPUE"),
         xaxis = list(title="Mean Submersed %"))



modsub<-lm(data=catd, `Mem-Trophy`~meanSub)
summary(modsub)
plot_ly() %>% 
  add_trace(x=catd$meanSub, y=catd$`Mem-Trophy`, type="scatter", text=catd$WaterBody) %>% 
  add_trace(x=catd$meanSub, y=as.numeric(predict(modsub)),  name="Trend",mode='lines') %>% 
  layout(title = 'Mean Submersed % vs Number of Mem-Trophy Bass',
         yaxis = list(title="Mem-Trophy CPUE"),
         xaxis = list(title="Mean Submersed %"))



modsub<-lm(data=catd, Legal~meanSub)
summary(modsub)
plot_ly() %>% 
  add_trace(x=catd$meanSub, y=catd$Legal, type="scatter", text=catd$WaterBody) %>% 
  add_trace(x=catd$meanSub, y=as.numeric(predict(modsub)), name="Trend",mode='lines') %>% 
  layout(title = 'Mean Submersed % vs Number of Legal Bass',
         yaxis = list(title="Legal CPUE"),
         xaxis = list(title="Mean Submersed %"))


## PROBABLY SHOULD RUN FOR ALL CLASSES. WOULD EXPECT SUBSTOCK TO INCREASE WITH SAV % WHILE BIGGER FISH DECREASE




#############################################################################################################
#############################################################################################################
#############################################################################################################
# Species richness per transect vs dominant SAV on transect

## here is the dataset of all the data where 'ALL' fish were targeted (community sampling)


#allf<-read.csv("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/R Scripts/PvsBass/allfish.csv")
allf <- read.csv("original-data/dew-fish-analysis/allfish.csv")

head(allf)

allf<-allf[allf$EffortSeconds==600,]  #only taking effort = 600 seconds


## get plant data summarized 
uniqtransects<-data.table(allf)[,  .(verifyall1 = length(unique(SchemaMasterID))) , 
                                by=list(WaterBodyID, WaterBody, County, Year,SchemaMasterID)]
uniqtransects$dominantSAV<-pl$Sub1Scientific[match(uniqtransects$SchemaMasterID, pl$SchemaMasterID)]  #put in the dominant SAV
uniqtransects$percentsubmersed<-as.numeric(pl$PercentSubmersed[match(uniqtransects$SchemaMasterID, pl$SchemaMasterID)]) # put in the percent SAV
uniqtransects$percentsubmersed[is.na(uniqtransects$percentsubmersed)]<-0  #put 0 for NA SAV

## now make a summary of the percent submersed veg by waterbody and year
sumpercsub<-data.table(uniqtransects)[, list( percentsubmersed=mean(percentsubmersed, na.rm=T) ), 
                                      by=list(WaterBodyID,  Year)]


## summary dataframe of species count by lake and year
lakesumall<-data.table(allf)[,  .(speciesCount = length(unique(SpeciesScientific))) , 
                             by=list(WaterBodyID, WaterBody, County, Year)]
#merge in mean submersed % by lake and year
fishlakeyearsumall<-merge(lakesumall, sumpercsub, by=c("WaterBodyID", "Year"), all.x=T)



head(fishlakeyearsumall)


#### Mean species richness vs % SAV
speciesmod<-lm(speciesCount~percentsubmersed, data=fishlakeyearsumall)  ## this is where you mentioned us looking into a GAM
summary(speciesmod)
plot_ly() %>% 
  add_trace(x = fishlakeyearsumall$percentsubmersed, y = fishlakeyearsumall$speciesCount, type="scatter", name="TotalWeightvsSubmersed", mode="markers",  # mode="lines+markers",
            text=paste0(fishlakeyearsumall$WaterBody, " - Transect Subm. Cover: ", round(fishlakeyearsumall$percentsubmersed)),
            marker = list(size = 10,
                          color = 'rgba(193, 182, 255, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
  add_trace(x=fishlakeyearsumall$percentsubmersed, y=as.numeric(predict(speciesmod)), name="Trend",mode='lines') %>% 
  
  layout(title = 'Percent SAV vs Species Richness by Lake and Year', showlegend=FALSE,
         yaxis = list(title="Number of Species"),
         xaxis = list(title="Total Submersed"))

############################
## At least how much SAV do we need to get similar species counts?
#### Mean species richness vs % SAV
fishlakeyearsumall<-fishlakeyearsumall[fishlakeyearsumall$percentsubmersed>=2,]  # what if we filter out everything under 2% SAV
speciesmod<-lm(speciesCount~percentsubmersed, data=fishlakeyearsumall)
summary(speciesmod)
plot_ly() %>% 
  add_trace(x = fishlakeyearsumall$percentsubmersed, y = fishlakeyearsumall$speciesCount, type="scatter", name="TotalWeightvsSubmersed", mode="markers",  # mode="lines+markers",
            text=paste0(fishlakeyearsumall$WaterBody, " - Transect Subm. Cover: ", round(fishlakeyearsumall$percentsubmersed)),
            marker = list(size = 10,
                          color = 'rgba(193, 182, 255, .9)',
                          line = list(color = 'rgba(152, 0, 0, .8)',
                                      width = 2))) %>% 
  
  add_trace(x=fishlakeyearsumall$percentsubmersed, y=as.numeric(predict(speciesmod)), name="Trend",mode='lines') %>% 
  
  layout(title = 'Percent SAV vs Species Richness', showlegend=FALSE,
         yaxis = list(title="Number of Species"),
         xaxis = list(title="Total Submersed"))
























