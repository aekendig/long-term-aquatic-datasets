
# LTM Plant frequencies



# find all the point files
dir<-list.dirs("S:/avm")

dirs<-dir[dir %like% "intercept_data|intercept-data"] ## take only directories with these two naming conventions
dirs
dirs<-dirs[dirs %like% "2015|2016|2017|2018|2019|2020|2021"] ## match the year folders to get rid of the outer directories
#dirs<-dirs[!dirs %like% "EmeraldaMarshArea3"]


# cycle through directories to get file names
fils<-NA
for (i in 1:length(dirs)){
  fils<-c(fils,list.files(dirs[i], pattern="lant", full.names=T))
  
}

pointfiles<-fils[!is.na(fils)]

pointfiles<-fils[!fils %like% "PrimaVista"]

############################################################################
library(readxl)

#pointdata<-pointfiles[1] #test


point<-ldply(pointfiles ,function(pointdata){

  
#poin<-read.csv(pointfiles[1])

  poin<-read_excel(pointdata, sheet=1) %>% as.data.frame()


colnames(poin)[tolower(names(poin)) %in% c("longitude","long","lng","x")]<-"X"
colnames(poin)[tolower(names(poin)) %in% c("latitude","lat","y")]<-"Y"

print(head(poin, 1))


NoAnames<-names(poin[ grepl("No_Access|NA_", names(poin)) ])  # names of columns that will represent a non-accessible site. 
totalSites<-nrow( poin[rowSums(poin[NoAnames] > 0, na.rm = T) == 0 ,] ) # number of sites that were actually accessible


plantCols<-names(poin)[! tolower(names(poin)) %in% c("lake", "date", "crew", "site", "x", "y", tolower(NoAnames)) ]
poin[plantCols][is.na(poin[plantCols])]<-0  # fill any NAs with 0's

plantfreq<-apply(poin[plantCols], 2 , function(x) length(which(x>0)))
freq3<-apply(poin[plantCols], 2 , function(x) length(which(x==3)))
freq2<-apply(poin[plantCols], 2 , function(x) length(which(x==2)))
freq1<-apply(poin[plantCols], 2 , function(x) length(which(x==1)))
pf<-data.frame(present=plantfreq, frequency=plantfreq/totalSites, 
               present3=freq3, frequency3=freq3/totalSites, 
               present2=freq2, frequency2=freq2/totalSites, 
               present1=freq1, frequency1=freq1/totalSites ,
               PCode=names(plantfreq))
pf

basenam<-gsub( "plantdata","",tolower( gsub(".*_", "", basename(pointdata)))) ## sub out everything after "_" and then sub out "PlantData"
year<-gsub("\\..*","",basenam)    # take away the file extension and we just have the year left

pf$year<-year  ## get this year from the folder, date column may be blank
pf$Lake<-poin$Lake[1]  #might want to get this from the main folder too...

pf

})




## key out the point with names, types, and origins
key<-read.csv("C:/Users/alex.dew/OneDrive - Florida Fish and Wildlife Conservation/rtables/shinyapps/whoml/www/ltmplants/biobasespecieskey_2019.csv")
key<-read_excel(pointfiles[2], sheet=2) %>% as.data.frame()


head(point)
nrow(point)
head(key)

point$scientific_name<-key$SCIENTIFIC_NAME[match(point$PCode, key$P_CODE)]
point$common_name<-key$COMMON_NAME[match(point$PCode, key$P_CODE)]

point$scientific_name[point$scientific_name=="NULL"]<-NA
point$scientific_name<-ifelse(is.na(point$scientific_name), point$common_name, point$scientific_name) # if no sci name, give it the common name in the sciname column

point$scientific_name<-ifelse(is.na(point$scientific_name), point$PCode, point$scientific_name)   #if still NA, give it the PCode value

point$origin<-key$ECO_TYPE[match(point$PCode, key$P_CODE)]
point$planttype<-key$TYPE[match(point$PCode, key$P_CODE)]

head(point)


point[which(is.na(point$scientific_name)),]
point[which(point$scientific_name=="NULL") ,]
point[which(point$origin=="NULL") ,]

## let's take out tussocks (floating islands) and filamentous algae
pointy<-point[which(!point$origin=="NULL") ,]


head(pointy)
nrow(pointy)

unique(pointy$Lake)
pointy$Lake[pointy$Lake=="Conway (North Lobe)"]<-'North Conway'
pointy$Lake[pointy$Lake=="NorthConway"]<-'North Conway'
pointy$Lake[pointy$Lake=="Conway (South Lobe)"]<-'South Conway'
pointy$Lake[pointy$Lake=="SouthConway"]<-'South Conway'

pointy$Lake[pointy$Lake %in% c("East Toho", "EastToho")]<-"East Tohopekaliga"

pointy$Lake[pointy$Lake=="Jackson" & pointy$year==2016]<-'Jackson (Highlands)'

pointy$Lake[pointy$Lake=="Jackson(Leon)" ]<-"Jackson (Leon)" 

pointy$Lake[pointy$Lake=="Rodman" ]<-"Rodman Reservoir"

pointy$Lake[pointy$Lake=="Jackson"][1:4]<-'Jackson (Highlands)'
pointy$Lake[pointy$Lake=="Jackson"]<-"Jackson (Leon)"

head(pointy)

write.csv(pointy, file="C:/users/alex.dew/downloads/LTM_Plant_Freqs.csv", row.names=F)
