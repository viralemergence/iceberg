###################################

#IUCN grabber: run back half of iucn habitat
#####################

# 0. Land use types

lutypes <- c("primf","primn","secdf","secdn","urban","c3ann","c4ann","c3per","c4per","c3nfx","pastr","range")

# 1. read in the list of species we're including
setwd('C:/Users/cjcar/Dropbox/ViralChatter')
spname <- gsub('_',' ',read.csv('FinalSpeciesList.csv')[,2])

# 2. read in their habitat types document

iucndat <- read.csv('Land cover data/IucnHabitatData.csv')
LUHdat<-as.list(rep(NA,length(spname))) 
convcodes <- read.csv('Land cover data/IUCN_LUH_conversion_table.csv')

for (i in 1:length(LUHdat)){
  
  #get iucn habitat codes for the species  
  iucncodes<-iucndat[iucndat$name %in% spname[i],"code"]
  
  if(length(iucncodes)==0) {
    LUHdat[[i]]<-NULL
    names(LUHdat)[[i]]<-as.character(spname[i])
    next
  }
  #get matching LUH land use types  
  LUH<-unique(convcodes$LUH[match(as.numeric(as.character(iucncodes)),convcodes$IUCN_hab)])
  
  # when there is no match make LUH = NA
  # A single NA value means the species lives only in IUCN categories that have no 
  # mapping to LUH data (e.g., caves or ocean habitats)
  suppressWarnings(if (is.na(LUH) & length(LUH)==1) {
    LUHdat[[i]]<-NA
    names(LUHdat)[[i]]<-as.character(spname[i])
    next
  })
  # when there is a match    
  LUHvector<-unlist(strsplit(as.character(LUH),".", fixed=TRUE))
  
  LUHdat[[i]]<-LUHvector
  names(LUHdat)[[i]]<-as.character(spname[i])
}

########################################################

# 3. create a for loop

noL.auto <- c()
drop.marine <- c()

#length(spname)
for (i in 1:length(spname)) {
  spnamei <- spname[i]
  LUHi <- LUHdat[spnamei]
    if (all(is.na(LUHi))) { 
      noL.auto <- c(noL.auto, spnamei)
    } else { 
      LUHi[[1]] <- LUHi[[1]][!is.na(LUHi[[1]])]
      if (all(LUHi[[1]]=='MARINE')) {
            drop.marine <- c(drop.marine,spnamei)
      } 
    }
}

