rm(list = ls())
setwd("/Volumes/proj/coracle/process/")
t0 = Sys.time()
write.table(t0,file="./runStatus.txt")
require(tidyverse)
require(feather)
library(widyr)
library(dplyr)

## import the latest version pmid-------
rpath <- "./latest_data/"

raw.pmid.info.extended <- NULL
if (file.exists(paste0(rpath, "pmid.info.extended.csv"))) {
  raw.pmid.info.extended <-
    read.csv2(paste0(rpath, "pmid.info.extended.csv"), stringsAsFactors = F)
}

## import feather node-------
dpath = "./data/covid19/"
node <- read_feather(paste0(dpath, "feather/nodes.feather"))
#mesh.map <- read_feather(paste0(dpath, "feather/mesh_map_info.feather"))

#if (length(setdiff(node$pmid, raw.pmid.info.extended$pmid)) > 0) {
  raw.cnet = NULL
  raw.relPMID = NULL
  if (file.exists(paste0(rpath, "cnet.csv"))) {
    raw.cnet <- read.csv2(paste0(rpath, "cnet.csv"), stringsAsFactors = F)
    raw.relPMID <-
      read.csv2(paste0(rpath, "relPMIDfull.csv"), stringsAsFactors = F)
  }
  
  raw.pmid.mesh = NULL
  raw.relMeSH = NULL
  if (file.exists(paste0(rpath, "relMeSH.csv"))) {
    raw.relMeSH <-
      read.csv2(paste0(rpath, "relMeSHfull.csv"), stringsAsFactors = F)
    raw.pmid.mesh <-
      read.csv2(paste0(rpath, "pmid.mesh.csv"), stringsAsFactors = F)
  }
  
  
  ## import other feather files-------
  nodems <-
    unique(read_feather(paste0(dpath, "feather/nodems.feather")))
  cnet <- read_feather(paste0(dpath, "feather/cnet.feather"))
  write.csv2(cnet,
             file = "./coracle_data/cnet.csv",
             quote = T,
             row.names = F)
  fdate <- read_feather(paste0(dpath, "feather/fdate.feather"))
  fdate = as.Date(as.POSIXct(fdate[[1]], 'GMT')) # update fdate to Date
  write.csv2(fdate,
             file = "./coracle_data/fdate.csv",
             quote = T,
             row.names = F,
             col.names=F)
  
  # update country names
  node$country[is.na(node$country)] = "Not available"
  node$country[node$country == "China (Republic : 1949- )"] = "China"
  node$country[node$country == "England" |
                 node$country == "Scotland"] = "United Kingdom"
  node$country[node$country == "Korea (South)" |
                 node$country == "Korea"] = "Korea, Republic of"
  node$country[node$country == "Iran"] = "Iran (Islamic Republic of)"
  node$country[node$country == "Russia (Federation)"] = "Russia"
  
  # update language
  lan <- readLines("./predata/languages.tsv")
  
  lan <- lan[2:length(lan)]
  strsplit(lan, "\t") %>% as.data.frame %>% t -> lan
  rownames(lan) = NULL
  colnames(lan) = c("URL", "code", "label")
  lan <- as.data.frame(lan)
  node$language <-
    factor(lan$label[match(as.character(node$language), as.character(lan$code))])
  rm(lan)
  # format date
  node$date <- as.Date(as.POSIXct(node$date, 'GMT'))
  
  ## output pmid.info.extended-------
  pmid.info.extended <- node
  write.csv2(
    pmid.info.extended,
    file = "./coracle_data/pmid.info.extended.csv",
    quote = T,
    row.names = F
  )
  
  pmid.info <- node[node$LitCovid == TRUE,] %>% as.data.frame()
  
  ## output pmid.pub.type------
  pmid.pub.type <-
    subset(nodems, pmid %in% pmid.info$pmid &
             type == "PublicationType")[, -2] %>% as.data.frame()
  write.csv2(
    pmid.pub.type,
    file = "./coracle_data/pmid.pub.type.csv",
    quote = T,
    row.names = F
  )
  
  ## output pmid.mesh------
  source('./checkMesh.R')
  if (!is.null(raw.pmid.mesh)){
    raw.pmid.mesh <- raw.pmid.mesh[raw.pmid.mesh$value!="not available",]
  }
  pmid.mesh.raw <-
    subset(nodems, pmid %in% pmid.info$pmid &
             type == "MESH")[, -2] %>% as.data.frame()
  load("./coracle_mesh/rmap.rdata")
  rmap.use <- rmap.raw[,c("PMID","Term")]
  colnames(rmap.use) <- colnames(pmid.mesh.raw)
  
  pmid.mesh.raw <- rbind.data.frame(pmid.mesh.raw,rmap.use)
  pmid.mesh.raw <- unique(pmid.mesh.raw)
  rm(rmap.raw,rmap.use)
  
  pmid.mesh.raw = subset(pmid.mesh.raw, !(pmid %in% raw.pmid.mesh$pmid))
  pmid.mesh = NULL
   if (!is.null(pmid.mesh.raw)){
     pmid.mesh.full <- checkMesh(pmid.mesh.raw)
     
     pmid.mesh <- unique(pmid.mesh.full[,c("pmid","MeshForm")])
     colnames(pmid.mesh)[2]="value"
     pmid.mesh = pmid.mesh[!is.na(pmid.mesh[,2]),]
 }
pmid.mesh <- rbind.data.frame(raw.pmid.mesh,pmid.mesh)
pmid.mesh <- pmid.mesh[!is.na(pmid.mesh[,1]),]
  
  pmid.mesh<-
    rbind(pmid.mesh,
          cbind.data.frame(
            pmid = setdiff(pmid.info$pmid, pmid.mesh$pmid),
            value = "not available"
          ))
  
  write.csv2(pmid.mesh,
             file = "./coracle_data/pmid.mesh.csv",
             quote = T,
             row.names = F)
  
  #source('/Volumes/proj/coracle/process/ham_dist.R')
  # relPMID-------
  relPMID= pairwise_count(cnet, target, source)
  #relPMID <- ham_dist(raw.cnet, cnet, raw.relPMID, "integer")
  colnames(relPMID) <- c("PMID1", "PMID2", "N")
  write.csv2(relPMID,
             file = "./latest/relPMIDfull.csv",
             quote = T,
             row.names = F)
  relPMID <- relPMID[relPMID$N > 2, ]
  write.csv2(relPMID,
             file = "./coracle_data/relPMID.csv",
             quote = T,
             row.names = F)
  
  # relMeSH----
  pmid.mesh <- subset(pmid.mesh, !(value %in%"not available"))
  relMeSH = pairwise_count(pmid.mesh, value, pmid)
  colnames(relMeSH) <- c("term1", "term2", "n")
  # 
  # relMeSH <- ham_dist(raw.pmid.mesh, pmid.mesh, raw.relMeSH, "mesh")
  # colnames(relMeSH) <- c("term1", "term2", "n")
  
  write.csv2(relMeSH,
             file = "./latest/relMeSHfull.csv",
             quote = T,
             row.names = F)
  relMeSH <- relMeSH[relMeSH$n > 2, ]
  write.csv2(relMeSH,
             file = "./coracle_data/relMeSH.csv",
             quote = T,
             row.names = F)
  
  ## copy to ./latest_data/-----
  if (file.info("./coracle_data/relMeSH.csv")$ctime > file.info("./coracle_data/pmid.info.extended.csv")$ctime) {
    # find the files that you want
    setwd("./coracle_data/")
    list.of.files <- list.files("./")
    
    # copy the files to the new folder
    file.copy(list.of.files, "../latest_data/", overwrite = T)
  }
  
  # git update-----
  setwd("/Volumes/proj/coracle/process/coracle_data/")
  x2 <- "git add *.csv"
  x3 <- paste0('git commit -m "', fdate,'"')
  x4 <- "git push -u origin master"
  system(x2)
  system(x3)
  system(x4)

  t1 = Sys.time()
  print(t1 - t0)
  
#}
