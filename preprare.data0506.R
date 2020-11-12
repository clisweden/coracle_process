rm(list = ls())
setwd("/Volumes/proj/coracle/process/")

require(hrbrthemes)

require(dplyr)
require(pbapply)
require(ggthemes)
require(igraph)
require(reshape2)

require(tidyverse)
require(feather)

## import feather-------
dpath="./data/covid19/"
node <- read_feather(paste0(dpath,"feather/nodes.feather"))
nodems <- unique(read_feather(paste0(dpath,"feather/nodems.feather")))
cnet <- read_feather(paste0(dpath,"feather/cnet.feather"))
fdate <- read_feather(paste0(dpath,"feather/fdate.feather"))
fdate = as.Date(as.POSIXct(fdate[[1]], 'GMT')) # update fdate to Date

# update country names
node$country[is.na(node$country)] = "Not available"
node$country[node$country=="China (Republic : 1949- )"] = "China"
node$country[node$country=="England" | node$country=="Scotland"] = "United Kingdom"
node$country[node$country=="Korea (South)" | node$country=="Korea"] = "Korea, Republic of"
node$country[node$country=="Iran"] = "Iran (Islamic Republic of)"

# update language
lan <- readLines("./predata/languages.tsv")
lan <- lan[2:length(lan)]
strsplit(lan,"\t") %>% as.data.frame %>% t -> lan
rownames(lan)=NULL
colnames(lan) = c("URL","code","label")
lan <- as.data.frame(lan)
node$language <- factor(lan$label[match(as.character(node$language),as.character(lan$code))])
rm(lan)
# format date
node$date <- as.Date(as.POSIXct(node$date, 'GMT'))

## output pmid.info.extented-------
pmid.info.extented <- node
write.csv2(pmid.info.extented,file="./coracle_data/pmid.info.extented.csv",quote=T,row.names = F)

pmid.info <- node[node$LitCovid==TRUE, ] %>% as.data.frame()

## output pmid.pub.type------
pmid.pub.type <- subset(nodems, pmid %in% pmid.info$pmid & type == "PublicationType")[,-2] %>% as.data.frame()
write.csv2(pmid.pub.type,file="./coracle_data/pmid.pub.type.csv",quote=T,row.names = F)

## output pmid.mesh------
pmid.mesh <- subset(nodems, pmid %in% pmid.info$pmid & type == "MESH")[,-2] %>% as.data.frame()

pmid.mesh <- rbind(pmid.mesh,cbind.data.frame(pmid=setdiff(pmid.info$pmid,pmid.mesh$pmid),value="not available"))
write.csv2(pmid.mesh,file="./coracle_data/pmid.mesh.csv",quote=T,row.names = F)

## haming distance
hamming2 <- function(X) {
  nn <- rownames(X)
  d <- matrix(0, nrow = nrow(X)*(nrow(X)-1)/2, ncol = 3)
  c = 0
  for ( i in 1:(nrow(X) - 1) ) {
    for ( j in (i + 1):nrow(X) ) {
      c = c+1
      if (sum((X[i,] + X[j,])==2)>0){
        d[c,] <- c(i,j,sum((X[i,] + X[j,])==2))
      }
    }
  }
  d <- d[d[,3]>0,]
  colnames(d) <- c("term1","term2","n")
  d <- as.data.frame(d)
  d[,1] <- rownames(X)[d[,1]]
  d[,2] <- rownames(X)[d[,2]]
  return(d)
}


## output cnet and relPMID
write.csv2(cnet,file="./coracle_data/cnet.csv",quote=T,row.names = F)
tmp0=cnet
tmp0$n=1
tmp <- spread(tmp0,target,n,fill=0) %>% t()
colnames(tmp) = tmp[1,]
tmp <- tmp[-1,]
class(tmp) = "numeric"
tmp2 <- hamming2(t(tmp))
relPMID <- tmp2
colnames(relPMID) <- c("PMID1","PMID2","N")  
write.csv2(relPMID,file="./coracle_data/relPMID.csv",quote=T,row.names = F)
rm(list=ls(pattern="tmp*"))

# generate relationships - rel
pmid.mesh$n=1
tmp <- spread(pmid.mesh,value,n,fill=0) %>% t()
colnames(tmp) = tmp[1,]
tmp <- tmp[-1,]
class(tmp) = "numeric"
tmp2 <- hamming2(tmp)

relMeSH = tmp2[order(tmp2$n,decreasing = T),]
write.csv2(relMeSH,file="./coracle_data/relMeSH.csv",quote=T,row.names = F)
rm(list=ls(pattern="tmp*"))



