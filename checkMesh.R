checkMesh <- function(raw) {
  #rm(list = ls())
  #setwd("/Volumes/proj/coracle/process/")
  ## import mesh.map data----
  if (!file.exists("./predata/mesh.map.rdata")) {
    library(feather)
    
    dpath = "./data/covid19/"
    mesh.map2020 <-
      read_feather(paste0(dpath, "feather/mesh_map_info/desc2020.feather"))
    temp <-
      read_feather(paste0(dpath, "feather/mesh_map_info/qual2020.feather"))
    temp2 <-
      read_feather(paste0(dpath, "feather/mesh_map_info/supp2020.feather"))
    mesh.map <- rbind(mesh.map2020, temp, temp2)
    mesh.map$termlow <- tolower(mesh.map$term)
    mesh.map$termlows <-
      gsub({
        "^ |,| |/|-|+|#|\""
      }, "", mesh.map$termlow)
    save(mesh.map, file = "./predata/mesh.map.rdata")
  }
  load("./predata/mesh.map.rdata")
  ## import checkedPMID.csv-----
  cPMID <- read.csv2("./coracle_mesh/checkedPMID.csv")
  # cPMID = matrix(c(1:3),nrow=3,ncol=1)
  # colnames(cPMID) <- "PMID"
  #write.csv2(cPMID,file="./latest_data/checkedPMID.csv",row.names = F)
  ## import mesh.annot.csv-----
  manu.mesh <-
    read.csv2("./coracle_mesh/manual.search.csv", stringsAsFactors = F)
  # mesh.annot <- read.csv2("./latest_data/mesh.annot.csv")
  load("./coracle_mesh/mesh.annot.rdata")
  #mesh.annot$CheckedBy <- "programe"
  manu.mesh.use <- manu.mesh[!is.na(manu.mesh$MeshID), ]
  mesh.annot <- rbind.data.frame(mesh.annot, manu.mesh.use)
  checkedMesh <-
    read.csv2("./coracle_mesh/checkedMesh.csv", stringsAsFactors = F)
  c1 <-
    as.character(mesh.annot$Term[mesh.annot$MeshMethod == "manual"])
  c2 <- union(manu.mesh$Term[nchar(manu.mesh$CheckedBy) > 0], c1)
  checkedMesh = union(c2, as.character(checkedMesh[, 1]))
  write.csv2(checkedMesh, file = "./coracle_mesh/checkedMesh.csv", row.names = F)
  rm(c1, c2, manu.mesh.use, manu.mesh)
  ## import pmid.mesh----
  # pmid.mesh <- read.csv2("https://raw.githubusercontent.com/clisweden/coracle_data/master/pmid.mesh.csv",stringsAsFactors = F)
  pmid.mesh = raw
  rm(raw)
  um.raw = as.data.frame(table(pmid.mesh$value), stringsAsFactors = F)
  colnames(um.raw)[1] = "Term"
  um.raw <- um.raw[um.raw$Term != "not available",]
  um.raw <- um.raw[order(um.raw$Freq, decreasing = T),]
  um <- um.raw$Term
  um <- setdiff(um, mesh.annot$Term)
  um.raw <- subset(um.raw, Term %in% um)
  # Map mesh from python to mesh.map------
  if (length(um) > 0) {
    ## clean up unique mesh term um
    tmp2 = which(grepl("<U+", um))
    tmp3 <- um[tmp2]
    tmp4 <- lapply(tmp3, function(a) {
      strsplit(a, "<|>")
    })
    tmp5 <- lapply(tmp4, function(a) {
      paste(a[[1]][which(sapply(a[[1]], function(b) {
        grepl("U+", b)
      }) == FALSE)], collapse = "")
    })
    um[tmp2] <- unlist(tmp5)
    rm(list = ls(pattern = "tmp*"))
    um2 <- gsub({
      "^ |,| |/|-|+|#|\""
    }, "", um)
    um.raw$MeshID <- NA
    um.raw$MeshMethod <- NA
    if (sum(is.element(um2, mesh.map$termlows)) > 0) {
      um.raw$MeshID <- mesh.map$id[match(um2, mesh.map$termlows)]
      um.raw$MeshMethod[!is.na(um.raw$MeshID)] = "pMap"
      um.raw$CheckedBy = "programe"
      mesh.annot <-
        rbind(mesh.annot, um.raw[!is.na(um.raw$MeshID),])
    }
    ## R search------
    
    msearch <- um.raw$Term[is.na(um.raw$MeshID)]
    if (length(msearch) > 0) {
      library(RISmed)
      pmid.s <-
        setdiff(subset(pmid.mesh, value %in% msearch)$pmid,
                unlist(cPMID))
      records <- EUtilsGet(pmid.s, db = "pubmed")
      cPMID <- union(unlist(cPMID), unlist(pmid.s))
      
      # save(records,file="rmap.rdata")
      # load("rmap.rdata")
      mesh.name <- Mesh(records)
      names(mesh.name) = pmid.s
      
      tmp0 = unlist(sapply(mesh.name, ncol))
      if (length(tmp0) > 0) {
        mesh.name2 <- mesh.name[names(tmp0)]
        mesh.name2 <- lapply(mesh.name2, unique)
        
        #if (length(mesh.name2) > 0) {
        r = mesh.name2
        temp4 = Reduce(rbind, Map(function(x, y)
          cbind(PMID = x, y), names(r), r))
        colnames(temp4)[1] = "PMID"
        colnames(temp4)[2] = "Term"
        rmap = temp4
        rmap$MeshID = mesh.map$id[match(tolower(rmap$Term), mesh.map$termlow)]
        rmap$MeshMethod[!is.na(rmap$MeshID)] = "rMap"
        
        if (length(!is.na(rmap$MeshID) > 0)) {
          load("./coracle_mesh/rmap.rdata")
          rmap$Freq = 0
          rmap$CheckedBy = "programe"
          rmap$PMID <- as.integer(as.character(rmap$PMID))
          rmap$Term <- as.character(rmap$Term)
          rmap$Type <- as.character(rmap$Type)
          tmp = rbind.data.frame(rmap.raw, rmap)
          rmap <- unique(tmp)
          rmap.raw <- rmap
          rm(tmp)
          save(rmap.raw, file = "./coracle_mesh/rmap.rdata")
          mesh.annot <-
            rbind(mesh.annot, unique(rmap[!is.na(rmap$MeshID), c("Term", "Freq", "MeshID", "MeshMethod", "CheckedBy")]))
          #um.raw <-
          #rbind(um.raw, unique(rmap[is.na(rmap$MeshID), c("Term", "Freq", "MeshID", "MeshMethod")]))
          pmid.mesh.r <-
            rmap[!is.na(rmap$MeshID), c("PMID", "Term")]
          colnames(pmid.mesh.r) = colnames(pmid.mesh)
          pmid.mesh <- rbind.data.frame(pmid.mesh, pmid.mesh.r)
        }
        
        #}
      }
      
    }
    mesh.annot <- unique(mesh.annot)
    save(mesh.annot, file = "./coracle_mesh/mesh.annot.rdata")
    pmid.mesh$MeshID <-
      mesh.annot$MeshID[match(pmid.mesh$value, mesh.annot$Term)]
    pmid.mesh$MeshForm <-
      mesh.map$term[match(pmid.mesh$MeshID, mesh.map$id)]
    
    cPMID <- union(unlist(cPMID), pmid.mesh$pmid)
    cPMID <- as.matrix(cPMID, nrow = length(cPMID), ncol = 1)
    #names(cPMID) = "PMID"
    write.csv2(cPMID, file = "./coracle_mesh/checkedPMID.csv", row.names = F)
    
    pmid.na <-
      pmid.mesh$pmid[is.na(pmid.mesh$MeshForm) &
                       pmid.mesh$value != "not available"]
    pmid.nna <- pmid.mesh$pmid[!is.na(pmid.mesh$MeshForm)]
    pmid.left <- setdiff(pmid.na, pmid.nna)
    
    tmp = subset(pmid.mesh, pmid %in% pmid.left)$value
    manu = as.data.frame(table(tmp))
    if (nrow(manu) > 0) {
      manu <- manu[order(manu$Freq, decreasing = T),]
      manu$MeshID = NA
      manu$MeshMethod = "manual"
      manu$CheckedBy = ""
      colnames(manu)[1] = "Term"
      manu <-
        manu[!is.element(as.character(manu$Term),
                         as.character(mesh.annot$Term)), ]
      manu <-
        manu[!is.element(as.character(manu$Term), checkedMesh), ]
      write.csv2(manu, file = "./coracle_mesh/manual.search.csv", row.names =F)
    }
    
    return(pmid.mesh)
  }
}
