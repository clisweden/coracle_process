ham_dist <- function(rawpair, newpair, rawdist, tag) {
  # newpair=pmid.mesh
  # rawpair=NULL
  # rawdist=NULL
  colnames(newpair) <- c("source", "target")
  if (!is.null(rawpair)) {
    colnames(rawpair) <- c("source", "target")
    
    colnames(rawdist) <- c("x", "y", "n")
    npair <- setdiff(newpair,rawpair)
    newdist <- rawdist
  } else{
    npair <- newpair
   # newdist <- rawdist
  }
  
  
  hamming_binary <- function(X, Y = NULL) {
    if (is.null(Y)) {
      D <- t(X) %*% X
    } else {
      t(X) %*% Y
    }
  }

  if (nrow(npair)>0){
    tid <- union(newpair$source, newpair$target)
    newdist = NULL
    
    if (length(tid) <= 20000 | is.null(rawpair)) {
      tmp0 = newpair
      tmp0$n = 1
      tmp <- spread(tmp0, target, n, fill = 0) %>% as.matrix()
      rownames(tmp) = tmp[, 1]
      tmp <- tmp[,-1]
      tmp2 <- hamming_binary(tmp)
    }
    else{
      unet = setdiff(newpair, rawpair)
      newid <- unique(unet$target)
      newdist <-
        rawdist[!is.element(rawdist$x, newid) &
                  !is.element(rawdist$y, newid),]
      
      
      tmp0 = subset(rawdist, x %in% newid | y %in% newid)
      tmpp <- union(tmp0$x, tmp0$y) %>% union(., newid)
      
      tmp0 = subset(newpair, target %in% tmpp)
      tmp0$n = 1
      tmp <- spread(tmp0, target, n, fill = 0) %>% as.matrix()
      rownames(tmp) = tmp[, 1]
      tmp <- tmp[,-1]
      tmp2 <-
        hamming_binary(tmp[, which(is.element(newid, colnames(tmp)))], tmp)
    }
    tmpx=which(tmp2>0,arr.ind = T)
    tmpx2=cbind.data.frame(x=rownames(tmp2)[tmpx[,1]],y=colnames(tmp2)[tmpx[,2]],n=tmp2[tmpx],stringsAsFactors =F)
    #tmpx2[,1] <- as.integer(tmpx2[,1])
    #tmpx2[,2] <- as.integer(tmpx2[,2])
    tmpx2 = tmpx2[(tmpx2[,1]!=tmpx2[,2]),]
   # tmp2 <- cbind.data.frame(x = rownames(tmp2), tmp2)
    #tmp3 <- gather(tmp2, colnames(tmp2)[-1], key = "y", value = "n")
    #tmp4 <- tmp3[tmp3$n > 0,]
    tmp4=tmpx2
    if (tag == "integer") {
      tmp4$x <- as.integer(as.character(tmp4$x))
      tmp4$y <- as.integer(tmp4$y)
    } else
      {tmp4$x <- as.character(tmp4$x)
    tmp5 = which(tmp4[, 1] > tmp4[, 2])
    tmp4[tmp5, 2:1] = tmp4[tmp5, 1:2]
    tmp5 = which(tmp4[, 1] == tmp4[, 2])
    tmp4 <- tmp4[-tmp5, ]
    tmp4 <- unique(tmp4)}
    
    newdist <- rbind.data.frame(tmp4, newdist)
    newdist <- unique(newdist)
    newdist <- newdist[newdist[,1]!=newdist[,2],]
    rownames(newdist) = NULL
  }
  
  return(newdist)
}
