library( mirt)

#eliminar todas las variables
rm(list = ls(all = TRUE))

data(LSAT7)
conv = function(LSAT7){
LsatMio=matrix(rep(0,1000),ncol=5,nrow=1000)
LSAT72 = as.matrix(LSAT7)
m = 1
for(i in 1:32){
  inds = LSAT7[i,6]
  for(j in 1:inds){
    LsatMio[m,] = LSAT72[i,-6]
    m = m + 1
  }
}
LsatMio 
}

LsatMio = conv(LSAT7)

write.table(LsatMio,file = "C:\\Users\\ALEJANDRO\\Documents\\Noahrm4\\lsatR.txt",sep= " ",
            row.names= FALSE, col.names= FALSE)

