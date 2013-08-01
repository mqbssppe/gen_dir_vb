#vb<-read.table("vb_result.m_alphas")
#vb<-vb[,2]
#n.pars<-8
#km<-kmeans(vb,centers = n.pars,iter.max = 1000)
#sizes<-numeric(n.pars)
#for (i in 1:n.pars)sizes[i]<-length(which(km$cluster==i))
#
#     write.table(km$cluster, file = "vb_clusters.txt", append = FALSE, quote = TRUE, sep = " ",
#                 eol = "\n", na = "NA", dec = ".", row.names = F,
#                 col.names = F, qmethod = c("escape", "double"),
#                 fileEncoding = "")



vb<-read.table("vb_result.m_alphas")
vb<-vb[,2]
K<-length(vb)
#x<-vb_plus/vb[1:(K-1)]

if (K<500){
clusters<-array(data = NA, dim = c(1,K-1))

for (i in 1:(K-1)){
clusters[1,i]<-i - 1
}

}else{

x<-vb[1:(K-1)]
index <- which(x==1)
if(length(index)>0){
x<-vb[1:(K-1)]/sum(vb)
y<-x[-index]
rest<-1:(K-1)
rest<-rest[-index]
vec <- c( 200 ) - 1
clusters<-array(data = NA, dim = c(length(vec),length(x)))
iter<-0
for (n.pars in vec ){
iter<-iter + 1
km<-kmeans(y,centers = n.pars,iter.max = 1000,nstart = 20)
clusters[iter,rest]<-km$cluster - 1
clusters[iter,index]<-rep(n.pars,length(index))
sizes<-numeric(n.pars)
for (i in 1:n.pars)sizes[i]<-length(which(km$cluster==i))
}
}else{
x<-vb[1:(K-1)]/sum(vb)
y<-x
rest<-1:(K-1)
vec <- c( 200 ) 
clusters<-array(data = NA, dim = c(length(vec),length(x)))
iter<-0
for (n.pars in vec ){
iter<-iter + 1
km<-kmeans(y,centers = n.pars,iter.max = 1000,nstart = 20)
clusters[iter,]<-km$cluster - 1
sizes<-numeric(n.pars)
for (i in 1:n.pars)sizes[i]<-length(which(km$cluster==i))
}
}



}



     write.table(clusters, file = "vb_clusters_new.txt", append = FALSE, quote = TRUE, sep = " ",
                 eol = "\n", na = "NA", dec = ".", row.names = F,
                 col.names = F, qmethod = c("escape", "double"),
                 fileEncoding = "")


