BiocManager::install("KEGGgraph")
library(KEGGgraph)

#Define path
kgfiles<-list.files("F:/SRC/Results/KEGG/Files",full.names=T)

#Read kgml files
KGMLGraphs=lapply(kgfiles,parseKGML2Graph,genesOnly=F)
length(KGMLGraphs)
head(KGMLGraphs)


l = lapply(KGMLGraphs,getKEGGedgeData)
length(l)

a = NULL
b = list()

listedgelists = list()

for(i in 1:length(l)){
  if(is.list(l[[i]]) & length(l[[i]]) == 0){
    listedgelists[[i]] = c(0,0,0)
  } else 
  {
    d = matrix(0,length(l[[i]]),3)
    d = as.data.frame(d)
    for(j in 1:length(l[[i]])){
      d[j,1] = l[[i]][[j]]@entry1ID
      d[j,2] = l[[i]][[j]]@entry2ID
      c1 = l[[i]][[j]]@subtype$subtype
      
      if(class(c1) == class(a)){
        d[j,3] =  "none"
      } else
      {
        d[j,3] = l[[i]][[j]]@subtype$subtype@name
      }
      
    }
    
    listedgelists[[i]] = d
  }
}

head(listedgelists[[1]])

table = do.call("rbind", listedgelists)
dim(table)







table1 = table[!duplicated(table),]

idx1 = which(table1[,3] == "activation")

idx2 = which(table1[,3] == "inhibition")

idx = c(idx1,idx2)

table2 = table1[idx,]

dim(table2)





idx3 = grep("^hsa" , table2[,1])

idx4 = grep("^hsa" , table2[,2])

idx5 = intersect(idx3,idx4)

table3 = table2[idx5,]

dim(table3)

sum(!grepl("^hsa" , table3[,2]))
sum(!grepl("^hsa" , table3[,1]))






kgid1 = as.character(table3[,1]) 

kgid2 = as.character(table3[,2])

geneid1 <- translateKEGGID2GeneID(kgid1)

geneid2 <- translateKEGGID2GeneID(kgid2)

any(is.na(geneid2))





require(org.Hs.eg.db)

genesymbol1 <- sapply(mget(geneid1, org.Hs.egSYMBOL, 
                           ifnotfound=NA), "[[",1)
length(genesymbol1)

any(is.na(genesymbol1))

genesymbol2 <- sapply(mget(geneid2, org.Hs.egSYMBOL, 
                           ifnotfound=NA), "[[",1)
length(genesymbol2)

any(is.na(genesymbol2))






edgelist = table3 

edgelist[,1] = genesymbol1

edgelist[,2] = genesymbol2

dim(edgelist)










id1 = paste0("E0000" , 1:9 )
id2 = paste0("E000" , 10:99)
id3 = paste0("E00" , 100:999)
id4 = paste0("E0" , 1000:length(edgelist[,1]))
ID = c(id1,id2,id3,id4)
edgelist = cbind(ID,edgelist)
colnames(edgelist) = c("ID" , "Gene1" , "Gene2" , "Interaction type")
edgelist = apply(edgelist , 2 , as.character)
edgelist = as.data.frame(edgelist , stringsAsFactors = F)
class(edgelist)


edgelist = edgelist[,c(2,4,3,1)]
edgelist[edgelist$`Interaction type` == "activation" , ][,2] = 1
edgelist[edgelist$`Interaction type` == "inhibition" , ][,2] = -1
edgelist$`Interaction type` = as.numeric(edgelist$`Interaction type`)
colnames(edgelist) = c("Source" , "Interaction" , "Target" , "Edge ID")
save(edgelist , file = "edgelist.RData")