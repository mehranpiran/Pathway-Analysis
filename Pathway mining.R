

library(igraph)
options(stringsAsFactors = F)
setwd("F:/SRC/KEGG Pathways")
load("KEGGedgelist.RData")


edgelist = edgelist[,c(2,4,3,1)]
edgelist[edgelist$`Interaction type` == "activation" , ][,2] = 1
edgelist[edgelist$`Interaction type` == "inhibition" , ][,2] = -1
edgelist$`Interaction type` = as.numeric(edgelist$`Interaction type`)
colnames(edgelist) = c("Source" , "Interaction" , "Target" , "Edge ID")
#save(edgelist , file = "edgelist.RData")


up = read.table("Up-regulated.txt" , sep = "\t")[,1]
down = read.table("Down-regulated.txt" , sep = "\t")[,1]
model.components = c(up,down)
g = graph_from_edgelist(as.matrix(edgelist[,c(1,3)]) , directed = T)
model.components = names(V(g))[which(names(V(g)) %in% model.components)]
model.components

duplicate <- function(x)
{
  x.vec <- apply(x, 1, paste, collapse = "")
  if(length(which(duplicated(x.vec))) > 0){
    return(x[-which(duplicated(x.vec)) , ])
  } else {return(x)}
}



all(model.components %in% unique(c(edgelist$Source , edgelist$Target)))


#edgelist[grep("^TJP" , edgelist[,1]) , ]
#edgelist[grep("^TJP" , edgelist[,3]) , ]



short_path = function(edgelist , model.components , input , N) {

g = graph_from_edgelist(as.matrix(edgelist[,c(1,3)]) , directed = T)
List.paths = list()
for(k in 1:length(model.components)){
  
  paths = all_shortest_paths(g, from = model.components[k], to = c(input) , mode = c("out"))
  index = which(lapply(lapply(paths$res , names) , length) <= N)
  
  if(length(index) > 0){
    
    
    paths1 = paths$res[index]
    p1 = lapply(paths1 , names)
    
    idx = c()
    for(i in 1:length(p1)){
      if(sum(p1[[1]] %in% model.components) > 1)
        idx = c(idx ,i)
    }
    
    if(length(idx)>0){
      p1 = p1[-idx]
    }
    
    if(length(p1) > 0){
      
      List1 = list()
      
      for(j in 1:length(p1)){
        p = p1[[j]]
        
        l = list()
        for(i in 1:(length(p)-1)){
          
          l[[i]] = (edgelist[edgelist[,1] == p[i] & edgelist[,3] == p[i+1],])
          
        }
        
        if(length(unique(unlist(lapply(lapply(l, dim) , function(x) x[1])))) != 1){
          a = unlist(lapply(lapply(l, dim) , function(x) x[1]))
          m = max(a)
          
          for(z in 1:length(l)){
            if((m - a[z])!= 0){
              mat = matrix(" " , m - a[z] , 4)
              colnames(mat)= colnames(l[[z]])
              l[[z]] = rbind(l[[z]],mat)
            }
          }
          List1[[j]] = do.call("cbind" , l)
          
        }else{
          List1[[j]] = do.call("cbind" , l)
        }
        
      }
      
      if(all(diff(do.call( "rbind" , lapply(List1 , dim))[,2])==0)){
        List.paths[[k]] = do.call("rbind" , List1)
      }else{
        List.paths[[k]] = List1  
      }
      
    }
    
  }
  
}

List.paths = List.paths[lapply(List.paths, length) > 0]

}


List.paths = short_path(edgelist ,  "SRC", model.components , 10)

List.paths = short_path(edgelist ,   model.components , "SRC", 10)

List.paths[[1]]

write.table(List.paths[[1]] , "DEGs_to_Src_pathways.txt" , sep = "\t"  , col.names = T , row.names = F , quote = F)

l.dataframe = List.paths[lapply(List.paths , class) == "data.frame"]
l.list =  List.paths[lapply(List.paths , class) == "list"]



if(length(l.dataframe) > 0){
  
  list2 = list()
  for(j in seq(4,10*4,4)){
    
    #assign(paste0("dframe" , j) , do.call("rbind" , l.dataframe[ which(unlist(lapply(lapply(l.dataframe , dim) , function(x) x[2])) == j) ]))
    
    list2[[j]] = do.call("rbind" , l.dataframe[ which(unlist(lapply(lapply(l.dataframe , dim) , function(x) x[2])) == j) ])
  }
  
  list2 = list2[lapply(list2, length) > 0]
  
}


length(list2)


d1 = duplicate(list2[[1]])
d2 = duplicate(list2[[2]])
d3 = duplicate(list2[[3]])
d4 = duplicate(list2[[4]])


write.table(d1 , "d1.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d2 , "d2.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d3 , "d3.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d4 , "d4.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)








if(length(l.list) > 0){
  
  
  for(j in 1:length(l.list)){
    l1 = list()
    idx1 = unlist(lapply(lapply(l.list[[j]], dim) , function(x) x[2]))
    idx2 = unique(unlist(lapply(lapply(l.list[[j]], dim) , function(x) x[2])))
    for(i in 1:length(idx2)){
      l1[[i]] = do.call("rbind" , l.list[[j]][idx1 == idx2[i]])
    }
    
    l.list[[j]] = l1
  }
  
  
  mx = max(unlist(lapply(l.list , length)))
  list1=list()
  for(i in 1:(mx+1)){
    list1[[i]] = list()
  }
  
  
  
  for(i in 1:mx){
    
    if(i < 3){
      
      for(j in c(4,8,12,16)){
        
        list1[[i]][[j]] = do.call("rbind" , lapply(l.list , function(x) x[[i]])[ unlist(lapply(lapply(lapply(l.list , function(x) x[[i]]) , dim) , function(x) x[2])) == j ])
        
      }
    }else if(i ==3){
      
      idx = which(lapply(l.list , length) == 3)
      for(j in c(4,8,12,16)){
        
        list1[[i]][[j]] = do.call("rbind" , lapply(l.list[idx] , function(x) x[[i]])[ unlist(lapply(lapply(lapply(l.list[idx] , function(x) x[[i]]) , dim) , function(x) x[2])) == j ])
        
      }
      
    }else if(i ==4){
      
      idx = which(lapply(l.list , length) == 3)
      for(j in c(4,8,12,16)){
        
        list1[[i]][[j]] = do.call("rbind" , lapply(l.list[idx] , function(x) x[[i]])[ unlist(lapply(lapply(lapply(l.list[idx] , function(x) x[[i]]) , dim) , function(x) x[2])) == j ])
        
      }
      
    }
    
  }
  
  list1 = list1[1:mx]
  
  for(i in 1:mx){
    list1[[i]] = list1[[i]][lapply(list1[[i]], length) > 0]
  }
  
  
  
  list1.1 = list()
  for(j in c(4,8,12,16)){
    list1.2 = list()
    for(i in 1:mx){
      list1.2[[i]] = as.data.frame(list1[[i]][unlist(lapply(lapply(list1[[i]],dim) , function(x) x[2]==j))])
    }
    list1.1[[j]] = do.call("rbind" , list1.2)
  }
  
  list1 = list1.1[lapply(list1.1 , length) > 0]
  
}

l.list[[2]]
l.list = l.list[[1]]


d1 = duplicate(l.list[[1]])
d2 = duplicate(l.list[[2]])
d3 = duplicate(l.list[[3]])
d4 = duplicate(l.list[[4]])
d5 = duplicate(l.list[[5]])
d6 = duplicate(l.list[[6]])

write.table(d1 , "d1.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d2 , "d2.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d3 , "d3.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d4 , "d4.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d5 , "d5.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)
write.table(d6 , "d6.txt" , sep = "\t" , col.names = T , quote = F , row.names = F)







