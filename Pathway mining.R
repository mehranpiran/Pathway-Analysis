

library(igraph)
options(stringsAsFactors = F)

# import KEGG edgelist
setwd("")
load("KEGG.edgelist.RData")



# model.components object contains all DEGs affected by Src overactivation
# select components that are present in KEGG edgelist
up = read.table("Up-regulated.txt" , sep = "\t")[,1]
down = read.table("Down-regulated.txt" , sep = "\t")[,1]
model.components = c(up,down)
g = graph_from_edgelist(as.matrix(edgelist[,c(1,3)]) , directed = T)
model.components = names(V(g))[which(names(V(g)) %in% model.components)]
model.components




all(model.components %in% unique(c(edgelist$Source , edgelist$Target)))
# True



# Writing the main function

# Function accepts four arguments:
# edgelist: edgelist must contain 4 columns including source genes, interaction type, target gene and any kind of annotation such as IDs (edgeID, Pubmed ID and so on).
# inputs: inputs are the genes you want to see their relashion with targets
# model.components: they are differentially expressed genes influenced by target genes
# N: this argument defines the maximum lenght of pathways
# there must be no intersect between inputs and model.components

short_path = function(edgelist , inputs , model.components , N) {

g = graph_from_edgelist(as.matrix(edgelist[,c(1,3)]) , directed = T)
List.paths = list()
for(k in 1:length(model.components)){
  
  paths = all_shortest_paths(g, from = c(inputs), to = model.components[k] , mode = c("out"))
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

# List.paths contains different pathways starts with SRC and ends with one of the SRC targets based on the edge information in KEGG edgelist
                            
# The rest of the codes are for ordering pathways from shortest to longest in tabular forms
                            
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

                                                                    

duplicate <- function(x)
{
  x.vec <- apply(x, 1, paste, collapse = "")
  if(length(which(duplicated(x.vec))) > 0){
    return(x[-which(duplicated(x.vec)) , ])
  } else {return(x)}
}

length(list2)
# 4
                                                                    
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

length(list1)
# 6


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







