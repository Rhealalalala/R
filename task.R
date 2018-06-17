tabuDiffStart=function(A)  ##A矩阵为样本空间
 {
   A <- as.matrix(read.table("/Users/xuehanqing/Desktop/Breast-FPKM-Sort-Mixture-Unique.txt"))
   n_gene = dim(A)[1]  ##number of genes
   n_sample = dim(A)[2]  ##number of samples
   p_value = rep(0,n_gene)
   tumor_ = rep(0,51)
   normal_ = rep(0,51)
   p_threshold = 0.05
   count = 0
##log2(x+1)
   A = logb( A+1, 2)

   for (i in 1:n_gene)
   {
      tumor_ = A[i,1:51]
      normal_ = A[i,52:102]
      p_value[i] = ks.test(tumor_,normal_)[2]
      if( p_value[i] < p_threshold )
        count = count+1
   }

##提取有明显差异的
   A_ = matrix(0,count,n_sample)
   index = which(p_value<p_threshold)
   A_ = A[index,]
   
##   return (list(p = p_value,a = A_))
## }

##z-score
   for (i in 1:count)
      scale(A_[i,], center=TRUE, scale=TRUE) 


   S_tumor = matrix(0,count,count)  ##tumor基因的相似度矩阵
   S_normal = matrix(0,count,count)  ##normal下基因的相似度矩阵
   S_mix = matrix(0,count,count)

   for( i in 1:count)
   {
      for( j in i:count)
      {
          S_tumor[i][j] = cor(A_[i,1:51],A_[j,1:51])
          S_normal[i][j] = cor(A_[i,52:102],A_[j,52:102])
          S_mix[i][j] = cor(A_[i,1:51],A_[j,52:102])
      }
   }

   for( i in 1:count)
   {
      for( j in 1:i)
      {
         S_tumor[i][j] = S_tumor[j][i]
         S_normal[i][j] = S_normal[j][i]
         S_mix[i][j] = S_mix[j][i]
      }
   }

   ##求样本均值
   Ave_tumor = rowMeans(A_[,1:51])
   Ave_normal = rowMeans(A_[,52:102])
   tumor <- matrix(0, ncol = 51, nrow = nrow(A_))
   normal <- matrix(0, ncol = 51, nrow = nrow(A_))

   ##线性回归
   for( i in 1:51)
   {
      model <- lm(A_[,51+i] ~ Ave_normal)
      normal[,i] <- model$residuals
      model <- lm(A_[,i] ~ Ave_tumor)
      tumor[,i] <- model$residuals
   }

   ##最大似然估计
   for( i in 1:51)
   {
      p>qnorm(0.95,mean(normal[i,]),var(normal[i,]))

      sigma_normal[i] = var(normal[i,])
      miu_tumor[i] = mean(tumor[i,])
      sigma_normal[i] = var(tumor[i,])
   }



   




   
   return(list(st = S_tumor, sn = S_normal, sm = S_mix))
 } 