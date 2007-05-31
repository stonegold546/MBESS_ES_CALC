`cor2cov` <-
function(cor.mat, sd)
{
if(dim(cor.mat)[1]!=dim(cor.mat)[2]) stop("'cor.mat' should be a square matrix")

n<-sqrt(length(cor.mat))
if(n!=length(sd)) stop("The length of 'sd' should be the same as the number of rows of 'cor.mat'")

if(length(sd[sd>0])!= n) stop("The elements in 'sd' shuold all be positive")

for(j in 1:n)
    {for(i in 1:n)
        {if(i==j)
            {if(cor.mat[i,j]!=1) stop ("The elements on the main diagonal of 'cor.mat' should all be 1")
            }
    
        if(cor.mat[i,j]!=cor.mat[j,i])
            {if (cor.mat[i,j]!=0 & cor.mat[j,i]!=0) stop("'cor.mat' should be either symmetric or triangular")
            if (cor.mat[i,j]!=0 & cor.mat[j,i]==0) cor.mat[j,i]<-cor.mat[i,j]
            if (cor.mat[j,i]!=0 & cor.mat[i,j]==0) cor.mat[i,j]<-cor.mat[j,i] 
            }  
        }  
    }

diag(sd)  %*% cor.mat  %*% diag(sd)
}

