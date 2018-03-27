###############################################################
# x1....x2/n = 1 in treatment group
# x2/n+1.....xn = 0 in control group
###############################################################
rm(list = ls())

n=50      #Samples
m=50     #Genes

alpha=0.1

N = 1000       #experiment numbers

beg = 0.1       #parameter for sigma_ksai_sd
endl = 0.9
byl = 0.1
num = (endl - beg)/byl+1

all_size_pro = rep(0,1+((endl-beg)/byl))          #record the proposed/holm method's size power
all_size_holm = rep(0,1+((endl-beg)/byl))
all_power_pro = rep(0,1+((endl-beg)/byl))
all_power_holm = rep(0,1+((endl-beg)/byl))



for(t in 1:num){
  
  sigma_ksai_sd = sqrt(beg + (t-1)*byl)
  sigma_epis_sd = sqrt(1-sigma_ksai_sd^2)
  
  
  sig_index_pro = matrix(0,nrow = N,ncol = m)     #record the significant   proposed method
  sig_index_holm = matrix(0,nrow = N,ncol = m)    #record the significant   holm
  
  for(b in 1:N) {
    
    obs_outc = matrix(0,ncol = n,nrow = m)         # observed matrix
    
    A = matrix(rnorm(m*n,0,sigma_epis_sd),nrow = m) 
    
    B = matrix(rnorm(n,0,sigma_ksai_sd),ncol = n,nrow = m,byrow = TRUE)
    
    f = function(x) 1.1*(x-45)/5
    beta = c(rep(0,45),f(seq(46,50,by=1)))    #vector beta
    
    X = c(rep(1,n/2),rep(0,n/2))
    
    C = beta%*%t(X)
    
    obs_outc = A+B+C     #obtain the observed matrix
    
    
    ###################################################################################
    ###################################################################################
    
    
    U = rowSums(obs_outc)    #
    
    V = apply(obs_outc,1, function(x){sum(x*x)})    #obtain the varianc e
    
    T = U*(1/V)*U      #obtain the test statistics
    
    P = 1-pchisq(T,df=1)    #obtain the p-value
    
    outc_sort = sort(P,index.return = TRUE)     #obtain the sort index
    
    p_sort = outc_sort$x      #sorted p-value
    
    T = T[outc_sort$ix]  #rearrange T ----> t(1) t(2).......
    
    
    
    
    ####################################################################################
    ####################################################################################
    
    
    k = 1000    #experiment numbers
    T_matrx = matrix(0,ncol = m , nrow = k)     #simulated test statistics matrix
    prob_H = rep(0,m)         
    
    for (i in 1:k) {
      G = rnorm(n)       # standard normal vector G
      
      temp = rep(0,m)
      
      U_hat = apply(obs_outc, 1 , function(x) sum(x*G))
      
      T_hat = U_hat*(1/V)*U_hat     
      
      for (j in 1:m) {
        temp[j] = as.numeric(max(T_hat[j:m])>T[j])
      }
      
      prob_H = prob_H + temp
      
      T_matrx[i,] = T_hat
      
    }
    
    prob_H = prob_H/k       #calculate Pr(max(T_hat)>t(j))
    
    
    ############################################################################
    ############################################################################
    
    i = 1
    
    while(prob_H[i]<=alpha&&i<=50){      # record the rejected genes  proposed method
      
      sig_index_pro[b,i] = outc_sort$ix[i]  
      
      i = i+1
      
    }
    
    
    
    ############################################################################
    ############################################################################
    
    i = 1
    
    while(p_sort[i]<=alpha/(m-i+1)&&i<=50){        # record the rejected genes holm
      
      sig_index_holm[b,i] = outc_sort$ix[i]
      
      i = i+1
      
    }
    
  }  
  FP_pro = apply(sig_index_pro, 1 , function(x){length(x[which(x>0&&x<=45)])>0})       #count numbers of false postive
  FP_pro = sum(FP_pro)
  
  TP_pro = apply(sig_index_pro, 1 , function(x){length(x[which(x>45)])>0})             #count numbers of true postive
  TP_pro = sum(TP_pro)
  
  size_pro = FP_pro/N
  power_pro = TP_pro/N
  
  
  
  FP_holm = apply(sig_index_holm, 1 , function(x){length(x[which(x>0&&x<=45)])>0})       #count numbers of false postive
  FP_holm = sum(FP_holm)
  
  TP_holm = apply(sig_index_holm, 1 , function(x){length(x[which(x>45)])>0})             #count numbers of true postive
  TP_holm = sum(TP_holm)
  
  size_holm = FP_holm/N
  power_holm = TP_holm/N
  
  
  all_size_pro[t] = size_pro
  all_power_pro[t] = power_pro
  all_size_holm[t] = size_holm
  all_power_holm[t] = power_holm
  
}         #the biggest loop over 


#############################################################################
#############################################################################

x = seq(from = beg , to = endl, by = byl)
plot(x,all_power_pro,xlab = "Intra-Class Correlation",ylab = "Size/Power",type = "l",ylim = c(0,1))
lines(x,all_size_pro)
lines(x,all_power_holm,lty = 2)
lines(x,all_size_holm,lty = 2)

legend("topleft",c("proposed","holm"),lty = c(1,2))







