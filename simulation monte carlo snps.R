###############################################################
#x1....x1000 are snps in 50 regions(each region has 20 snps)
#There is a linkage equilibrium among the regions, and a LD within each region
# x970 and x990 are two snps that are correlated to the disease

###############################################################

#set the parameter

rm(list = ls())

n = 100      #subjects
m = 1000     #number of snps
alpha = 0.1
N = 10
k = 500

beg = 0       #LD = 0,0.05,0.1,0.15,0.2
endl = 0.2
byl = 0.2
num = (endl - beg)/byl+1

all_size_pro = rep(0,num)          #record the proposed/holm method's size/power
all_size_holm = rep(0,num)
all_power_pro = rep(0,num)
all_power_holm = rep(0,num)


############################################################################################
############################################################################################

for(LD_num in 1:num){
  
  LD = (LD_num-1)*byl

  sig_index_pro = matrix(0,nrow = N,ncol = m)     #record the significant   proposed method
  sig_index_holm = matrix(0,nrow = N,ncol = m)    #record the significant   holm
  
  
  f = function(x){      #"0" represents the major allele   "1" represents the minor allele
    if(x==0){y = rbinom(1,1,0.3-LD/0.7)}
    if(x==1){y = rbinom(1,1,0.3+LD/0.3)}
    return(y)
  }
  f = Vectorize(f)
  
  ########################################################################################
  ########################################################################################
  
  for(b in 1:N) {
    
    # generate 1000 biallelic snps
    
    Y = rep(0,n)          #document the subjects' disease condition
    
    case_snp_minor = rep(0,m)            #document the chi_square test statistics
    case_snp_major = rep(0,m)
    control_snp_minor = rep(0,m)
    control_snp_major = rep(0,m)
    
    for(sub in 1:n){
      A = matrix(0,50,20)
      A[,1] = rbinom(50,1,0.3)
      
      for (t in 2:20) {
        A[,t] = f(A[,t-1])
      }
      A = sapply(t(A),unlist)
      
      
      B = matrix(0,50,20)
      B[,1] = rbinom(50,1,0.3)
      
      for (t in 2:20) {
        B[,t] = f(B[,t-1])
      }
      B = sapply(t(B),unlist)
      
      X = (A+B)>0      #(A+B) is the number of the subject's minor allele
       
      Y[sub] = rbinom(1,1,1/(1+exp(3-0.8*X[970]-1*X[990])))     #document the subjects' disease condition
      
      
      ##################################################################################
      ##################################################################################
      
      #define the chi square test statistics
      
      if(Y[sub]==1){ case_snp_minor = case_snp_minor + (A+B)
                     case_snp_major = case_snp_major + rep(2,1000) - (A+B)                
      }
      if(Y[sub]==0){ control_snp_minor = control_snp_minor + (A+B)
                     control_snp_major = control_snp_major + rep(2,1000) - (A+B)
      }
      
    }
    
    
    case_snp_sum = case_snp_minor + case_snp_major                      
    control_snp_sum = control_snp_minor + control_snp_major
    minor_snp_sum = case_snp_minor + control_snp_minor
    major_snp_sum = case_snp_major + control_snp_major
    
    expected_case_minor = case_snp_sum*minor_snp_sum/(2*n)
    expected_case_major = case_snp_sum*major_snp_sum/(2*n)
    expected_control_minor = control_snp_sum*minor_snp_sum/(2*n)
    expected_control_major = control_snp_sum*major_snp_sum/(2*n)
    
    
    ###################################################################################
    ###################################################################################
    
    T = (case_snp_minor-expected_case_minor)^2/expected_case_minor + (case_snp_major - expected_case_major)^2/expected_case_major + (control_snp_minor - expected_control_minor)^2/expected_control_minor + (control_snp_major - expected_control_major)^2/expected_control_major      #obtain the test statistics
    
    P = 1-pchisq(T,df=1)    #obtain the p-value
    
    outc_sort = sort(P,index.return = TRUE)     #obtain the sort index
    
    p_sort = outc_sort$x      #sorted p-value
    
    T = T[outc_sort$ix]  #rearrange T ----> t(1) t(2).......
    
    
    
    
    ####################################################################################
    ####################################################################################
    
    #regenerate the data and the test statistics T_hat
    
    T_matrx = matrix(0,ncol = m , nrow = k)     #simulated test statistics matrix
    prob_H = rep(0,m)         
    
    for (i in 1:k) {     # 1 use the hypergeometric on case_minor_hat  2 use the binomial on case_minor_hat and control_minor_hat
      case_minor_hat = rbinom(length(case_snp_sum), case_snp_sum, 0.5)      #rhyper(m, case_snp_sum, control_snp_sum, minor_snp_sum)         rbinom(length(case_snp_sum), case_snp_sum, 0.5)
      case_major_hat = case_snp_sum - case_minor_hat
      control_minor_hat = rbinom(length(case_snp_sum), control_snp_sum, 0.5)    #minor_snp_sum - case_minor_hat       rbinom(length(case_snp_sum), control_snp_sum, 0.5)            #minor_snp_sum - case_minor_hat
      control_major_hat = control_snp_sum - control_minor_hat
      
      minor_sum_hat = case_minor_hat + control_minor_hat
      major_sum_hat = case_major_hat + control_major_hat
      
      expected_case_minor_hat = case_snp_sum*minor_sum_hat/(2*n)
      expected_case_major_hat = case_snp_sum*major_sum_hat/(2*n)
      expected_control_minor_hat = control_snp_sum*minor_sum_hat/(2*n)
      expected_control_major_hat = control_snp_sum*major_sum_hat/(2*n)
      
      
      temp = rep(0,m)
      
      T_hat = (case_minor_hat-expected_case_minor_hat)^2/expected_case_minor_hat + (case_major_hat - expected_case_major_hat)^2/expected_case_major_hat +(control_minor_hat - expected_control_minor_hat)^2/expected_control_minor_hat + (control_major_hat - expected_control_major_hat)^2/expected_control_major_hat  
      
     
       #calculate Pr(max(T_hat)>t(j))
      
      for (j in 1:m) {             
        temp[j] = as.numeric(max(T_hat[j:m])>T[j])
      }
      
      prob_H = prob_H + temp
      
      T_matrx[i,] = T_hat
      
    }
    
    prob_H = prob_H/k       
    
    
    ############################################################################
    ############################################################################
    
    i = 1
    
    while(prob_H[i]<=alpha&&i<=m){      # record the rejected snps proposed method
      
      sig_index_pro[b,i] = outc_sort$ix[i]  
      
      i = i+1
      
    }
    
    
    
    ############################################################################
    ############################################################################
    
    i = 1
    
    while(p_sort[i]<=alpha/(m-i+1)&&i<=m){        # record the rejected genes holm
      
      sig_index_holm[b,i] = outc_sort$ix[i]
      
      i = i+1
      
    }

   }
  
   FP_pro = apply(sig_index_pro, 1 , function(x){length(x[which(x>0&&x<=960)])>0})       #count numbers of false postive
   FP_pro = sum(FP_pro)

   TP_pro = apply(sig_index_pro, 1 , function(x){length(x[which(x>960)])>0})             #count numbers of true postive
   TP_pro = sum(TP_pro)

   size_pro = FP_pro/N
   power_pro = TP_pro/N



   FP_holm = apply(sig_index_holm, 1 , function(x){length(x[which(x>0&&x<=960)])>0})       #count numbers of false postive
   FP_holm = sum(FP_holm)

   TP_holm = apply(sig_index_holm, 1 , function(x){length(x[which(x>960)])>0})             #count numbers of true postive
   TP_holm = sum(TP_holm)

   size_holm = FP_holm/N
   power_holm = TP_holm/N


   all_size_pro[LD_num] = size_pro
   all_power_pro[LD_num] = power_pro
   all_size_holm[LD_num] = size_holm
   all_power_holm[LD_num] = power_holm
   
}         #the biggest loop over


 #############################################################################
 #############################################################################

 x = seq(from = beg , to = endl, by = byl)
 plot(x,all_power_pro,xlab = "Linkage Disequilibrium",ylab = "Size/Power",type = "l",ylim = c(0,1))
 lines(x,all_size_pro,col = 2)
 lines(x,all_power_holm,lty = 2)
 lines(x,all_size_holm,lty = 2,col = 2)

 legend("topleft",c("proposed","holm"),lty = c(1,2))







