# example 4.1
x<-matrix(c(0,0,0,0,0,0,1,2,0,2,
            1,0,0,1,1,1,2,1,0,1,
            1,0,0,1,0,0,1,1,1,1,
            1,1,1,1,0,1,2,1,1,1,
            0,1,1,1,1,1,2,1,0,1),nrow=5,byrow=T)
ones<-matrix(rep(1,nrow(x)),nrow=nrow(x))
sigma_g<-1
sigma_e<-1
miu<-1

###假设只有第一个标记有效应，服从N(0,1)分布，先抽样出标记效应
snp_eff<-rnorm(1,0,1) #0.3667677

###再抽样出随机残差的效应，其服从N(0,1)分布
residual<-rnorm(nrow(x),0,1)

###模拟出表型值
y<-ones*miu+x[,1]*snp_eff+residual

lhs<-rbind(cbind(5,t(ones)%*%x),cbind(t(x)%*%ones,crossprod(x)+diag(10)))
rhs<-rbind(t(ones)%*%y,t(x)%*%y)

estimate<-solve(lhs,rhs)

# example 4.2

n<-nrow(x)

p<-apply(x,2,function(x)sum(x)/(2*n))
sum_pq<-sum(2*p*(1-p))
p<-matrix(rep(2*p,each=5),nrow=5)
w<-x-p
gmat<-tcrossprod(w)/sum_pq
gmat2<-gmat+diag(n)*0.01
ginv<-solve(gmat2)
z<-diag(n)
lhs<-rbind(cbind(5,crossprod(ones,z)),cbind(crossprod(z,ones),crossprod(z)+ginv*1/(1*sum_pq)))
rhs<-rbind(crossprod(ones,y),crossprod(z,y))
estimate<-solve(lhs,rhs)

# ###利用约束极大似然(REML)估计参数
# loglike<-function(theta){
#   v<-x%*%t(x)*theta[1] + diag(nrow(x))*theta[2]
#   v.inv<-solve(v)
#   p<-v.inv-v.inv%*%ones%*%solve(t(ones)%*%v.inv%*%ones)%*%t(ones)%*%v.inv
#   p1<-t(y)%*%p%*%y
#   p2<-log(det(v))
#   p3<-log(det(t(ones)%*%v.inv%*%ones))
#   logvalue<--0.5*(p1+p2+p3)
#   return(-logvalue)
# }
# 
# theta<-c(1,1)
# par<-optim(theta,fn=loglike,method="L-BFGS-B",hessian=T,lower=1e-10,upper=1e+10)
# theta<-par$par
# (vara<-theta[1])
# (vare<-theta[2])
# (lambda<-vara/vare)


expect<-ones*miu
sigma<-tcrossprod(x[,1])*sigma_g+diag(nrow(x))*sigma_e

y<-rnorm(10,expect,sigma)
