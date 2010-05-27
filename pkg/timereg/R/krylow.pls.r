krylow.pls<-function(D,d,dim)
{
R=d;  Sxxsxy=R;
if (dim>=2)
for (i in 2:dim)
{
Sxxsxy=D %*% Sxxsxy  ;
R=cbind(R,Sxxsxy);
}
beta= R %*% solve(t(R) %*% D %*% R) %*% t(R) %*% d
beta=  beta;
return(list(beta=beta))
}
