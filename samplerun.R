# sample run 
require('quantreg')
set.seed(123)
n=300
x=rbinom(n, 2, 0.2)
y=rnorm(n, mean=0, sd=1)
z=cbind(rbinom(n, 1, 0.3), rnorm(n, mean=2, sd=2))
tau=c( 0.25, 0.5, 0.75)
dat <- data.frame(x=x,z=z,y=y)

rst <- QRank(y~x+z,dat,tau)
print(rst) # same with the old print out

rst <- QRank(y~factor(x)*z.1+z.2,dat,tau) # with factor and interaction
print(rst)

rst <- QRank(y~.-1,dat,tau) # all variables, no intercept
print(rst)

# wrong input
rst <- QRank(y~x+z,tau)
rst <- QRank(y~x+z,dat,tau = c(0.5,0.7,1.2))
print(rst)


# plot 
plot_Qrank(y~.,'z.2',dat) # z.2 is continuous, plot at 0.25,0.5,0.75 quantiles
plot_Qrank(y~.,'x',dat) # x is treated as continuous, however x only has three values, so plot at 0.25,0.5,0.75 quantiles of unique values 
plot_Qrank(y~factor(x)+z.1+z.2,'x',dat) # x is treated as categorical, plot at all levels.
plot_Qrank(y~factor(x)+z.1+z.2,'x',dat,tau = seq(0.02, 0.98, by=0.05), ylab='hahah',main='title!') # change tau

plot_Qrank(y~.,'z.2',dat,linequantile = c(0.01,0.1,0.5,0.9,0.99)) # user defined quantile for z.2. I always keep the mid line black, and other lines grey.

plot_Qrank(y~.,'z.2',dat,linequantile = c(0.1,0.5,0.9), subjectID = c(1,15,19)) # user defined subjects, all lines are black.

plot_Qrank(y~x+z.2,'z.2',dat,linequantile = c(0.1,0.5,0.9), newdata = data.frame(x=c(0,1,2),z.2=c(-1,2,5))) # user defined data, all lines are black.

#wrong input
ndat <- dat
ndat$x <- factor(dat$x)
plot_Qrank(y~x+z.1,'x',ndat)
plot_Qrank(y~factor(x)+z.1,'x',ndat) # that is right! 

ndat <- dat
ndat[2:3,2] <- NA # subject 2,3 have Na
plot_Qrank(y~x+z.1,'x',ndat,subjectID = c(2,3,8,10))
