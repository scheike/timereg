context("Reshaping data")

test_that("fast reshape I", {
  m <- lvm(y~x)
  set.seed(1)
  d <- sim(m,100)
  ## Just a stratified analysis
  e <- estimate(list("Group 1"=m,"Group 2"=m),list(d,d))
  expect_equivalent(coef(e)[[1]][1:2,1],coef(lm(y~x,d)))
  expect_equivalent(coef(e)[[2]][1:2,1],coef(lm(y~x,d)))
})


fast.reshape(fast.reshape(d,var=c("y","z","w")),id="id",var=c("y","z","w"))


library(mets)
x <- matrix(1:10,5,2)
x[3,2] <- 8
x[3,2] <- NA
cluster <- c(1,1,2,2,3)
x <- cbind(x,cluster)
x
ud <- fast.reshape(data.frame(x),"cluster")
ud
###
out=cluster.index(cluster)
out
###
ud <- faster.reshape(x,cluster)
ud
ud <- faster.reshape(data.frame(x),cluster)
ud
###
colnames(x) <- c("y1","y2","cluster")
x
ud <- fast.reshape(data.frame(x),"cluster")
ud
###
num <- c(2,1,1,2,2)
x
out <- faster.reshape(x,cluster)
out
out <- faster.reshape(x,cluster,num=num)
out

### problem med 

