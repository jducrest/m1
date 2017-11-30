#1.create a single queue
install.packages('queueing')
library(queueing)
q1 <- QueueingModel(NewInput.MM1(lambda=5,mu=7))
print(summary(q1))

#2.create a Jacson Network
n1 <- NewInput.MM1(lambda=8,mu=14,n=0)
n2 <- NewInput.MM1(lambda=0,mu=9,n=0)
n3 <- NewInput.MM1(lambda=6,mu=17,n=0)
n4 <- NewInput.MM1(lambda=0,mu=7,n=0)
m  <- c(0, 0.2, 0.56,0.24,1,0,0,0,0,0,0,0,0,0,0,0)
prob <- matrix(data=m, nrow=4, ncol=4,byrow=TRUE)
ojn1 <- NewInput.OJN(prob,n1,n2,n3,n4)

#3.compare MM1,MM2 and chelou network
n1 <- NewInput.MM1(lambda=8,mu=14)
n2 <- NewInput.MMC(lambda=8,mu=7,c=2)
print(summary(QueueingModel(n1)))
print(summary(QueueingModel(n2)))
n3 <- NewInput.MM1(lambda=8,mu=Inf,n=0)
n4 <- NewInput.MM1(lambda=0,mu=7,n=0)
n5 <- NewInput.MM1(lambda=0,mu=7,n=0)
m  <- c(0, 0.5, 0.5,0,0,0,0,0,0)
prob <- matrix(data=m, nrow=3, ncol=3,byrow=TRUE)
print(prob)
ojn1 <- NewInput.OJN(prob,n3,n4,n5)
print(summary(QueueingModel(ojn1)))
help(NewInput.MM1)
