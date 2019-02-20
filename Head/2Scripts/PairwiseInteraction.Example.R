rm(list=ls(all=TRUE))

VecA = c(1,1,0,0,0)
VecB = c(0,0,1,1,1)
POP = data.frame(VecA,VecB) # row names = individual numbers
# permutation is unlinking of the mutation from the background
ResVec = c()
for (i in 1:100000)
{
Per = POP; Per$VecA = sample(Per$VecA); Per$VecB = sample(Per$VecB);
N = nrow(Per[Per$VecA == 1 & Per$VecB == 1,])
ResVec = c(ResVec,N)
}
summary(ResVec)
mean(ResVec)/5 # 0.239 ~ 0.24


