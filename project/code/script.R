generate_hmm<-function(init,trans1,trans2,size)
{
  states=c(0,1)  
  res = matrix(0,nrow=2,ncol=size)
  res[1,1] = sample(states,1,prob=init)
  res[2,1] = sample(states,1,prob=trans2[res[1,1]+1,])
  for(i in 2:size)
  {
    res[1,i] = sample(states,1,prob=trans1[res[1,i-1]+1,])
    res[2,i] = sample(states,1,prob=trans2[res[1,i]+1,])
  }
  return(res)
}

trans1 = matrix(data=c(0.6,0.4,0.7,0.3),nrow=2,ncol=2,byrow=TRUE)
trans2 = matrix(data=c(0.1,0.9,0.8,0.2),nrow=2,ncol=2,byrow=TRUE)
init = c(0.5,0.5)

hmm = generate_hmm(init, trans1, trans2,10)

