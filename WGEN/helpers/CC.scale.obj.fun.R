
CC.scale.obj.fun <- function(param,q.old,mu.old,param.old,qq,perc.q,perc.mu,emission.type) {
  
  if (emission.type=="gamma") {
    
    shape.old <- param.old[1]
    rate.old  <- param.old[2]
    
    #perturb the gamma parameters in a way that gaurentees the perc.mu change
    shape.new <- shape.old*perc.mu*param
    rate.new <- rate.old*param
    q.new <- qgamma(qq,shape.new,rate.new)
    #calcualte error between percent change and target percent change
    e.q <- (q.new/q.old - perc.q)  
    return(e.q^2)
  }
  
  if (emission.type=="mixed.exp") {

    l1.old <- param.old[1]
    l2.old  <- param.old[2]
    m.old  <- param.old[3]
    
    #perturb the 2 exponential parameters
    l1.new <- l1.old + param[1]  #l1.old*param[1]
    l2.new <- l2.old + param[2]  #l2.old*param[2]
    m.new <- m.old + param[3]

    mu.new <- m.new/l1.new+(1-m.new)/l2.new
    q.new <- qmixedexp(qq,l1.new,l2.new,m.new)
    #calcualte error between percent change and target percent change
    e.mu <- (mu.new/mu.old - perc.mu)  
    e.q <- (q.new/q.old - perc.q)  
    return(e.mu^2 + e.q^2)
    
  }  
  
}