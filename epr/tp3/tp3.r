M <- 14
rexp2 <- function(n, lambda)
{
  n <- n
  next_package <- 0
  cur_package_size <- 0
  buffer_size <- 1
  time <- (1:1000000)
  buffer_by_time <- (1:1000000)
  i <- 1
  while(n != 0 || buffer_size!=0 )
  {
    if(next_package == 0 && n != 0)
    {
        next_package <- ceiling(-log(runif(1,min=0,max=1))/lambda)
        buffer_size <- buffer_size + 1
        n <- n - 1
    }
    if(cur_package_size == 0 && buffer_size != 0)
    {
        cur_package_size <- ceiling(runif(1,min=0,max=M))
        buffer_size <- buffer_size - 1
    }
    cur_package_size <- max(cur_package_size - 1,0)
    next_package <- max(next_package - 1,0)
    
    buffer_by_time[i] <- buffer_size
    i <- i+1
    
  }
  plot(time[1:i-1],buffer_by_time[1:i-1])
  
}
rexp2(10000,0.13)

