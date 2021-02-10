


suppressMessages(library(CVXR)) #convex optimization


#given a sequence of infected, compute matrix A
#I, sequence of infected
#T, generative period, 0, 1, 2, T
fillMatrix <- function (I, T) {
  
  
  #think on D-indices from T-1 to |I|-2
  
  #init matrix A
  A = matrix(0, nrow = 1, ncol = T + 1)
  #first part
  A[1,1:T] = I[seq(from=T, to=1 , by= -1)]
  #diagonal
  A[1,T + 1] = - I[T + 1]
  
  
  return (A)
}


#compute the generative period values and Rt
generative <- function (A, T){   
  #place holder for the generative function 
  w <- Variable(T)
  R <- Variable(1)
  #creat unknow variable
  beta <- vstack(w,R)
  
  #define constraints
  constr <- list (sum(w) == 1, w >= 0, R >= 0)
  
  #define objective function
  obj <- sum_squares(A %*% beta)
  #state minimization problem 
  prob <- Problem(Minimize(obj), constr)
  #compute the minimum
  result <- solve(prob)
  #print(result$status)
  #print(result$value)
  
  x = result$getValue(beta)
  error = (A %*% x)^2
  
  res = list(x = x, error = error)
  #print(result$getValue(beta))
  return (res)
}