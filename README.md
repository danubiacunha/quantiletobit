# quantiletobit
Fit for the log-symmetric quantile tobit regression model

The inputs of the function "tobitlogsymqreg.fit" are:

+ x: explanatory variables for Q
+ w: explanatory variables for phi
+ y: dependent variable
+ q: quantile of interest
+ xi:  extra parameter
+ status: vector of censored 0 or uncensored 1 values
+ cens: censoring point
+ link.Q: link for Q 
+ link.phi: link for phi
+ family: family ("Normal" for log-normal,"Student" for log-t,"Powerexp" for log-PE,"Sinh-normal" for EBS)
+ method: otimization method ("BFGS","Nelder-Mead","SANN","CG")
                          

