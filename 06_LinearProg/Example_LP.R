# Use of Simplex Algorithm in R
# Should return the same results as runing same input data in Excel Solver (uses Simplex Algorithm)

# To clean up the memory of your current R session run the following line
rm(list=ls(all=TRUE))

# Solver
library(lpSolve)
# defining parameters
obj.fun <- c(7,884.09, 
             6,993.19,
             6,052.00,
             7,884.09,
             7,884.09,
             7,884.09,
             7,884.09,
             7,884.09,
             7,884.09,
             6,990.26,
             7,884.09)

constr <- matrix(c(6,993.19,
                   6,052.00,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   6,990.26,
                   7,884.09,
                   6,993.19,
                   6,052.00,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   7,884.09,
                   6,990.26,
                   7,884.09), ncol=1, byrow=TRUE)

constr.dir <- c("<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=",
                "<=")
                
rhs <- c(6,993.19,
         6,052.00,
         9,115.57,
         8,138.25,
         11,629.08,
         7,884.09,
         9,729.12,
         11,983.58,
         6,990.26,
         9,248.15,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09,
         7,884.09)

# solving model
prod.sol <- lp("max", obj.fun, constr, constr.dir, rhs, compute.sens = TRUE)

# accessing R output
prod.sol$objval # objective function value
prod.sol$solution # decision variables value



































































