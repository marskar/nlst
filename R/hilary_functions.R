

# Limit scientific notation
options(scipen=5)

# Function to create "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# Function to extract coefficient SEs
se.coef <- function(model)(sqrt(diag(vcov(model))))

# Function to nicely print ORs and CIs
print.ORs <- function(model)(round(cbind(OR=exp(coef(model)), lower=exp(coef(model)-1.96*se.coef(model)),
                                         upper=exp(coef(model)+1.96*se.coef(model))), 3))

# Function to print the range of a vector excluding outliers
range_without_outliers <- function(vec) {
  vec <- as.vector(vec)
  fences <- c(quantile(vec, 0.25) - 1.5*(quantile(vec, 0.75) - quantile(vec, 0.25)), quantile(vec, 0.75) + 1.5*(quantile(vec, 0.75) - quantile(vec, 0.25)))
  subset <- vec[vec>fences[1] & vec<fences[2]]
  c(min(subset), max(subset))
}


# The following are notes

# to view a subset of columns of a data frame [using dplyr] do: View(select(df, row1, row2))

# dput(): which allows you to dump an object in the form of R code.

# with() sets up an environment within which the R expression is evaluated. 
  # within() does the same thing but allows you to modify the data object used to create the environment.

# condition <- runif(1) > 0.5
  # if(condition) x <- 1 else x <- 2
  # OR you can do:
  # x <- if(condition) 1 else 2