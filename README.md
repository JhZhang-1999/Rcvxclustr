## Rcvxclustr

**Robust Convex Clustering Algorithm Implemented in R**

***

### Description

This package implements the robust convex clustering problem proposed by [Liu, Sun and Tan (2019)](https://arxiv.org/abs/1906.09581v2). Classical approaches towards convex clustering solves a convex optimization problem with the cost function being a squared loss plus a fusion penalty that encourages the estimated centroids for observations in the sames cluster to be identical. Those approaches are not robust to arbitrary outliers, and when data are contaminated, they fail to identify the correct cluster relationships. This proposed robust convex clustering algorithm, applying Huber loss and a modified weight function, performs well in cases with outliers. It does not break down until more than half of the observations are arbitrary outliers. 

In a convex clustering problem, the input is a data matrix *X* of dimention *n*&times;*p*, with *n* is the number of samples and *p* the number of features. Convex clustering algorithm estimates a centroid matrix *U* of the same size as *X*, with *X<sub>i</sub>* and *X<sub>j</sub>* being in the same cluster if and only if *U<sub>i</sub>=U<sub>j</sub>*. 

This proposed algorithm is a modified version of that proposed by [Chi and Lange (2015)](https://arxiv.org/abs/1304.0499), implemented in package `cvxclustr`. The two packages are compared in [this paper](https://arxiv.org/abs/1906.09581v2) and can be reproduced. 

### Installation

Install `Rcvxclustr` from GitHub: 
```r
install.packages("devtools")
library(devtools)
devtools::install_github("JhZhang-1999/Rcvxclustr")
library(Rcvxclustr)
```

Helper functions can be accessed by typing `help(function name)` in R command. 

### Functions

Two main functions are implemented, and other functions in the package are dependencies. 

- `robust_weights`: This function implements the new weight function applied in fusion penalty, proposed in [the paper](https://arxiv.org/abs/1906.09581v2). It returns the *n &times;n* weight matrix for data *X* with *n* rows (number of samples). 

- `robustcvxclust`: This function implements the alternating direction method of multipliers algorithm for solving the convex clustering objective function. It returns the matrix of centriod differences, i.e. matrix *V*, with *U<sub>i</sub>-U<sub>i'</sub>=V<sub>ii'</sub>* for all *i<i'*. 

- `create_adjacency_matrix`: This function


### Examples

In the simulation, there are *n=20* observations that belong to two distinct non-overlapping clusters. The data matrix *X* is generated according to the model *X<sub>i</sub>=U<sub>1</sub>+&epsilon;<sub>i</sub>* when *i* belongs to the first cluster, and *X<sub>i</sub>=U<sub>2</sub>+&epsilon;<sub>i</sub>* when *i* belongs to the second cluster. Here *U<sub>1</sub>* and *U<sub>2</sub>* subject to different multivariate Gaussian distribution with identical covariance matrix and different means. 











