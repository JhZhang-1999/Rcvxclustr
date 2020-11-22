## Rcvxcluster

**Robust Convex Clustering Algorithm Implemented in R**

***

### Description

This package implements the robust convex clustering problem proposed by [Liu, Sun and Tan (2019)](https://arxiv.org/abs/1906.09581v2). Classical approaches towards convex clustering solves a convex optimization problem with the cost function being a squared loss plus a fusion penalty that encourages the estimated centroids for observations in the sames cluster to be identical. Those approaches are not robust to arbitrary outliers, and when data are contaminated, they fail to identify the correct cluster relationships. This proposed robust convex clustering algorithm, applying Huber loss and a modified weight function, performs well in cases with outliers. It does not break down until more than half of the observations are arbitrary outliers. 





