# VBmixtures

Model-based clustering by Dirichlet Process or Finite Mixture models. Mean Field approximation is utilised for inference.

Four types of model are coded: Gaussian mixtures, Beta, Binomial and Poisson. Each section provides:
- simulator of mixture data
- initialisation of the variational algorithm
- the main variational algorithm

The discriminative accuracy function and the forward selection algorithm to compute the accuracy of the feature sets are also given for the Beta model. These functions can be modified accordingly for each model type by changing the likelihood.





