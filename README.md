# sgdGMF
An R package for efficient estimation of generalized matrix factorization (GMF) models [[1,2,3]](#1,#2,#3).
The package implements the adaptive stochastic gradient descent with block- and coordinate-wise sub-sampling strategies proposed in [[4]](#4).
Additionally, sgdGMF implements the alternated iterative re-weighted least squares [[1,3]](#1,#3) and diagonal-Hessian quasi-Newton [[1]](#1) algorithms.

## References
<a id="1">[1]</a>
Collins, M., Dasgupta, S., Schapire, R.E. (2001).
A generalization of principal components analysis to the exponential family.
Advances in neural information processing systems, 14.

<a id="2">[2]</a>
Kidzinski, L., Hui, F.K.C., Warton, D.I., Hastie, T.J. (2022).
Generalized Matrix Factorization: efficient algorithms for fitting generalized linear latent variable models to large data arrays.
Journal of Machine Learning Research, 23(291): 1--29.

<a id="3">[3]</a>
Wang, L., Carvalho, L. (2023).
Deviance matrix factorization.
Electronic Journal of Statistics, 17(2): 3762--3810.

<a id="4">[4]</a>
Castiglione, C., Segers, A., Clement, L, Risso, D. (2024).
Stochastic gradient descent estimation of generalized matrix factorization models with application to single-cell RNA sequencing data.
arXiv preprint: arXiv:2412.20509.

