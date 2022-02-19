# **Python unconstraint optimization package :**

This package contains the implementation of all the introduced algorithms in the course of unconstraint optimization for the ``1^{st}`` year AI engineering student at UM5-ENSIAS, with in addition, the possibility to plot the results of each optimization.

This package requires the following Python 3 libraries :
- numpy
- time
- math
- numdifftools
- scipy
- matplotlib
- mpl_toolkits 
- copy

The module `matrix.py` contains some elementary operations on matrices :
- sum
- product
- transposed
- rows swap
- rows dilatation
- rows transvection
- zeros matrix generator
- identity matrix generator
- printing a matrix

The module `gauss.py` contains the implementation of the operations :
- Gauss-Jordan pivot
- Inverse of a matrix using Gauss
- Solving a linear system with Gauss pivot

The module `deco.py` containing the code of the decomposition algorithms :
- Lower-Upper decomposition
- Cholesky's decomposition
and also, 2 functions to solve a linear system with the 2 above decompositions.

The module `gradient.py` contains the code of 3 optimization algorithms for multivariable functions :
- Gradient descent
- Conjugate gradient
- Newton gradient

As for the monovariable unimodal functions, the subfolder unidim contains two modules (`elimination.py` and `interpolation.py`) with 10 different methods of optimization :
- unconstrainted search with fixed steps
- unconstrainted search with accelerated steps
- brute force 
- binary search
- interval halving method
- Fibonacci method
- Golden section method
- Newton-Rapson method
- Quasi-Newton method
- Secant method

This package is developped with 3 main programs : `main.py`, `main2.py` and `test_compare.py`. They are built-in to help using some of the methods on random functions and data implemented in the file `functions.py`.
