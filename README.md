# Ksparse Clustering
### Clustering with removal of noisy features using alternating minimization over the l1-ball
- Cyprien Gilet, Michel Barlaud, Jean-Baptiste Caillau, Marie Deprez

We provide here our unsupervised clustering method with removal of noisy features in high dimensional space in Matlab. 
The problem is to estimate both labels and a sparse projection matrix of weights. 
To address this combinatorial non-convex problem maintaining a strict control on the sparsity of this matrix of weights, we propose an alternating minimization of the Frobenius norm criterion.
We provide a new efficient algorithm named k-sparse which alternates k-means with projection-gradient minimization. 
The projection-gradient step is a method of splitting type, with exact projection on the $\ell^1$ ball to promote sparsity.

### This folder contains:
- The script "script_ksparse_simu_real_data.m" we used for all our experiments.
- A script "script_ksparse_datareal_silh_eta.m" which provides an example of how to efficiently select the best parameter $\eta$.
- A script "script_ksparse_general.m" which can be easily used for any other databases.
- A subfolder "data_prepared" containing some databases to experiment.
- A subfolder "functions" containing all the functions.
