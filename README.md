# normalized-cuts

Here is a java implementation of the Normalized Cuts image segmentation algorithm as proposed in this [paper](https://people.eecs.berkeley.edu/~malik/papers/SM-ncut.pdf) by Jianbo Shi and Jitendra Malik. This is a recursive spatial clustering algorithm, which takes influence from spectral clustering.

## Structure

Package normCuts stores an abstract java class NormalizedCuts.java which stores all logic shared between any implementation of the algorithm. There are three methods that must be implemented when implementing this abstract class: 

	- setSimilarityMatrix()
	- getFeatureFactor()
	- getWeight()

These methods allow an implementation of the algorithm to be created for your specific need. For example we have an implementation to be used for grayscale images, NormalizedCutsGray. These methods are where we can define the similarity between any two points in our dataset, ie pixels.

## Lanczos Algorithm

We end up with some massive matrices in this algorithm. For example if an image is 50x100 pixels. Then the similarity matrix in the first recursive call will be of size (50x100) x (50x100). This would result in 25000000 double entries in our matrix. We are required to eigendecompose a matrix of the same size, and for large image sizes this becomes infeasible. Hence, we include an implementation of the Lanczos algorithm, which is used for calculating the k most extreme eigenvalues and eigenvectors, this code is found in normCuts.LanczosSolver. The wikipedia page for the Lanczos algorithm can be found [here](https://en.wikipedia.org/wiki/Lanczos_algorithm). It exploits increased eigencdecomposition efficiency for tridiagonal matrices.

## Dependencies

We are dependent on the apache commons math library. 

	- commons-math3-3.6.1-tools.jar
	- commons-math3-3.6.1.jar
	
We also have some methods that use java swing to display a simple GUI of the clustering at each iteration of the algorithm. I am planning to change this so that we will use imageJ to do this, so in the future we will be dependent on that.

## Limitations

We are limited in the maximum size of image. The apache commons math sparse matrix class OpenMapRealMatrix limits us to a maximum size, meaning the max image size we can perform the algorithm on is 215x215. I would like to create a sparse matrix class of my own which is optimised to symmetric matrices as in our case which would allow our algorithm to run on larger images. 

