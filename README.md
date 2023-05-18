ETV is a Github repository that provides MATLAB demos of two algorithms described in the paper "Enhanced total variation minimization for stable image reconstruction" by Congpei An, Hao-Ning Wu, and Xiaoming Yuan. The software demonstrates Algorithm 1 for image denoising and Algorithm 2 for image reconstruction in this paper, which can be accessed at this [link](https://iopscience.iop.org/article/10.1088/1361-6420/acd4e1).

# Demo 1: Image Denoising
The demo_denoising.m script reproduces Figure 1 from the paper, showcasing the denoising capability of the enhanced TV model. It utilizes two functions:
* denoiseTV.m: This function solves the TV denoising model using the split Bregman method.
* denoiseETV.m: This function solves the enhanced TV denoising model using the DCA+ADMM approach (Algorithm 1 in the paper).

# Demo 2: Image Reconstruction
The demo_reconstruction.m script reproduces the third column of Figure 5 in the paper, illustrating the reconstruction capability of the enhanced TV model. Additionally, a lightweight version of the demo, demo_reconstruction_lightweight.m, is provided, which focuses on reconstructing a 64-by-64 image.

Four functions are involved in the reconstruction demos:
* MRITV.m: This function solves the TV reconstruction model using the split Bregman method.
* MRIL12.m: This function solves the weighted anisotropic and isotropic TV reconstruction model using the DCA+split Bregman approach.
* MRIETV.m: This function solves the anisotropic enhanced TV reconstruction model using the DCA+ADMM approach (Algorithm 2 in the paper).
* MRIETVisotropic.m: This function solves the isotropic enhanced TV reconstruction model using the DCA+ADMM approach.

# Handling Noisy Measurements
In situations where the measurements are noisy, with a noise level denoted as $\tau$, our Algorithm 2 requires a projection onto a ball centered at 0 with a radius of $\tau$. To implement this projection, the function project_L2.m from [The Proximity Operator Repository](http://proximity-operator.net/index.html) needs to be utilized.
