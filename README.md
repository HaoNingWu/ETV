# ETV
This is a collection of two MATLAB demos of Algorithms 1 (image denoising) and 2 (image reconstruction) in the paper "Enhanced total variation minimization for stable image reconstruction" by Congpei An, Hao-Ning Wu, and Xiaoming Yuan.

# Demo 1: Image denoising
The demo, demo_denoising.m, reproduces Figure 1 in our paper, examining the denoising ability of the enhanced TV model. 
Two functions are involved: 
denoiseTV.m: Solving the TV denoising model by the split Bregman
denoiseETV.m: Solving the enhanced TV denoising model by DCA+ADMM (Algorithm 1 in our paper)

# Demo 2: Image reconstruction
The demo, demo_reconstruction.m, reproduces the third column of Figure 5 in our paper, showing the reconstruction ability of the enhanced TV model. 
A lightweight demo, demo_reconstruction_lightweight.m, is also provided, which only reconstructs a 64-by-64 image.
Four functions are involved:
MRITV.m: Solving the TV reconstruction model by the split Bregman
MRIL12.m: Solving the weighted anisotropic and isotropic TV reconstruction model by DCA+split Bregman
MRIETV.m: Solving the anisotropic enhanced TV reconstruction model by DCA+ADMM (Algorithm 2 in our paper)
MRIETVisotropic.m: Solving the isotropic enhanced TV reconstruction model by DCA+ADMM

# When the measurements are noisy...
If the measurements are noisy with noise level $\tau$, our Algorithm 2 requires a projection onto the ball centered at 0 with radius $\tau$. 

To implement such a projection, the function project_L2.m in [The Proximity Operator Repository](http://proximity-operator.net/) is needed.
