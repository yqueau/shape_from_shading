# Codes_Shape_from_shading
Matlab codes for shape-from-shading. Several variants are implemented: 
- ADMM-based variational shape from shading with general camera (orthographic or perspective) and lighting (spherical harmonics), see [1] 
- Lax-Friedriechs solving of the eikonal case (orthographic camera, frontal directional lighting), cf. Equation (8) in [2]  
- Semi-Lagrangian solver for the eikonal case (orthographic camera, frontal directional lighting), see [3] 
- Semi-Lagrangian solver for the perspective eikonal case (perspective camera, frontal directional lighting), see [4]

## Introduction

These codes can be used to solve the shape-from-shading (SfS) problem (estimate shape, given a single image). Main features:
- possibility to add a shape prior in order to guide the solution (useful for instance in RGB-D sensing)
- minimal surface regularization to smooth out the residual noise
- handles second-order spherical harmonics lighting
- handles orthographic or perspective camera
- handles grey or RGB images
 
Note: the classic eikonal SfS can also be achieved as a special case.


## Demos

The following two demo files accompanying [1] are provided: 

- `demo_1_lena_eikonal.m` : classic SfS (greylevel image, orthographic camera, frontal lighting) applied to the standard Lena image

- `demo_2_vase_SH2.m` : refinement of the depth map obtained with a RGB-D sensor. Source of the dataset: https://github.com/pengsongyou/SRmeetsPS

The three other demo files illustrate alternative PDE-based methods based semi-Lagrangian schemes, when lighting is frontal and directional, see [2,3,4] for details.   

## Contents

The main fuctions for ref [1] are in the Toolbox/ folder:
- `generic_sfs.m`: main SfS code
- `theta_fun.m`: cost function with respect to the surface gradient values
- `estimate_lighting.m`: can be used to estimate spherical harmonics lighting, given an image and a shape estimate
- `make_gradient.m`: finite differences stencils on a non-rectangular grid
- `export_obj2.m`: to produce a .obj file readable with meshlab


## Dependencies

- minFunc: need first be compiled: go to Toolbox/minFunc and run mexAll.m script 
(source: https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) 

- CMG (recommended for faster results): http://www.cs.cmu.edu/~jkoutis/cmg.html 

## References

[1] "A Variational Approach to Shape-from-shading Under Natural Illumination", Y. Quéau et al., EMMCVPR 2017

[2] "A comprehensive introduction to photometric 3D-reconstruction", J.-D. Durou et al., 2020

[3] "An algorithm for the global solution of the Shape-fromShading model", M. Falcone and M. Sagona, ICIAP 1997 

[4] "Some remarks on perspective shape-from-shading models", E. Cristiani et al., SSVM 2007
