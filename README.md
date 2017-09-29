# Codes_Shape_from_shading
Matlab codes for ADMM-based variational shape from shading with spherical harmonics lighting 

## Introduction

These codes can be used to solve the shape-from-shading (SfS) problem (estimate shape, given a single image). Main features:
- possibility to add a shape prior in order to guide the solution (useful for instance in RGB-D sensing)
- minimal surface regularization to smooth out the residual noise
- handles second-order spherical harmonics lighting
- handles orthographic or perspective camera
- handles grey or RGB images
 
Note: the classic eikonal SfS can also be achieved as a special case.


## Demos

The following demo files are provided: 

- `demo_1_lena_eikonal.m` : classic SfS (greylevel image, orthographic camera, frontal lighting) applied to the standard Lena image

- `demo_2_vase_SH2.m` : refinement of the depth map obtained with a RGB-D sensor. Source of the dataset: https://github.com/pengsongyou/SRmeetsPS


## Contents

The main fuctions are in the Toolbox/ folder:
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
