%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of viscosity solution to perspective shape-from-shading under
% vertical lighting
%
% Based on the paper: 
% E. Cristiani, M. Falcone and A. Seghini
% "Some remarks on perspective shape-from-shading models"
% In Proc. SSVM 2007, pp. 276-287
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath('Toolbox/')
addpath('Data/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('bunny.mat');
scale_factor = 1; % Set to <1 to downscale data for faster result
u = -imresize(u,scale_factor); % Ground truth depth map
mask = u>-Inf; % Pixels to be reconstructed
mask_eroded = imerode(mask,strel('square',3)); % Eroded mask
boundary = mask-mask_eroded; % Boundary of the eroded mask 
K = K*scale_factor; % Intrinsics matrix

% Simulate image using finite differences and ground truth depth map
[y,x] = meshgrid(1:size(mask,2),1:size(mask,1));
x = x-K(2,3);
y = y-K(1,3);
problem.I = min(0.98,max(0.02,simulate_image(x,y,u,-K(1,1),mask,[0;0;1],1)));

% Load image
%~ problem.I = min(0.98,max(0.02,imresize(I,scale_factor)));

% Load intrinsics
problem.f = -K(1,1); % Focal length
problem.x0 = K(2,3); % x coordinates of the principal point (x is "bottom")
problem.y0 = K(1,3); % y coordinates of the principal point (y is "right")

% Create mask: here we force u=ground truth on the boundaries
% and Create g (force u = g wherever mask=0)
problem.mask = mask; % Take the whole mask
problem.g = min(u(problem.mask>0))*ones(size(problem.I)); % Fix depth at boundary to max ground truth value

%~ problem.mask = mask_eroded; % Take the eroded mask
%~ problem.g = max(u(problem.mask>0))*ones(size(problem.I)); % Fix depth at max ground truth value
%~ problem.g(boundary>0) = u(boundary>0); % On the boundary fix ground truth values

figure(145)

subplot(2,2,1)
imshow(uint8(255*problem.I))
axis image
axis off
title('$$I$$','Interpreter','Latex','Fontsize',18);

subplot(2,2,3)
imagesc(problem.g)
colorbar
axis image
axis off
title('$$g$$','Interpreter','Latex','Fontsize',18);

subplot(2,2,4)
imagesc(problem.mask)
colorbar
axis image
axis off
title('Mask','Interpreter','Latex','Fontsize',18);

drawnow



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters.delta = 1/max(size(I,1),size(I,2)); % Grid size
parameters.h = 0.999*parameters.delta; % Stepsize
parameters.xi = 1e-8; % Stopping criterion
parameters.maxit = 1e6; % Stopping criterion
parameters.N_theta = 16; % Number of angles for discretizing the unit ball
parameters.display = 2; % 0: no display, 1: numerics, 2: 3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u,x,y] = eikonal_sl_perspective(problem,parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u(problem.mask==0) = NaN;

figure('Name','Simulated image')
surf(x,y,u,ones(size(u)));
shading flat
colormap gray
axis equal
%~ axis off
camproj('perspective')
campos([0 0 0])
camtarget([0 0 -1])
camup([-1 0 0])
camva = rad2deg(2*atan(0.5*size(problem.I,2)/(problem.f)));
material([0 1 0]);
camlight('headlight','infinite')
camlight('headlight','infinite')
xlabel('$$x(X,Y)$$','Interpreter','Latex','Fontsize',18)
ylabel('$$y(X,Y)$$','Interpreter','Latex','Fontsize',18)
print -dpng frontal_view

figure('Name','3D-reconstruction')
surf(x,y,u,ones(size(u)));
shading flat
colormap gray
axis equal
camproj('perspective')
campos([0 0 0])
camtarget([0 0 -1])
camup([-1 0 0])
camva = rad2deg(2*atan(0.5*size(problem.I,2)/(problem.f)));
material([0 1 0]);
camlight('headlight','infinite')
camlight('headlight','infinite')
campos([2 0 0])
xtarget = round(problem.x0);
ytarget = round(problem.y0);
camtarget([x(xtarget,ytarget) y(xtarget,ytarget) u(xtarget,ytarget)])
xlabel('$$x(X,Y)$$','Interpreter','Latex','Fontsize',18)
ylabel('$$y(X,Y)$$','Interpreter','Latex','Fontsize',18)
zlabel('$$u(X,Y)$$','Interpreter','Latex','Fontsize',18)
%~ title('3D-reconstruction','Interpreter','Latex','Fontsize',18,'Location','West')
print -dpng lateral_view

