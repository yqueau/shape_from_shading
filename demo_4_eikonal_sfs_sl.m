%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo of viscosity solution to orthographic shape-from-shading under
% vertical lighting
%
% Based on the paper: 
% M. Falcone and M. Sagona
% "An algorithm for the global solution of the Shape-fromShading model" 
% in Proc. ICIAP 1997, pp. 596–603.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath('Toolbox/')
addpath('Data/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load image 
I = double(imresize(rgb2gray(imread('lena.png')),[128 128]))./255;
I = min(0.99,max(I,0.01));

% Second member of eikonal PDE
problem.I = I;
problem.f = sqrt(1./(I.^2)-1);

% Create g (force u = g wherever mask=0)
problem.g = zeros(size(I));

% Create mask: here we force u=0 on the boundaries, and on 5 points inside
problem.mask = ones(size(I));
problem.mask(1,:) = 0;
problem.mask(end,:) = 0;
problem.mask(:,1) = 0;
problem.mask(:,end) = 0;
problem.mask(floor(0.5*size(I,1)),floor(0.5*size(I,2))) = 0;
problem.mask(floor(0.25*size(I,1)),floor(0.25*size(I,2))) = 0;
problem.mask(floor(0.25*size(I,1)),floor(0.75*size(I,2))) = 0;
problem.mask(floor(0.75*size(I,1)),floor(0.75*size(I,2))) = 0;
problem.mask(floor(0.75*size(I,1)),floor(0.25*size(I,2))) = 0;

figure(145)

subplot(2,2,1)
imshow(uint8(255*I))
axis image
axis off
title('$$I$$','Interpreter','Latex','Fontsize',18);

subplot(2,2,2)
imagesc(problem.f)
colorbar
axis image
axis off
title('$$f$$','Interpreter','Latex','Fontsize',18);

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
parameters.h = parameters.delta./max(1./problem.f(:)); % Stepsize
parameters.epsilon = 0.2; % For truncation of f
parameters.xi = 1e-8; % Stopping criterion
parameters.maxit = 1e6; % Stopping criterion
parameters.N_theta = 16; % Number of angles for discretizing the unit ball
parameters.display = 1; % 0: no display, 1: numerics, 2: 3D


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = eikonal_sl(problem,parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,x] = meshgrid(1:size(problem.f,2),1:size(problem.f,1));
x = parameters.delta*x;
y = parameters.delta*y;

figure(146)
surf(x,y,u,ones(size(u)));
shading flat
colormap gray
axis equal
axis off
material([0 1 0]);
view(90,90);
camlight('headlight','infinite')
camlight('headlight','infinite')
title('Simulated image','Interpreter','Latex','Fontsize',18)

figure(147)
surf(x,y,u,ones(size(u)));
shading flat
colormap gray
axis equal
view(90,90);
camlight('headlight','infinite')
camlight('headlight','infinite')
view(70,30);
xlabel('$$X$$','Interpreter','Latex','Fontsize',18)
ylabel('$$Y$$','Interpreter','Latex','Fontsize',18)
zlabel('$$u(X,Y)$$','Interpreter','Latex','Fontsize',18)
title('3D-reconstruction','Interpreter','Latex','Fontsize',18)
