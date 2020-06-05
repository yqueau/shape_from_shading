    function [u,x,y] = eikonal_sl_perspective(problem,parameters)
%EIKONAL_SL_PERSPECTIVE Semi-Lagrangian solver for the Eikonal-like equation |\nabla u| = f abs(u + p.\nabla u)
% with f = sqrt((1/I^2-1)/focal^2)
% 
% Required inputs:
% 	problem.I (nrows x ncols): the 2D image
%   problem.g (nrows x ncols): the boundary condition (used only on mask)
%   problem.mask (nrows x ncols): mask of interior points (u is set to problem.g wherever mask==0)
%   problem.f: focal length
%   problem.x0: x coordinates of the principal point
%   problem.y0: y coordinates of the principal point
%   
%   parameters.delta: Grid size (e.g. 1.size(f,2))
%   parameters.h: stepsize for descent (e.g. delta./max(1./f(:)))
%   parameters.xi: stopping criterion (e.g. 1e-8)
%   parameters.maxit: stopping criterion (e.g. 1e6) 
%   parameters.N_theta: number of angles for discretizing the unit ball (e.g. 8)
%   parameters.display: 0: no display, 1: numerics, 2: show live 3D-reconstruction 

	problem.g = problem.g*parameters.delta;
	problem.f = problem.f*parameters.delta;
	problem.x0 = problem.x0*parameters.delta;
	problem.y0 = problem.y0*parameters.delta;
	problem.I_tilde = sqrt((1./(problem.I.^2)-1)./(problem.f^2));

	%% Initialization
	% Shape
	u = problem.g;
	% Residual
	residual = Inf;
	% Iteration
	it = 0;

	%% A few shortcuts...
	% Indices of pixels inside the mask
	imask = find(problem.mask>0);
	% Points on the unit ball
	a = zeros(2,parameters.N_theta);
	for i = 1:parameters.N_theta
		theta_i = i*2*pi/parameters.N_theta;
		a(1,i) = cos(theta_i);
		a(2,i) = sin(theta_i);
	end
	% Allocate memory for the directional derivative search
	u_dir = repmat(u,[1 1 parameters.N_theta]);
	% Grid coordinates
	[y,x] = meshgrid(1:size(problem.I_tilde,2),1:size(problem.I_tilde,1));
	x = parameters.delta*x;
	y = parameters.delta*y;
	
	%% Main loop
	while(residual>parameters.xi & it<parameters.maxit)
		u_old = u;
		it = it+1;
		
		% Update v 
		for i = 1:parameters.N_theta
			inc_xi = (x-problem.x0)+a(1,i)./problem.I_tilde;
			inc_yi = (y-problem.y0)+a(2,i)./problem.I_tilde;
			x_i = min(max(x(:)),max(min(x(:)),x+parameters.h*inc_xi));
			y_i = min(max(y(:)),max(min(y(:)),y+parameters.h*inc_yi));
			u_i = interp2(y,x,u,y_i,x_i); 
			u_i(problem.mask==0) = u_dir(problem.mask==0); 
			u_dir(:,:,i) = u_i; 
		end
		min_ui = min(u_dir,[],3);	
		u(imask) = min_ui(imask)./(parameters.h+1);
		
		% Simulate image based on current estimate
		I_simulated = simulate_image((x-problem.x0),(y-problem.y0),u,problem.f,problem.mask,[0;0;1],parameters.delta);
		MAD = mean(abs(problem.I(imask)-I_simulated(imask)));
		
		% Evaluate residual
		residual = max(abs(u(imask)-u_old(imask)));
		
		% Display some information
		if(parameters.display==1)
			disp(sprintf('It %d - residual : %.9f - MAD : %.4f',it,residual,MAD));
		end
		if(parameters.display==2)
			disp(sprintf('It %d - residual : %.9f - MAD : %.4f',it,residual,MAD));
			
			xx = ((u.*(x-problem.x0)./problem.f))./parameters.delta;
			yy = ((u.*(y-problem.y0)./problem.f))./parameters.delta;
			zz = u./parameters.delta;
			zz(problem.mask==0) = NaN;
			
		
			
			figure(146)
			surf(xx,yy,zz,ones(size(zz)));
			shading flat
			colormap gray
			axis equal
			axis off
			camproj('perspective')
			campos([0 0 0])
			camtarget([0 0 -1])
			camup([-1 0 0])
			camva = rad2deg(2*atan(0.5*size(problem.I,2)/(problem.f./parameters.delta)));
			material([0 1 0]);
			camlight('headlight','infinite')
			camlight('headlight','infinite')
			title('Simulated image','Interpreter','Latex','Fontsize',18)

			figure(147)
			surf(xx,yy,zz,ones(size(zz)));
			shading flat
			colormap gray
			axis equal
			camproj('perspective')
			campos([0 0 0])
			camtarget([0 0 -1])
			camup([-1 0 0])
			camva = rad2deg(2*atan(0.5*size(problem.I,2)/(problem.f./parameters.delta)));
			material([0 1 0]);
			camlight('headlight','infinite')
			camlight('headlight','infinite')
			campos([2 0 0])
			xtarget = round(problem.x0/parameters.delta);
			ytarget = round(problem.y0/parameters.delta);
			camtarget([xx(xtarget,ytarget) yy(xtarget,ytarget) zz(xtarget,ytarget)])
			xlabel('$$x(X,Y)$$','Interpreter','Latex','Fontsize',18)
			ylabel('$$y(X,Y)$$','Interpreter','Latex','Fontsize',18)
			zlabel('$$u(X,Y)$$','Interpreter','Latex','Fontsize',18)
			title('3D-reconstruction','Interpreter','Latex','Fontsize',18)
			
			figure(148)
			imagesc(problem.I,[0 1])
			colormap gray
			axis image
			title('Input')
			figure(149)
			imagesc(I_simulated,[0 1])
			colormap gray
			axis image
			title('Simulation')
			figure(150)
			imagesc(abs(problem.I-I_simulated))
			axis image
			colorbar
			title('Absolute diff')
			
		end 
	end
	x = (u.*(x-problem.x0)./problem.f)./parameters.delta;
	y = (u.*(y-problem.y0)./problem.f)./parameters.delta;
	u = u./parameters.delta;
end
