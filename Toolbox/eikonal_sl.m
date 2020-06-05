function u = eikonal_sl(problem,parameters)
%EIKONAL_SL Semi-Lagrangian solver for the Eikonal equation |\nabla u| = f
% 
% Required inputs:
% 	problem.f (nrows x ncols): the 2D velocity map
%   problem.g (nrows x ncols): the boundary condition (used only on mask)
%   problem.mask (nrows x ncols): mask of interior points (u is set to problem.g wherever mask==0)
%   
%   parameters.delta: Grid size (e.g. 1.size(f,2))
%   parameters.h: stepsize for descent (e.g. delta./max(1./f(:)))
%   parameters.epsilon: truncation for f (e.g. = 0.2)
%   parameters.xi: stopping criterion (e.g. 1e-8)
%   parameters.maxit: stopping criterion (e.g. 1e6) 
%   parameters.N_theta: number of angles for discretizing the unit ball (e.g. 8)
%   parameters.display: 0: no display, 1: numerics, 2: show live 3D-reconstruction 

	%% Truncate f to avoid ambiguities
	f = max(problem.f,parameters.epsilon);

	%% Initialization
	% Shape
	v = 1-exp(-problem.g);
	v(problem.mask>0) = 0;
	u = -log(1-v);
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
	v_dir = repmat(v,[1 1 parameters.N_theta]);
	% Grid coordinates
	[y,x] = meshgrid(1:size(problem.f,2),1:size(problem.f,1));
	x = parameters.delta*x;
	y = parameters.delta*y;

	%% Main loop
	while(residual>parameters.xi & it<parameters.maxit)
		u_old = u;
		v_old = v;
		it = it+1;
		
		% Update v 
		for i = 1:parameters.N_theta
			x_i = min(max(x(:)),max(min(x(:)),x+parameters.h*a(1,i)./problem.f));
			y_i = min(max(y(:)),max(min(y(:)),y+parameters.h*a(2,i)./problem.f));
			v_i = interp2(y,x,v,y_i,x_i); 
			v_i(problem.mask==0) = v_dir(problem.mask==0); 
			v_dir(:,:,i) = v_i; 
		end
		min_vi = min(v_dir,[],3);	
		v(imask) = 1+exp(-parameters.h)*(min_vi(imask)-1);
		
		% Update u
		u = -log(1-v);
		
		% Simulate image based on current estimate
		I_simulated = simulate_image_ortho(u,problem.mask,[0;0;1],parameters.delta);
		MAD = mean(abs(problem.I(imask)-I_simulated(imask)));
		RMSE = sqrt(mean((problem.I(imask)-I_simulated(imask)).^2));
		
		% Evaluate residual
		residual = max(v(imask)-v_old(imask));
		
		% Display some information
		if(parameters.display==1)
			disp(sprintf('It %d - residual : %.9f - MAD : %.4f - RMSE : %.4f',it,residual,MAD,RMSE));
		end
		if(parameters.display==2)
			disp(sprintf('It %d - residual : %.9f - MAD : %.4f - RMSE : %.4f',it,residual,MAD,RMSE));
			
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
end
