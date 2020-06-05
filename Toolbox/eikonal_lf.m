function u = eikonal_sl(problem,parameters)
%EIKONAL_SL Lax-Friedrichs solver for the Eikonal equation |\nabla u| = f
% 
% Required inputs:
% 	problem.f (nrows x ncols): the 2D velocity map
%   problem.g (nrows x ncols): the boundary condition (used only on mask)
%   problem.mask (nrows x ncols): mask of interior points (u is set to problem.g wherever mask==0)
%   
%   parameters.delta: Grid size (e.g. 1/size(f,2))
%   parameters.h: viscosity (e.g. 1)
%   parameters.epsilon: truncation for f (e.g. = 0.2)
%   parameters.xi: stopping criterion (e.g. 1e-8)
%   parameters.maxit: stopping criterion (e.g. 1e6) 

	%% Truncate f to avoid ambiguities
	f = max(problem.f,parameters.epsilon);

	%% Initialization
	u = problem.g; 
	% Residual
	residual = Inf;
	% Iteration
	it = 0;

	%% A few shortcuts...
	% Indices of pixels inside the mask
	imask = find(problem.mask>0);
	% Grid coordinates
	[y,x] = meshgrid(1:size(problem.f,2),1:size(problem.f,1));
	x = parameters.delta*x;
	y = parameters.delta*y;

	%% Main loop
	while(residual>parameters.xi & it<parameters.maxit)
		u_old = u;
		it = it+1;
		
		% Average of u, over the four neighbors
		u_average = 0.5*(u([2:end end],:)+u([1 1:end-1],:)...
					    +u(:,[2:end end])+u(:,[1 1:end-1]))/parameters.delta;

		% Norm of Grad(u), with centered finite differences
		uy = 0.5*(u([2:end end],:)-u([1 1:end-1],:))/parameters.delta;
		ux = 0.5*(u(:,[2:end end])-u(:,[1 1:end-1]))/parameters.delta;
		norm_grad = sqrt(ux.^2+uy.^2); 
		
		% Lax-Friedrichs update
		u(imask) = 0.5*parameters.delta*(parameters.h*u_average(imask)-(norm_grad(imask)-problem.f(imask)))/parameters.h;
		
		% Simulate image based on current estimate
		I_simulated = simulate_image_ortho(u,problem.mask,[0;0;1],parameters.delta);
		MAD = mean(abs(problem.I(imask)-I_simulated(imask)));
		RMSE = sqrt(mean((problem.I(imask)-I_simulated(imask)).^2));
		
		% Evaluate residual
		residual = max(u(imask)-u_old(imask));
		
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
