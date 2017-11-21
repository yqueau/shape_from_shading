function [z,N,XYZ,N_SH] = generic_sfs(data,params,options);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Preprocessings
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Subsampling
	if(options.ratio>1)
		data.I = data.I(1:options.ratio:end,1:options.ratio:end,:);
		data.mask = data.mask(1:options.ratio:end,1:options.ratio:end);
		data.mask_z0 = data.mask_z0(1:options.ratio:end,1:options.ratio:end);
		data.rho = data.rho(1:options.ratio:end,1:options.ratio:end,:);
		data.z0 = data.z0(1:options.ratio:end,1:options.ratio:end);
		data.z_init = data.z_init(1:options.ratio:end,1:options.ratio:end);
		data.K(1:2,:) = data.K(1:2,:)./options.ratio;
	end
	clear data.ratio;
	
	% Initialization of depth
	z = data.z_init;
	
	% Masked pixels
	imask = find(data.mask>0);
	imask_z0 = find(data.mask_z0>0);
	[foo,imask_z0] = ismember(imask_z0,imask);
	imask_z0 = imask_z0(find(foo>0));
	mat_z0 = speye(length(imask));
	mat_z0 = mat_z0(imask_z0,:);
	clear data.mask_z0 data.mask

	% Auxiliary variables
	npix = length(imask);
	npix_z0 = length(imask_z0);
	[nrows,ncols,nchannels] = size(data.I);
	if(data.K(1,3)>0)
		[xx,yy] = meshgrid(1:ncols,1:nrows);
		xx = xx(imask);
		xx = xx-data.K(1,3);
		yy = yy(imask);
		yy = yy-data.K(2,3);
		z = log(z(imask));
		data.z0 = log(data.z0(imask));
	else
		xx = zeros(npix,1);
		yy = zeros(npix,1);
		data.K(1,1) = 1;
		data.K(2,2) = 1;
		data.K(1,3) = 0;
		z = z(imask);
		data.z0 = data.z0(imask);
	end
	G = make_gradient(data.mask); % Finite differences stencils
	Dx = G(1:2:end-1,:);
	Dy = G(2:2:end,:);
	clear G

	% Set BFGS options
	opts_minfunc = [];
	opts_minfunc.display = 'iter';
	opts_minfunc.MaxIter = options.maxit_bfgs;
	opts_minfunc.optTol = options.tolX_bfgs;
	opts_minfunc.progTol = options.tolFun_bfgs;
	opts_minfunc.DerivativeCheck = options.check_grad;
	opts_minfunc.numDiff = options.use_jac;
	clear options.maxit_bfgs options.tolX_bfgs options.tolFun_bfgs options.check_grad options.use_jac

	% Vectorize data
	data.I = reshape(data.I,[nrows * ncols,nchannels]);
	data.I = data.I(imask,:);
	data.rho = reshape(data.rho,[nrows * ncols,nchannels]);
	data.rho = data.rho(imask,:);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Initialization
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Initial gradient
	zx = Dx*z;
	zy = Dy*z;
	zx(isnan(zx)) = 0;
	zy(isnan(zy)) = 0;
	z(isnan(z)) = 1;

	% Initial augmented normals
	N = zeros(npix,9); 
	dz = max(eps,sqrt((data.K(1,1)*zx).^2+(data.K(2,2)*zy).^2+(-1-xx.*zx-yy.*zy).^2));
	N(:,1) = data.K(1,1)*zx./dz;
	N(:,2) = data.K(2,2)*zy./dz;
	N(:,3) = (-1-xx.*zx-yy.*zy)./dz;
	N(:,4) = 1;
	N(:,5) = N(:,1).*N(:,2);
	N(:,6) = N(:,1).*N(:,3);
	N(:,7) = N(:,2).*N(:,3);
	N(:,8) = N(:,1).^2-N(:,2).^2;
	N(:,9) = 3*N(:,3).^2-1;

	% Initial A and b fields
	A = zeros(npix,nchannels,2); 
	B = zeros(npix,nchannels,1);
	for ch = 1:nchannels
		A(:,ch,1) = data.rho(:,ch).*(data.K(1,1)*data.s(ch,1)-xx*data.s(ch,3))./dz;
		A(:,ch,2) = data.rho(:,ch).*(data.K(2,2)*data.s(ch,2)-yy*data.s(ch,3))./dz;

		B(:,ch) = data.I(:,ch)+data.rho(:,ch)*data.s(ch,3)./dz-data.rho(:,ch).*sum(bsxfun(@times,data.s(ch,4:9),N(:,4:9)),2);
	end

	% Initial dual variables for ADMM
	theta_p = zx;
	theta_q = zy;
	u_p = zeros(npix,1);
	u_q = zeros(npix,1);
	rk = NaN;
	sk = NaN;
	tab_primal = [norm(rk)];
	tab_dual = [norm(sk)];

	
	% Initial energy: shading + prior + smoothness
	energy = 0;
	% Shading term
	for ch = 1:nchannels
		energy = energy + 0.5*params.lambda*sum((A(:,ch,1).*zx+A(:,ch,2).*zy-B(:,ch)).^2);
	end
	% Prior term
	energy = energy + 0.5*params.mu*sum((z(imask_z0)-data.z0(imask_z0)).^2);
	% Smoothness term
	energy = energy+params.nu*sum(dz);
	tab_energy = [energy];
	disp(sprintf('It. 0 - energy : %.6f',energy));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Display
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if(options.display)

		% Fig 2 : shape
		figure(2)
		Ndx = zeros(size(data.mask));Ndx(imask) = N(:,1);
		Ndy = zeros(size(data.mask));Ndy(imask) = N(:,2);
		Ndz = zeros(size(data.mask));Ndz(imask) = -N(:,3);
		Nd = 0.5*(1+cat(3,Ndx,Ndy,Ndz));
		imshow(min(1,max(0,Nd)))
		axis equal
		axis off
		title('Normals','Interpreter','Latex','Fontsize',18)

		drawnow	
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Main loop
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for it = 1:options.maxit

		z_before = z;
		theta_before = [theta_p;theta_q];

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Nonlinear theta update
		theta = minFunc(@(theta)theta_fun(theta,zx,zy,u_p,u_q,data.I,data.s,data.rho,options.beta,xx,yy,data.K,imask,params.lambda,params.nu),theta_before,opts_minfunc);
		theta_p = theta(1:npix);
		theta_q = theta(npix+1:2*npix);


		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Linear z update
		mat_z = [];
		sec_z = [];
		% Prior part
		mat_z = [mat_z;sqrt(0.5*params.mu)*mat_z0];
		sec_z = [sec_z;sqrt(0.5*params.mu)*data.z0(imask_z0)];
		% ADMM part
		mat_z = [mat_z;sqrt(0.5*options.beta)*Dx;sqrt(0.5*options.beta)*Dy];
		sec_z = [sec_z;sqrt(0.5*options.beta)*(theta_p-u_p);sqrt(0.5*options.beta)*(theta_q-u_q)];
		% Normal equations
		sec_z = mat_z'*sec_z;
		mat_z = mat_z'*mat_z;	
		% Preconditioning
		if(strcmp(options.precond_pcg,'none'))
			M1 = [];
			M2 = [];
		elseif(strcmp(options.precond_pcg,'ichol'))
			M1 = ichol(mat_z);
			M2 = transpose(M1);
		elseif(strcmp(options.precond_pcg,'cmg'))
			M1 = cmg_sdd(mat_z);
			M2 = [];
		end
		% PCG
		z = pcg(mat_z,sec_z,options.tolFun_pcg,options.maxit_pcg,M1,M2,z_before);
		% Update gradient
		zx = Dx*z;
		zy = Dy*z;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Dual update
		u_p = u_p + zx - theta_p;
		u_q = u_q + zy - theta_q;

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% dz, N updates
		dz = max(eps,sqrt((data.K(1,1)*zx).^2+(data.K(2,2)*zy).^2+(-1-xx.*zx-yy.*zy).^2));
		N(:,1) = data.K(1,1)*zx./dz;
		N(:,2) = data.K(2,2)*zy./dz;
		N(:,3) = (-1-xx.*zx-yy.*zy)./dz;
		N(:,4) = 1;
		N(:,5) = N(:,1).*N(:,2);
		N(:,6) = N(:,1).*N(:,3);
		N(:,7) = N(:,2).*N(:,3);
		N(:,8) = N(:,1).^2-N(:,2).^2;
		N(:,9) = 3*N(:,3).^2-1;
	
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% A,B updates
		for ch = 1:nchannels
			A(:,ch,1) = data.rho(:,ch).*(data.K(1,1)*data.s(ch,1)-xx*data.s(ch,3))./dz;
			A(:,ch,2) = data.rho(:,ch).*(data.K(2,2)*data.s(ch,2)-yy*data.s(ch,3))./dz;
			B(:,ch) = data.I(:,ch)+data.rho(:,ch).*(data.s(ch,3))./dz-data.rho(:,ch).*sum(bsxfun(@times,data.s(ch,4:9),N(:,4:9)),2);
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% New energy
		energy = 0;
		% Shading term
		for ch = 1:nchannels
			energy = energy + 0.5*params.lambda*sum((A(:,ch,1).*zx+A(:,ch,2).*zy-B(:,ch)).^2);
		end
		% Prior term
		energy = energy + 0.5*params.mu*sum((z(imask_z0)-data.z0(imask_z0)).^2);
		% Smoothness term
		energy = energy+params.nu*sum(dz);
		tab_energy = [tab_energy,energy];
				
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Primal and dual residuals
		rk = [zx-theta_p;zy-theta_q];
		sk = [-options.beta*(theta_p-theta_before(1:npix));-options.beta*(theta_q-theta_before(npix+1:2*npix))];
		resPrim = norm(rk);
		resDual = norm(sk);
		relResPrim = resPrim./max(norm([zx;zy]),norm([theta_p;theta_q]));
		relResDual = resDual./norm([u_p;u_q]);
		tab_primal = [tab_primal,relResPrim];
		tab_dual = [tab_dual,relResDual];

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Display
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if(options.display)
			
			% Fig 2 : shape
			figure(2)
			Ndx = zeros(size(data.mask));Ndx(imask) = N(:,1);
			Ndy = zeros(size(data.mask));Ndy(imask) = N(:,2);
			Ndz = zeros(size(data.mask));Ndz(imask) = -N(:,3);
			Nd = 0.5*(1+cat(3,Ndx,Ndy,Ndz));
			imshow(min(1,max(0,Nd)))
			axis equal
			axis off
			title('Normals','Interpreter','Latex','Fontsize',18)


			% Fig 3 : energy
			figure(3)
			subplot(1,3,1)
			semilogy(0:it,tab_energy)
			title('Energy','Interpreter','Latex','Fontsize',18)
			subplot(1,3,2)
			semilogy(0:it,tab_primal)
			title('Relative primal residual','Interpreter','Latex','Fontsize',18)
			subplot(1,3,3)
			semilogy(0:it,tab_dual)
			title('Relative dual residual','Interpreter','Latex','Fontsize',18)
		
			drawnow	
		end	

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% CV Test
		resX = norm(z-z_before)./norm(z_before);
		resFun = abs(tab_energy(end)-tab_energy(end-1))./abs(tab_energy(end-1));
		disp(sprintf('It. %d - energy: %.6f - resF: %.6f - resX: %.6f - resP: %.6f - resD: %.6f - beta: %.6f',it,energy,resFun,resX,relResPrim,relResDual,options.beta));
		% Update penalty
		if(resPrim/resDual > options.tau)
			options.beta = options.eta*options.beta;
			u_p = u_p./options.eta;
			u_q = u_q./options.eta;
			dontstop = 1;
		elseif(resDual/resPrim > options.tau)
			options.beta = options.beta/options.eta;			
			u_p = u_p.*options.eta;
			u_q = u_q.*options.eta;
			dontstop = 1;
		else
			dontstop = 0;
		end
		% Test
		if((resFun<options.tolFun & it>options.minit) | (relResPrim < options.tolEps & relResDual < options.tolEps & ~dontstop & it > options.minit) | (resX <options.tolX & it>options.minit))
			break;
		end
	end

	% Final normal map
	N_SH = N;
	N_final = zeros(nrows*ncols,3);
	N_final(imask,:) = N(:,1:3);
	N = reshape(N_final,[nrows ncols 3]);

	% Final depth map
	z_final = zeros(nrows,ncols);
	if(data.K(1,3)>0)
		z_final(imask) = exp(z);
	else
		z_final(imask) = z;
	end
	z = z_final;
	z(data.mask==0) = NaN;

	% Final point cloud
	XYZ = NaN*ones(nrows*ncols,3);
	[xx,yy] = meshgrid(1:ncols,1:nrows);
	xx = xx(imask);
	yy = yy(imask);
	if(data.K(1,3)>0)
		xx = (xx-data.K(1,3))./data.K(1,1);
		yy = (yy-data.K(2,3))./data.K(2,2);
		XYZ(imask,1) = z(imask).*xx;
		XYZ(imask,2) = z(imask).*yy;
		XYZ(imask,3) = z(imask);
	else
		XYZ(imask,1) = xx;
		XYZ(imask,2) = yy;
		XYZ(imask,3) = z(imask);
	end
	XYZ = reshape(XYZ,[nrows ncols 3]);
end




