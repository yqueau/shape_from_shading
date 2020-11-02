function I_simulated = simulate_image(x,y,u,f,mask,omega,delta)

	I_simulated = zeros*ones(size(mask,1),size(mask,2),4);
	
	u(mask==0) = NaN;
	
	% Forward x, forward y
	uy = (u(:,[2:end end])-u)./delta;
	ux = (u([2:end end],:)-u)./delta;
	N = cat(3,-ux,-uy,u/f+x.*ux/f+y.*uy/f);
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,1) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	
	% Forward x, backward y
	uy = (u(:,[2:end end])-u)./delta;
	ux = (u-u([1 1:end-1],:))./delta;
	N = cat(3,-ux,-uy,u/f+x.*ux/f+y.*uy/f);
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,2) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	% Backward x, backward y
	uy = (u-u(:,[1 1:end-1]))./delta;
	ux = (u-u([1 1:end-1],:))./delta;
	N = cat(3,-ux,-uy,u/f+x.*ux/f+y.*uy/f);
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,3) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	% Backward x, forward y
	uy = (u-u(:,[1 1:end-1]))./delta;
	ux = (u([2:end end],:)-u)./delta;
	N = cat(3,-ux,-uy,u/f+x.*ux/f+y.*uy/f);
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,4) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	I_simulated(isnan(I_simulated)) = Inf;
	I_simulated = min(I_simulated,[],3);
	I_simulated(I_simulated==Inf) = 0;
end
