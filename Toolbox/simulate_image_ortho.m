function I_simulated = simulate_image_portho(u,mask,omega,delta)

	I_simulated = zeros*ones(size(mask,1),size(mask,2),4);
	
	u(mask==0) = NaN;
	
	% Forward x, forward y
	ux = (u(:,[2:end end])-u)./delta;
	uy = (u([2:end end],:)-u)./delta;
	N = cat(3,-uy,-ux,ones(size(mask)));
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,1) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	
	% Forward x, backward y
	ux = (u(:,[2:end end])-u)./delta;
	uy = (u-u([1 1:end-1],:))./delta;
	N = cat(3,-uy,-ux,ones(size(mask)));
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,2) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	% Backward x, backward y
	ux = (u-u(:,[1 1:end-1]))./delta;
	uy = (u-u([1 1:end-1],:))./delta;
	N = cat(3,-uy,-ux,ones(size(mask)));
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,3) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	% Backward x, forward y
	ux = (u-u(:,[1 1:end-1]))./delta;
	uy = (u([2:end end],:)-u)./delta;
	N = cat(3,-uy,-ux,ones(size(mask)));
	N = bsxfun(@rdivide,N,sqrt(sum(N.^2,3)));
	I_simulated(:,:,4) = omega(1).*N(:,:,1)+omega(2).*N(:,:,2)+omega(3).*N(:,:,3);
	
	I_simulated(isnan(I_simulated)) = Inf;
	I_simulated = min(I_simulated,[],3);
	I_simulated(I_simulated==Inf) = 0;
end
