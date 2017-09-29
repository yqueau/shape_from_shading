function s = estimate_lighting(I,rho,z,mask,K,nb_harmonics)

	disp('Estimating lighting');

	% Some variables
	mask(z<=0) = 0;
	imask = find(mask>0);
	npix = length(imask);
	nrows = size(I,1);
	ncols = size(I,2);
	nchannels = size(I,3);
	if(size(K,1)>0)
		[xx,yy] = meshgrid(1:ncols,1:nrows);
		xx = xx(imask);
		xx = xx-K(1,3);
		yy = yy(imask);
		yy = yy-K(2,3);
		z = log(z);
	else
		xx = zeros(npix,1);
		yy = zeros(npix,1);
		K(1,1) = 1;
		K(2,2) = 1;
	end


	% Compute depth gradient
	G = make_gradient(mask); % Finite differences stencils
	Dx = G(1:2:end-1,:);
	Dy = G(2:2:end,:);
	zx = Dx*z(imask);
	zy = Dy*z(imask);

	% Normals
	N = zeros(npix,9); 
	dz = max(eps,sqrt((K(1,1)*zx).^2+(K(2,2)*zy).^2+(-1-xx.*zx-yy.*zy).^2));
	N(:,1) = K(1,1)*zx./dz;
	N(:,2) = K(2,2)*zy./dz;
	N(:,3) = (-1-xx.*zx-yy.*zy)./dz;
	N(:,4) = 1;
	N(:,5) = N(:,1).*N(:,2);
	N(:,6) = N(:,1).*N(:,3);
	N(:,7) = N(:,2).*N(:,3);
	N(:,8) = N(:,1).^2-N(:,2).^2;
	N(:,9) = 3*N(:,3).^2-1;

	% Estimate lighting
	for ch = 1:nchannels
		Ich = I(:,:,ch);
		rhoch = rho(:,:,ch);
		Ich = Ich(imask);
		rhoch = rhoch(imask);
		s(ch,1:nb_harmonics) = transpose(bsxfun(@times,rhoch,N(:,1:nb_harmonics)) \ Ich); 
	end
	s(:,nb_harmonics+1:9) = 0;
end
