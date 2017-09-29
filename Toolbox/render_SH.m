function I = render_SH(s,n);

	nchannels = size(s,1);

	[Nx,Ny] = meshgrid(1:n,1:n);
	r = 0.5*n;
	Nx = Nx-r;
	Ny = Ny-r;
	mask = r^2-Nx.^2-Ny.^2 > 0;
	Nz = sqrt(r^2-Nx.^2-Ny.^2);
	Nz(mask==0) = NaN;
	Nz = -Nz;
	norm = sqrt(Nx.^2+Ny.^2+Nz.^2);
	Nx = Nx./norm;
	Ny = Ny./norm;
	Nz = Nz./norm;
	Nx(mask==0) = NaN;
	Ny(mask==0) = NaN;
	
	N = cat(3,Nx,Ny,Nz);

	% Make image 
	I = zeros(n,n,nchannels);
	NxNy = N(:,:,1).*N(:,:,2);
	NxNz = N(:,:,1).*N(:,:,3);
	NyNz = N(:,:,2).*N(:,:,3);
	for ch = 1:nchannels
		%~ I(:,:,ch) = ...
		 %~ max(0,rho(:,:,ch).*(s(ch,1)*N(:,:,1)+s(ch,2)*N(:,:,2)+s(ch,3)*N(:,:,3)...
		 %~ +s(ch,4)...
		 %~ +s(ch,5)*NxNy+s(ch,6)*NxNz+s(ch,7)*NyNz+s(ch,8)*(N(:,:,1).^2-N(:,:,2).^2)+s(ch,9)*(3*N(:,:,3).^2-1)));
		I(:,:,ch) = (s(ch,1)*N(:,:,1)+s(ch,2)*N(:,:,2)+s(ch,3)*N(:,:,3)...
		 +s(ch,4)...
		 +s(ch,5)*NxNy+s(ch,6)*NxNz+s(ch,7)*NyNz+s(ch,8)*(N(:,:,1).^2-N(:,:,2).^2)+s(ch,9)*(3*N(:,:,3).^2-1));
	end
	
end
