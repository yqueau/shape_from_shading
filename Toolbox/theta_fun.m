function [cost,jac] = theta_pq_update_3(theta,zx,zy,u_p,u_q,I,s,rho,beta,xx,yy,K,imask,sfs_pen,area_pen)

	npix = length(imask);
	nchannels = size(I,2);
	
	theta_p = theta(1:npix);
	theta_q = theta(npix+1:2*npix);
	
	dz = sqrt((K(1,1)*theta_p).^2+(K(2,2)*theta_q).^2+(-1-xx.*theta_p-yy.*theta_q).^2);
	d2 = dz.^2;
	d4 = dz.^4;
	N1 = K(1,1)*theta_p./dz;
	N2 = K(2,2)*theta_q./dz;
	N3 = (-1-xx.*theta_p-yy.*theta_q)./dz;
	N4 = ones(npix,1);
	N5 = N1.*N2;
	N6 = N1.*N3;
	N7 = N2.*N3;
	N8 = N1.^2-N2.^2;
	N9 = 3*N3.^2-1;

	dz_p = (K(1,1)*K(1,1)*theta_p-xx.*(-1-xx.*theta_p-yy.*theta_q))./dz;
	dz_q = (K(2,2)*K(2,2)*theta_q-yy.*(-1-xx.*theta_p-yy.*theta_q))./dz;
		
	cost = 0;
	if(nargout>1)
		jac = zeros(2*npix,1);
		
		N1_p = (K(1,1)*dz-K(1,1)*theta_p.*dz_p)./d2;
		N2_p = -(K(2,2)*theta_q.*dz_p)./d2;
		N3_p = (-xx.*dz-(-1-xx.*theta_p-yy.*theta_q).*dz_p)./d2;
		N4_p = zeros(npix,1);
		N5_p = N1_p.*N2+N1.*N2_p;
		N6_p = N1_p.*N3+N1.*N3_p;
		N7_p = N2_p.*N3+N2.*N3_p;
		N8_p = 2*N1.*N1_p-2*N2.*N2_p;
		N9_p = 6*N3.*N3_p;
		
		N1_q = -(K(1,1)*theta_p.*dz_q)./d2;
		N2_q = (K(2,2)*dz-K(2,2)*theta_q.*dz_q)./d2;
		N3_q = (-yy.*dz-(-1-xx.*theta_p-yy.*theta_q).*dz_q)./d2;
		N4_q = zeros(npix,1);
		N5_q = N1_q.*N2+N1.*N2_q;
		N6_q = N1_q.*N3+N1.*N3_q;
		N7_q = N2_q.*N3+N2.*N3_q;
		N8_q = 2*N1.*N1_q-2*N2.*N2_q;
		N9_q = 6*N3.*N3_q;
	end
	for ch = 1:nchannels
		cost_ch = rho(:,ch).*(s(ch,1)*N1+s(ch,2)*N2+s(ch,3)*N3+s(ch,4)*N4+s(ch,5)*N5+s(ch,6)*N6+s(ch,7)*N7+s(ch,8)*N8+s(ch,9)*N9)-I(:,ch);		
		if(nargout>1)
			DFDP = rho(:,ch).*(s(ch,1)*N1_p+s(ch,2)*N2_p+s(ch,3)*N3_p+s(ch,4)*N4_p+s(ch,5)*N5_p+s(ch,6)*N6_p+s(ch,7)*N7_p+s(ch,8)*N8_p+s(ch,9)*N9_p);
			DFDQ = rho(:,ch).*(s(ch,1)*N1_q+s(ch,2)*N2_q+s(ch,3)*N3_q+s(ch,4)*N4_q+s(ch,5)*N5_q+s(ch,6)*N6_q+s(ch,7)*N7_q+s(ch,8)*N8_q+s(ch,9)*N9_q);
			jac = jac+[sfs_pen*DFDP.*cost_ch;sfs_pen*DFDQ.*cost_ch];
		end
		cost = cost+0.5*sfs_pen*sum((cost_ch).^2);
	end
	cost = cost+0.5*beta*sum((zx-theta_p+u_p).^2)+0.5*beta*sum((zy-theta_q+u_q).^2)+area_pen.*sum(dz);
	if(nargout>1)
		jac = jac-[beta*(zx-theta_p+u_p);beta*(zy-theta_q+u_q)];
		jac = jac+area_pen.*[dz_p;dz_q];
	end
end

