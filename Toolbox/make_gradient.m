% Functions for computing the gradient operator on non-rectangular domains
function [M,imask] = make_gradient(mask)

	% Compute forward (Dxp and Dyp) and backward (Dxm and Dym) operators
	[Dyp,Dym,Dxp,Dxm,Sup,Sum,Svp,Svm,Omega,index_matrix,imask] = gradient_operators(mask);
	[nrows,ncols] = size(mask);

	% When there is no bottom neighbor, replace by backward (or by 0 if no top)
	Dy = Dyp;
	no_bottom = find(~Omega(:,:,1));
	no_bottom = nonzeros(index_matrix(no_bottom));
	Dy(no_bottom,:) = Dym(no_bottom,:);

	% Same for the x direction (right / left)
	Dx = Dxp;
	no_right = find(~Omega(:,:,3));
	no_right = nonzeros(index_matrix(no_right));
	Dx(no_right,:) = Dxm(no_right,:);

	M = sparse([],[],[],2*size(Dx,1),size(Dx,2),2*length(imask));
	M(1:2:end-1,:) = Dx;
	M(2:2:end,:) = Dy;
end

function [Dup,Dum,Dvp,Dvm,Sup,Sum,Svp,Svm,Omega,index_matrix,imask] = gradient_operators(mask)

	[nrows,ncols] = size(mask);
	Omega_padded = padarray(mask,[1 1],0);

	% Pixels who have bottom neighbor in mask
	Omega(:,:,1) = mask.*Omega_padded(3:end,2:end-1);
	% Pixels who have top neighbor in mask
	Omega(:,:,2) = mask.*Omega_padded(1:end-2,2:end-1);
	% Pixels who have right neighbor in mask
	Omega(:,:,3) = mask.*Omega_padded(2:end-1,3:end);
	% Pixels who have left neighbor in mask
	Omega(:,:,4) = mask.*Omega_padded(2:end-1,1:end-2);
	

	imask = find(mask>0);
	index_matrix = zeros(nrows,ncols);
	index_matrix(imask) = 1:length(imask);

	% Dv matrix
	% When there is a neighbor on the right : forward differences
	idx_c = find(Omega(:,:,3)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc,yc+1));
	indices_right = indices_right(:);
	II = indices_centre;
	JJ = indices_right;
	KK = ones(length(indices_centre),1);
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;-ones(length(indices_centre),1)];
	
	Dvp = sparse(II,JJ,KK,length(imask),length(imask));
	Svp = speye(length(imask));
	Svp = Svp(index_matrix(idx_c),:);

	% When there is a neighbor on the left : backward differences
	idx_c = find(Omega(:,:,4)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc,yc-1));
	indices_right = indices_right(:);
	II = [indices_centre];
	JJ = [indices_right];
	KK = [-ones(length(indices_centre),1)];
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;ones(length(indices_centre),1)];
	
	Dvm = sparse(II,JJ,KK,length(imask),length(imask));
	Svm = speye(length(imask));
	Svm = Svm(index_matrix(idx_c),:);



	% Du matrix
	% When there is a neighbor on the bottom : forward differences
	idx_c = find(Omega(:,:,1)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc+1,yc));
	indices_right = indices_right(:);
	II = indices_centre;
	JJ = indices_right;
	KK = ones(length(indices_centre),1);
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;-ones(length(indices_centre),1)];
	
	Dup = sparse(II,JJ,KK,length(imask),length(imask));
	Sup = speye(length(imask));
	Sup = Sup(index_matrix(idx_c),:);

	% When there is a neighbor on the top : backward differences
	idx_c = find(Omega(:,:,2)>0);
	[xc,yc] = ind2sub(size(mask),idx_c);
	indices_centre = index_matrix(idx_c);	
	indices_right = index_matrix(sub2ind(size(mask),xc-1,yc));
	indices_right = indices_right(:);
	II = [indices_centre];
	JJ = [indices_right];
	KK = [-ones(length(indices_centre),1)];
	II = [II;indices_centre];
	JJ = [JJ;indices_centre];
	KK = [KK;ones(length(indices_centre),1)];
	
	Dum = sparse(II,JJ,KK,length(imask),length(imask));
	Sum = speye(length(imask));
	Sum = Sum(index_matrix(idx_c),:);	

end
