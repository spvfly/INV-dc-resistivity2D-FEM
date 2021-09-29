function model = model_constrains(model)
%-----------------------------------------
% [wx,wz,ws] = first_order_reg(model);
[wx,wz,ws] = second_order_reg(model);
model.alp_s = 5e-3;
wb = ones(model.num_param,1);
wb(model.bc)=1;
wb = sparse(diag(wb));
model.wx = wx*wb;
model.wz = wz*wb;
model.ws = ws*wb;
%-----------------------------------------
model.Rx = sparse(eye(size(model.wx))); 
model.Rz = sparse(eye(size(model.wz)));
model.Rs = sparse(eye(size(model.wz)));
%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [cx,cz,cs] = first_order_reg(model)
% DCIP2D
ind = 1:model.num_param;
indm = reshape(ind,model.num_block_z,model.num_block_x);
dx = abs(model.param_xz(:,1)-model.param_xz(:,5));
dz = abs(model.param_xz(:,2)-model.param_xz(:,6));

%  --x1--x2----
cx = zeros(model.num_param,model.num_param);
indx1 = indm(:,1:end-1); indx1 = indx1(:);
indx2 = indx1+model.num_block_z;
dx1 = dx(indx1); dx2 = dx(indx2);

dz1 = dz(indx1);
fx = sqrt(2*dz1./(dx1+dx2));
ind = sub2ind(size(cx),indx1,indx1);
cx(ind) = -fx; 
ind = sub2ind(size(cx),indx1,indx2);
cx(ind) = fx;

cx = sparse(cx);
%----------z1-z2----------------------------
% cz
cz = zeros(model.num_param,model.num_param);
indz1 = indm(1:end-1,:); indz1 = indz1(:);
indz2 = indz1+1;  
dz1 = dz(indz1); dz2 = dz(indz2);

dx1 = dx(indz1);
fz = sqrt(2*dx1./(dz1+dz2));
ind = sub2ind(size(cz),indz1,indz1);
cz(ind) = -fz;
ind = sub2ind(size(cz),indz1,indz2);
cz(ind) = fz;

cz = sparse(cz);

cs = reshape(sqrt(dx.*dz),[],1);
cs = sparse(diag(cs));
end

function [cx,cz,cs] = second_order_reg(model)
%
ind = 1:model.num_param;
indm = reshape(ind,model.num_block_z,model.num_block_x);
dx = abs(model.param_xz(:,1)-model.param_xz(:,5));
dz = abs(model.param_xz(:,2)-model.param_xz(:,6));
A = sqrt(dx.*dz);
%  --x1--x2--x3--
cx = zeros(model.num_param,model.num_param);

indx1 = indm(:,1:end-2); indx1 = indx1(:);
indx2 = indm(:,2:end-1); indx2 = indx2(:);
indx3 = indm(:,3:end);   indx3 = indx3(:);
%  --x1--x2--x3--
dx1 = dx(indx1); dx2 = dx(indx2); dx3 = dx(indx3);
Ax = A(indx2);
a = 4.*Ax./((dx1+2.*dx2+dx3).*(dx1+dx2)); 
ind = sub2ind(size(cx),indx2,indx1);
cx(ind) = a;
e = 4.*Ax./((dx1+2.*dx2+dx3).*(dx2+dx3));
ind = sub2ind(size(cx),indx2,indx3);
cx(ind) = e;

ind = sub2ind(size(cx),indx2,indx2);
cx(ind) = -(a+e);

% % ×ó±ß
indx2 = indm(:,1); 
indx3 = indm(:,2);
dx2 = dx(indx2); dx3 = dx(indx3);
Ax = A(indx2);
% 
e = 2.*Ax./(dx2+dx3);
ind = sub2ind(size(cx),indx2,indx3);
cx(ind) = e;
ind = sub2ind(size(cx),indx2,indx2);
cx(ind) = -e;
% ÓÒ±ß
indx1 = indm(:,end-1);
indx2 = indm(:,end); 
Ax = A(indx2);
dx1 = dx(indx1); dx2 = dx(indx2);
% 
a = 2.*Ax./(dx1+dx2);
ind = sub2ind(size(cx),indx2,indx1);
cx(ind) = a;
ind = sub2ind(size(cx),indx2,indx2);
cx(ind) = -a;

cx = sparse(cx);
%----------z1-z2----------------------------
% cz
cz = zeros(model.num_param,model.num_param);

indz1 = indm(1:end-2,:); indz1 = indz1(:);
indz2 = indm(2:end-1,:); indz2 = indz2(:);
indz3 = indm(3:end,:);   indz3 = indz3(:);
dz1 = dz(indz1); dz2 = dz(indz2); dz3 = dz(indz3);
Az = A(indz2);
b = 4.*Az./((dz1+2.*dz2+dz3).*(dz1+dz2));
ind = sub2ind(size(cz),indz2,indz1);
cz(ind) = b;
d = 4.*Az./((dz1+2.*dz2+dz3).*(dz2+dz3));
ind = sub2ind(size(cz),indz2,indz3);
cz(ind) = d;

ind = sub2ind(size(cz),indz2,indz2);
cz(ind) = -(b+d);
% 
% ÉÏ±ß
indz2 = indm(1,:); 
indz3 = indz2 + 1;
Az = A(indz2);
dz2 = dz(indz2); dz3 = dz(indz3);
d = 2.*Az./(dz2+dz3);
ind = sub2ind(size(cz),indz2,indz3);
cz(ind) = d;
ind = sub2ind(size(cz),indz2,indz2);
cz(ind) = -d;

% ÏÂ±ß
indz2 = indm(end,:);  
indz1 = indz2-1;
Az = A(indz2); 
dz1 = dz(indz1); dz2 = dz(indz2);
b = 2.*Az./(dz1+dz2);
ind = sub2ind(size(cz),indz2,indz1);
cz(ind) = b;
ind = sub2ind(size(cz),indz2,indz2);
cz(ind) = -b;
cz = sparse(cz);

cs = reshape(A,[],1);
cs = sparse(diag(cs));
end





