function model = model_constrains(model,mesh)
%-----------------------------------------
tmp = zeros(mesh.numTri,mesh.numTri);
elemIdx = model.p2c;
neigh2tri = mesh.neigh2tri;
idx = (neigh2tri==0);
neigh2tri(idx)=mesh.numTri;
node = mesh.node';
paramX = node(mesh.tri2node',1);
paramX = reshape(paramX,[],3);
paramZ = node(mesh.tri2node',2);
paramZ = reshape(paramZ,[],3);
xc = mean(paramX,2); zc = mean(paramZ,2);
for i = 1:length(elemIdx)
    ci = elemIdx(i); 
    xci = xc(ci); zci = zc(ci);
    Ai = sqrt(mesh.area(ci));
%     Ai = 1;
    neighs = reshape(neigh2tri(:,ci),1,[]);
    xni = xc(neighs); zni = zc(neighs);
    rni = sqrt((xni-xci).^2 + (zni-zci).^2);
    tmp(ci,neighs) = Ai./rni;
end
c = tmp(elemIdx,elemIdx);
for i = 1:size(c,1)
   c(i,i) =  -sum(c(i,:));
end
model.wc = sparse(c);
ws = sparse(eye(size(model.wc)));
model.ws = ws;
%=============================================
model.Rc = sparse(eye(size(model.wc))); 
model.Rs = sparse(eye(size(model.ws))); 
%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%%%
end





