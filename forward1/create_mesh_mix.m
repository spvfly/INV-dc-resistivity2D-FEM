function [model,mesh] = create_mesh_mix(input,model)
% base-x
dd = min(diff(model.base_x))/2;
base_x = min(model.dimension_x):dd:max(model.dimension_x);
base_x = unique([model.dimension_x,base_x]);
% base-z
if isempty(input.bore_elec)
    base_z = (model.dimension_z);
    % 外扩
    L = max(input.elec_xz(:,1))-min(input.elec_xz(:,1));
else
    dd = min(diff(model.base_z))/2;
    base_z0 = max(model.base_z):dd:min(model.base_z);
    base_z = unique([model.dimension_z,base_z0]);
    base_z = sort(base_z,'descend');
    Lx = max(input.elec_xz(:,1))-min(input.elec_xz(:,1));
    Lz = max(input.bore_elec(:,2))-min(input.bore_elec(:,2));
    L  = max(Lx,Lz);
end
ext_x = 5*L;  ext_z = 5*L;
% 大网格的范围
limit = [min(base_x)-ext_x,max(base_x)+ext_x,max(base_z),min(base_z)-ext_z];
%--------------------------------------------------
% 网格扩展扩展  用非结构网格镶嵌
%    ---->
% 1------------2                         m+ccc------------m+n
% |            |                          |                |
% |            3                         m+cc              |
% |            |                          |                |
% ...         ...                       ...              ...
% |            |                          |                |
% |            m-m+1--...----------------m+c               |
% |                                                        |
% m+n+2--------------------------------------------------m+n+1
%      <-----------<------------------------

node(1,1) = limit(1);
node(1,2) = limit(3);
% 小网格左边
zi = base_z;
x = min(base_x)*ones(size(zi));
node = [node;x(:),zi(:)];

xi = base_x(2:end-1);
z = min(base_z)*ones(size(xi));
node = [node;xi(:),z(:)];

zi = fliplr(base_z);
x = max(base_x)*ones(size(zi));
node= [node;x(:),zi(:)];

tmp = length(node)+1;
node(tmp,1) = limit(2);
node(tmp,2) = limit(3);
tmp = tmp + 1;
node(tmp, 1) = limit(2);
node(tmp, 2) = limit(4);
tmp = tmp + 1;
node(tmp,1) = limit(1);
node(tmp,2) = limit(4);
%-----------------------------------------------------------
edge = [];
for i = 1:length(node)-1
    edge = [edge; i,i+1];
end
edge = [edge;length(node),1];
faces{1} = [1:length(edge)];
warning off
addpath ('mesh2d');
options.output=false;
[VERT ,TRIA,fnum]=meshfaces(node,edge,faces,[],options);
warning on
bnode = node(2:end-3,:);
bind = [];
for i = 1:length(bnode)
   ind = find(VERT(:,1)==bnode(i,1)&VERT(:,2)==bnode(i,2));
   bind = [bind;ind];
end
%-----------------------------------------------------------
%  网格四边形的 （反演模型区域内）
%             1---4
%             |\/ |
%             | 5 | 
%             |/\ |
%             2---3
loc_base_x = base_x(2:end-1); loc_base_z = base_z(1:end-1);
[tmpx,tmpz] = meshgrid(loc_base_x,loc_base_z);
node = [VERT;tmpx(:),tmpz(:)];
%-----按照节点编号对齐
ind = length(VERT)+(1:numel(tmpx));
indm = reshape(ind,size(tmpx));

lb = 1:length(loc_base_z);
indm = [bind(lb),indm];
rb = (length(bind)-length(loc_base_z)+1):length(bind);
indm = [indm,flipud(bind(rb))];
bb = lb(end)+(1:length(base_x));
indm = [indm;bind(bb)'];
[r,c] = size(indm);
tmpind = length(node)+1;
faces = TRIA;
for ri = 1:r-1
    for cj = 1:c-1
        n1i = indm(ri,cj);     n4i = indm(ri,cj+1);
        n2i = indm(ri+1,cj);   n3i = indm(ri+1,cj+1);
        n5i = tmpind; 
        tmp = mean( node([n1i,n2i,n3i,n4i],:),1);
        node = [node;tmp];
        ele1 = [n1i,n2i,n5i];
        ele2 = [n2i,n3i,n5i];
        ele3 = [n3i,n4i,n5i];
        ele4 = [n4i,n1i,n5i];
        faces = [faces; ele1; ele2; ele3; ele4];
        tmpind = tmpind+1;
    end
end
mesh.node = node';
mesh.numNode = length(mesh.node);
mesh.tri2node = faces';
mesh.numTri = length(faces);
%------------------------------------------------------------------%
%        link model and mesh
%------------------------------------------------------------------%
elec_x = input.elec_xz(:,1); 
elec_z = input.elec_xz(:,2);
for j=1:input.num_elec
    ind=find(elec_x(j)==mesh.node(1,:) & elec_z(j)==mesh.node(2,:));
    if ~isempty(ind)
        mesh.mes_node(j)=ind;
    end
end

if length(mesh.mes_node)~= input.num_elec
    tt = sprintf(['function [model,mesh]=create_mesh(model)==>\n',...
            '\t function [model,mesh] = create_mesh_normal(model)==>\n',...
            'Something went wrong']);
    errordlg(tt);
    pause
end

e1=mesh.tri2node(1,:); e2=mesh.tri2node(2,:); e3=mesh.tri2node(3,:); 
xe=(mesh.node(1,e1) + mesh.node(1,e2) + mesh.node(1,e3))/3;
ze=(mesh.node(2,e1) + mesh.node(2,e2) + mesh.node(2,e3))/3;
cell2param_ind = zeros(mesh.numTri,1);
for j=1:model.num_param
     minx = min(model.param_xz(j,1:2:8)); minz = min(model.param_xz(j,2:2:8));
     maxx = max(model.param_xz(j,1:2:8)); maxz = max(model.param_xz(j,2:2:8));
     ind=((xe>=minx) & (xe<=maxx) & (ze<=maxz) & (ze>=minz));
     cell2param_ind(ind)=j; %
end

mesh.c2p  = cell2param_ind;
model.p2c = find(cell2param_ind~=0);
mesh.rho(1:mesh.numTri) = model.mean_res;
mesh = calc_edge(mesh,model);
for i = 1:model.num_param
    ind = (mesh.c2p==i);
    mesh.rho(ind) = model.res_param(i);
end
%%%%%%%%%%%%%%function end!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function mesh = calc_edge(mesh,model)
% The following array associates          node1 *
% the three edges of a triangle     edge3     /   \  edge2
% with the three local vertex                /     \
% indices (see description above).    node2 *-------* node3
%                                            edge1
edgesByVert = [ 2 3; 3 1; 1 2]; 
edges = mesh.tri2node(edgesByVert', :); % (2 x 3)*nTri 每个单元三条边
edges = reshape(edges, 2, 3 * size(mesh.tri2node, 2)); % 2 x （3*nTri）
edges = sort(edges, 1);  % 按行排序
[edges, ~, idx] = unique(edges', 'rows'); %去重，得到边的集合，和相对应的序号
mesh.edge2node = edges'; 
% tri2edge=[edge1;edge2;edge3]
mesh.tri2edge = reshape(idx, 3, size(mesh.tri2node, 2)); 
%========find boundary edge=============
edgeId = mesh.tri2edge(:);
counts = accumarray(edgeId,1);
bdEdgeId = find(counts==1)';
% 除去top自由边界,利用节点标记挑选边界
edges = mesh.edge2node(:,bdEdgeId);
node1 = edges(1,:);  node2 = edges(2,:);
x1 = mesh.node(1,node1); x2 = mesh.node(1,node2);
xc = 0.5*(x1+x2);
z1 = mesh.node(2,node1); z2 = mesh.node(2,node2);
zc = 0.5*(z1+z2);
I1 = find(zc~=0); % 去除top
I2 = find(xc~=min(model.dimension_x) & xc~=max(model.dimension_x));
I = intersect(I1,I2);
mesh.bdEdgeId = bdEdgeId(I);
tmp = mesh.edge2node(:,mesh.bdEdgeId);
mesh.bdNodeId = unique(tmp(:));
end
