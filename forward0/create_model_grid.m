function [model] = create_model_grid(input)
%--------------------------
% model dimension-x
if isempty(input.surf_elec)
    % 只有井中电极
    b_x = unique(input.bore_elec(:,1));
    fprintf('%d boreholes\n',length(b_x));
    x0 = b_x(1):input.unit_a:b_x(end); 
    x1 = x0(1:end-1)+diff(x0)/2; % x1 = [];
    base_x = unique([x0,x1]); 
else
   elec_x = unique(input.elec_xz(:,1))';
   x0 = elec_x(1):input.unit_a:elec_x(end);
   x1 = x0(1:end-1)+diff(x0)/2; % x1 = [];
   base_x = unique([x0,x1,elec_x]);  
end
model.base_x = base_x;
space = 1:3;
left = base_x(1)-fliplr(space)*input.unit_a;
right = base_x(end)+space*input.unit_a;
model.dimension_x = unique([left,base_x,right]); 

% model dimension-z
if isempty(input.bore_elec) % 
    depth_n = zeros(1,100);
    max_d = max(input.unit_a.*input.n)/2*1.05;
    depth_n(1)=0;
    depth_n(2)=input.unit_a*0.25;
    thick = depth_n(2)-depth_n(1);
    tmp_d = max(depth_n);
    id = 3;
    while tmp_d < max_d
        if thick<input.unit_a*5
            thick=depth_n(id-1)-depth_n(id-2);
        end
        depth_n(id)=depth_n(id-1)+thick*1.1;
        id = id+1;
        tmp_d = max(depth_n);
    end
    depth_n = unique(depth_n);
else
   b_z = unique([input.bore_elec(:,2)',0]);
   depth_n = sort(abs(b_z),'ascend'); 
   depth_n = reshape(depth_n,1,[]);
end
thick = max(diff(depth_n));
model.base_z = -depth_n;
space = 1:3;    %[1 3 5 7 9 12 15]; %:7;  %
bottom = depth_n(end)+space*thick;
model.dimension_z = -1*unique([depth_n,bottom]);
model.num_layer = length(depth_n);
% 
[model] = create_model(model);

% mean median mode??
model.mean_res = mean(input.real_data); 
%----------------------------------------------------
model.res_param(1:model.num_param,1) = model.mean_res;
model.res_param_itr = model.res_param;    
%%%%%%%%%%%%%%%%%666%%%%%%%%%%%%%%%%%%%%%%%%
end


function [model] = create_model(model)

model.num_block_x = length(model.dimension_x)-1;
model.num_block_z = length(model.dimension_z)-1;
model.num_param = model.num_block_x*model.num_block_z;
[xxx,zzz] = meshgrid(model.dimension_x,model.dimension_z);
tmp = 1:numel(xxx); ind_m = reshape(tmp,size(xxx));
faces = [];
verts = [xxx(:),zzz(:)];
[r,c] = size(xxx);
tmp = 1;
for cj = 1:c-1
    for ri = 1:r-1
        node1 = ind_m(ri,cj);
        node2 = ind_m(ri,cj+1);
        node3 = ind_m(ri+1,cj+1);
        node4 = ind_m(ri+1,cj);
        faces = [faces;node1,node2,node3,node4];
        params(tmp,:) = [verts(node1,:), verts(node2,:),...
                     verts(node3,:), verts(node4,:)];
        tmp = tmp+1;
    end
end
model.param_xz = params;
model.faces = faces;
model.verts = verts;

xc = mean(params(:,1:2:8),2);
zc = mean(params(:,2:2:8),2);
left = find(xc==min(xc)); right = find(xc==max(xc));
bott = find(zc==min(zc)); 
model.bc = unique([left;right;bott]); %;top]);
%====================================================
% 构建一个网格，为利用插值作图 做准备
[model.grid_x,model.grid_z] = meshgrid(model.base_x,model.base_z);
end




