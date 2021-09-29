function model = link_model_mesh(input,mesh) 
mean_rho = mean(input.real_data);
model.mean_res = mean_rho;
model.p2c = find(mesh.subdm==1);
num_param = length(model.p2c);
model.num_param = num_param;
model.res_param(1:num_param,1) = model.mean_res;
model.res_param_itr = model.res_param;
model.mref = [];
% =========================================
% about plot resitivity contour 
cellIdx = mesh.tri2node(:,model.p2c)';
node_x = mesh.node(1,:)';
paramX = node_x(cellIdx);
node_z = mesh.node(2,:)';
paramZ = node_z(cellIdx);

xc = mean(paramX,2);
zc = mean(paramZ,2);
model.param_xz = [xc,zc];
% gridx,gridz
min_x = min(input.elec_xz(:,1)); max_x = max(input.elec_xz(:,1));
L = abs(max_x-min_x);
xStep = L/200;
baseX = min_x:xStep:max_x;
%
if isempty(input.bore_elec)
    dep = L*0.17;
    zStep = dep/100;
    baseZ = sort((-dep:zStep:0),'descend'); 
else
    dd = min(input.bore_elec(:,2))/100;
    baseZ = 0:dd:min(input.bore_elec(:,2));
end

[grid_x,grid_z] = meshgrid(baseX,baseZ);
model.grid_x = grid_x; 
model.grid_z = grid_z; 

end

