function PlotModel(model,input)
faces = model.faces;
verts = model.verts;
empty = (model.res_param==0);
model.res_param(empty) = model.mean_res;
patch('faces',faces,'vertices',verts,'facecolor','flat','edgecolor','k',...
    'FaceVertexCData',model.res_param);

x1 = min(model.base_x); x2 = max(model.base_x);
z1 = max(model.base_z); z2 = min(model.base_z);
pp = [x1,z1; x1,z2; x2,z2; x2,z1];
hold on
plot(pp(:,1),pp(:,2),'r','LineWidth',2);
plot(input.elec_xz(:,1),input.elec_xz(:,2),'o')
colorbar
axis equal;
ext = 1;
xlim([min(model.dimension_x)-ext, max(model.dimension_x)+ext]);
% 注意，画图时深度是负值
ylim([min(model.dimension_z)-ext, max(model.dimension_z)+ext]);
end