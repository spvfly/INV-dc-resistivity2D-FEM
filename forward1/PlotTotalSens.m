function PlotTotalSens(model,fem)
R = sqrt(diag(fem.array_jacobian'*fem.array_jacobian));
xc = mean(model.param_xz(:,1:2:8),2); %xc = model.paramXZ(:,1);
zc = mean(model.param_xz(:,2:2:8),2); %zc = model.paramXZ(:,2);
ZI=TriScatteredInterp(xc,zc,R);
ZII = ZI(model.grid_x,model.grid_z);
contourf(model.grid_x,model.grid_z,ZII,17,'EdgeColor','none');
colormap('jet');
end
