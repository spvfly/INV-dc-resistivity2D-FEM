function PlotTotalSens(model,fem)
R = sqrt(diag(fem.array_jacobian'*fem.array_jacobian));
xc = model.paramXZc(:,1);
zc = model.paramXZc(:,2);
ZI=TriScatteredInterp(xc,zc,R,'natural');
ZII = ZI(model.gridX,model.gridZ);
contourf(model.gridX,model.gridZ,ZII)
% contourf(model.gridX,model.gridZ,ZII,17,'EdgeColor','none');
colormap('jet');
end

% function PlotTotalSens(model,mesh,fem)
% R = sqrt(diag(fem.array_jacobian'*fem.array_jacobian));
% faces = mesh.tri2node(:,model.p2c)';
% verts = mesh.node';
% patch('faces',faces,'vertices',verts,'facecolor','flat','edgecolor','k',...
%     'FaceVertexCData',R);
% % contourf(model.gridX,model.gridZ,ZII,17,'EdgeColor','none');
% colormap('jet');
% end