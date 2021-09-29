function PlotMesh(mesh)
TRIA = mesh.tri2node';
VERT = mesh.node';
% EDGE = mesh.edge';
patch('faces',TRIA(:,1:3),'vertices',VERT, ...
    'facecolor','w', ...
    'edgecolor',[.2,.2,.2]) ;
% hold on; 
probe = mesh.node(:,mesh.mes_node);
hold on
plot(probe(1,:),probe(2,:),'ro');
minx = min(VERT(:,1)); maxx = max(VERT(:,1));
minz = min(VERT(:,2)); maxz = max(VERT(:,2));
xlim([minx-5 maxx+5]);
zlim([minz-5 maxz+5]);
axis equal;
colorbar off;
end