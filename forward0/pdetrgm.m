function fem = pdetrgm(mesh)
% 单元分析

% Corner point indices  角点索引
% 网格编码为逆时针
%                       I<-1
%                         *       
%                       /   \    
%                      /     \
%                     *-------*
%                   J->2       3->M  
I=mesh.tri2node(1,:);
J=mesh.tri2node(2,:);
M=mesh.tri2node(3,:);
% L = a + bx + cy  广义坐标直接构造形函数
bi=mesh.node(2,J)-mesh.node(2,M); ci=mesh.node(1,M)-mesh.node(1,J);
bj=mesh.node(2,M)-mesh.node(2,I); cj=mesh.node(1,I)-mesh.node(1,M);
bm=mesh.node(2,I)-mesh.node(2,J); cm=mesh.node(1,J)-mesh.node(1,I);
% Calculate area of element Delta^e
area_element =(bj.*cm-bm.*cj)./2;
fem.area = area_element; 
fem.b1 = bi; fem.c1 = ci;
fem.b2 = bj; fem.c2 = cj;
fem.b3 = bm; fem.c3 = cm;
% 
fem.kIndx = [I',I';I',J';I',M';...
             J',I';J',J';J',M';...
             M',I';M',J';M',M'];
%-------------------------------------------
% Nb 边界上的有限元
edge = mesh.edge2node(:,mesh.bdEdgeId);
fem.BCi = repmat(reshape(edge, 2, 1, []), [1, 2, 1]);
fem.BCj = repmat(reshape(edge, 1, 2, []), [2, 1, 1]);
L_ref = [2,1;1,2];
numEdge = length(edge);
L = zeros(2,2,numEdge);
for i = 1:numEdge
    n1 = edge(1,i); n2 = edge(2,i);
    x1 = mesh.node(1,n1); z1 = mesh.node(2,n1);
    x2 = mesh.node(1,n2); z2 = mesh.node(2,n2);
    dL = sqrt( (x1-x2)^2 + (z1-z2)^2 );
    L(:,:,i) = dL/6*L_ref;
end
fem.Le = L;
end

