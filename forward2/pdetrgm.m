function fem = pdetrgm(mesh)
%
Bk = zeros(2, 2, mesh.numTri);
Binv = zeros(2, 2, mesh.numTri);
detBk = zeros(mesh.numTri, 1);
area = zeros(mesh.numTri, 1);

get_Bk = @(i) [ ...
    mesh.node(1, mesh.tri2node(1, i)) - mesh.node(1, mesh.tri2node(3, i)) mesh.node(2, mesh.tri2node(1, i)) - mesh.node(2, mesh.tri2node(3, i)); ...
    mesh.node(1, mesh.tri2node(2, i)) - mesh.node(1, mesh.tri2node(3, i)) mesh.node(2, mesh.tri2node(2, i)) - mesh.node(2, mesh.tri2node(3, i)) ...
    ]; 
for ii = 1:mesh.numTri
    Bk(:, :, ii) = get_Bk(ii);
    Binv(:, :, ii) = inv(Bk(:, :, ii)); 
    detBk(ii) = det(Bk(:, :, ii));  
    if detBk(ii) < 0 
        mesh.tri2node([2, 3], ii) = mesh.tri2node([3, 2], ii);
        fprintf('change tri2node\n');
    end
    area(ii) = 0.5 * abs(detBk(ii));
end
fem.Bk = Bk;
fem.Binv = Binv;
fem.detBk = detBk;
fem.area = reshape(area,1,[]);

% fem.dimension = 2;   fem.order = 1;
% Size of element matrices is nD x nD: nD为每个三角单元的节点数
% nD = (fem.dimension + 1) + (fem.order - 1) * fem.dimension * (fem.dimension + 1) / 2;
nD = 3;
% N(L1,L2) ==> NL = (L1,L2,1-L1-L2)
% N'*N==>积分解析公式
M_ref = 1 / 24 * ...
        [2 1 1; ...
        1 2 1; ...
        1 1 2];
% gradN(L1,L2)
gradNL = [1,0,-1; 0,1,-1];

nEl = mesh.numTri;% Number of triangular elements

% Allocate storage for all element matrices in nD x nD x nEl arrays
K = zeros(nD, nD, nEl);  M = zeros(nD, nD, nEl);
% !!!!!!!!!!!!!!!!
fem.Ci = repmat(reshape(mesh.tri2node, nD, 1, []), [1, nD, 1]); % 节点的横向索引
fem.Cj = repmat(reshape(mesh.tri2node, 1, nD, []), [nD, 1, 1]); % 节点的纵向索引
for i = 1:nEl
    K_k = zeros(nD, nD);
    d = 2*area(i);
    Binv = fem.Binv(:, :, i);

    B = Binv * gradNL;  % grad(N_1,2,3)
    % Gauss积分点 L1=1/3,L2=1/3,w=1/2;
    K_k = K_k + 0.5 * (B' * B); 
 
    K(:, :, i) = K_k   * d;
    M(:, :, i) = M_ref * d; 
end
fem.Ke1s = K;
fem.Ke2s = M;

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

% n1 = edge(1,i); n2 = edge(2,i);
% x1 = mesh.node(1,n1); z1 = mesh.node(2,n1);
% x2 = mesh.node(1,n2); z2 = mesh.node(2,n2);
% dL = sqrt( (x1-x2)^2 + (z1-z2)^2 );
% src = mean(mesh.node(:,mesh.mes_node),2);
% eec = [0.5*(x1+x2);0.5*(z1+z2)];
% dx = abs(src(1)-eec(1));
% dz = abs(src(2)-eec(2));
% r = sqrt( dx^2 + dz^2 );
% cos = max(dx/r,dz/r);
% kr = fem.k(kit)*r;
% if kr>600,kr=600;end
% bk1 = besselk(1,kr)*cos; %K1
% bk0 = besselk(0,kr); %K0



% if strcmp(outputformat, 'matrices')
%     if verbose
%         fprintf('  Assembly of system matrices...');
%         t0 = tic();
%     end
%     
%     % Multiply with parameter
%     fem.Ke = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, K, reshape(1 ./ fem.mu, 1, 1, []))), nDOFS, nDOFS);
%     fem.Me = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, M, reshape(fem.sigma, 1, 1, []))), nDOFS, nDOFS);
%     
%     fem.Kh = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, K, reshape(1 ./ fem.sigma, 1, 1, []))), nDOFS, nDOFS);
%     fem.Mh = sparse(Ci(:), Cj(:), tools.asColumn(bsxfun(@times, M, reshape(fem.mu, 1, 1, []))), nDOFS, nDOFS);
% end