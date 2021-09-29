function fem=compute_femAndJ(input,model,mesh,fem,jc_flag,itr)
for i = 1:model.num_param
    ind = (mesh.c2p==i);
    mesh.rho(ind) = model.res_param_itr(i);
end
% matrix preallocation 矩阵预分配
array_jac=zeros(input.num_mes,model.num_param);
jac_tmp = zeros(input.num_mes,length(model.p2c));
num_elecs = length(mesh.mes_node);
aaa=zeros(num_elecs,mesh.numNode);
%----------有限元正演+完全的jacob矩阵计算----------------
for kit=1:length(fem.k)
    %-----------画进度条
    h2 = waitbar(kit/length(fem.k));
    drawnow;
    %------组成刚度矩阵并求有限元解
    [tmp_aaa]=fem_calc_kit(kit,mesh,fem);
    fem.tmp_aaa{kit}=tmp_aaa;
    %-----Inverse fourier transform 时间域的有限元解
    aaa = aaa+fem.g(kit).*tmp_aaa/2;%/pi; 
    %---求Jacobian矩阵(这是按照正演公式求的)
    if jc_flag==2 || itr==0 
       % 注意，还在波数域下！！！
       trans_jac1=jacob_calc_kit(kit,input,model,mesh,tmp_aaa,fem);
       % Inverse fourier transform
       jac_tmp=jac_tmp+fem.g(kit)*trans_jac1/2; %/pi;
    end
end
fem.aaa = aaa;  % 数值模拟解（电位矩阵）
close(h2);
%-----------------------------------------------
%-------计算视电阻率结果
app_res_itr = array_app_res(input,mesh,aaa);

% 计算均匀半空间的解析解
if itr==0
    fem = homog_field(input,model,mesh,fem);
    app_res_ana = array_app_res(input,mesh,fem.Ap);
    % 视电阻率的修正的系数
    fem.normalized = app_res_ana./app_res_itr; %%
end
fem.app_res_itr = app_res_itr;
% -------jacob矩阵，jacob==>log(jacob)----------------
if (jc_flag==2 || itr==0)
    c2p = nonzeros(mesh.c2p);
    for i = 1:input.num_mes
        jacb_mes = accumarray(c2p,jac_tmp(i,:));
        jacb_mes = reshape(jacb_mes,1,[]);
        app_res = fem.app_res_itr(i);    % 视电阻率
        sig = 1./model.res_param';  % 模型电阻率（真）
        array_jac(i,:) = jacb_mes.*(sig/app_res);
    end
    fem.array_jacobian = array_jac;
end

if(jc_flag==1 && itr>0)
    fem=quasi_newton(model,fem);
end

%%%%%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function fem=quasi_newton(model,fem)
    % 李天成. 电阻率成像技术的二维三维正反演研究[D]. 
    % 中国地质大学（北京）, 2008.
    dd = model.dx1'*model.dx1;
    dm=log10(fem.app_res_itr) - log10(fem.app_res) - fem.array_jacobian*model.dx1;
    fem.array_jacobian = fem.array_jacobian + dm*model.dx1'/dd;
end
%%
% *********************************************************************
% *   fem_calc( )构建整体刚度矩阵 并 求解线性方程组                    *
% *-------------------------------------------------------------------*
% * This function forms the global stiffness matrix.                  *
% *                                                                   *
% * This function calculate the stiffnes term for every element.      *
% * It is valid only for linear triangular elements                   *
% *********************************************************************/
function [tmp_aaa]=fem_calc_kit(kit,mesh,fem)
% Start to assemble all area elements in Omega 
sig = 1./mesh.rho;
%有限元第一项系数
p1 = reshape(sig,1,1,[]);
% 有限元第二项系数
p2 = sig.*fem.k(kit).^2;
p2 = reshape(p2,1,1,[]);
nDOFS = mesh.numNode;
asColumn = @(x) reshape(x,[],1);
K1 = sparse(fem.Ci(:), fem.Cj(:), ...
    asColumn(bsxfun(@times, fem.Ke1s, p1)),nDOFS, nDOFS);
K2 = sparse(fem.Ci(:), fem.Cj(:), ...
    asColumn(bsxfun(@times, fem.Ke2s, p2)),nDOFS, nDOFS);
K = K1 + K2;

% Start with boundary conditions. 
if fem.bc_option==1
    tmp_aaa = dirichlet_bc(K,mesh);
elseif fem.bc_option==2
    tmp_aaa = mixed_boundary_conditions(K,mesh,fem,kit);
end
%%%%%%%%%%%%%%%%%%%%function end!!666%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%
% ********************************************************************
% *                          homogeneous dirichlet BC( )             *
% *------------------------------------------------------------------*
function tmp_aaa = dirichlet_bc(k,mesh)
% Apply the essential boundary conditions (Dirichlet)
% -------------------------------------------------------------------%
% 源项
num_elecs = length(mesh.mes_node);
f=zeros(mesh.numNode,num_elecs);
current=1;
% 这块应该参照徐世浙书中求解dirichlet边界时的求法
for t=1:num_elecs
    source=mesh.mes_node(t);
    f(source,t)=f(source,t)+current;
end
% 构成 如下形式的矩阵
%                 node_bc
%            k11 ... 0 ... k1n
%            ... ... 0 ... ...
% node_bc    0.. ... 1 ... ..0
%            ... ... 0 ... ...
%            kn1 ... 0 ... knn
k(mesh.bdNodeId,:)=0;
k(:,mesh.bdNodeId)=0;
Ind = sub2ind(size(k),mesh.bdNodeId,mesh.bdNodeId);
k(Ind)=1; %1e10;
tmp_aaa=k\f;
tmp_aaa=full(tmp_aaa)';
%%%%%%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%
end
%%
%**************************************************************************
%*                 mixed_boundary_conditions( )                           *
%*------------------------------------------------------------------------*
function tmp_aaa = mixed_boundary_conditions(K,mesh,fem,kit)
% 这种情况下，不同的电源有不同的刚度矩阵，所以对不同电源点逐个计算其求解。
% For each source, I need the original k, before the bc
num_elecs = length(mesh.mes_node);
tmp_aaa=zeros(mesh.numNode,num_elecs);
% 源项
current=1;
% 这块应该参照徐世浙书中求解dirichlet边界时的求法
sig = 1/mesh.rho(end);
nDofs = mesh.numNode;
for t=1:num_elecs
    src=mesh.mes_node(t);
    srcx = mesh.node(1,src); srcz = mesh.node(2,src);
    f=zeros(mesh.numNode,1);
    n1 = mesh.edge2node(1,mesh.bdEdgeId);
    n2 = mesh.edge2node(2,mesh.bdEdgeId);
    bxc = 0.5*(mesh.node(1,n1)+mesh.node(1,n2));
    bzc = 0.5*(mesh.node(2,n1)+mesh.node(2,n2));
    rx = abs(bxc-srcx); rz = abs(bzc-srcz);
    r = sqrt( rx.^2 + rz.^2);
    cos = max([rx./r;rz./r],[],1);
    kr = fem.k(kit)*r;
    if kr>600,kr=600;end
    bk1 = besselk(1,kr);
    bk0 = besselk(0,kr);
    gg  = (sig*fem.k(kit))*(cos.*bk1./bk0);
    gg  = reshape(gg,1,1,[]);
    L   = reshape(bsxfun(@times, fem.Le, gg),[],1);
    KL  = sparse(fem.BCi(:),fem.BCj(:),L,nDofs,nDofs); 
    k = K+KL;
    f(src)=f(src)+current;
    tmp_aaa(:,t)=k\f;    
end
tmp_aaa=full(tmp_aaa)';
%%%%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%%
end

%********************************************************************
%*                          fast_JAC_CONTROL( )                     *
%*------------------------------------------------------------------*
%* This function controls the calculation of the jacobian .         *
%* it calls a series of functions in order to calculate the         *
%* transformed jacobian                                             *
%******************************************************************** 
function jac = jacob_calc_kit(kit,input,model,mesh,aaa,fem)
% 波数域下的jacobian矩阵，大小为input.num_mes*model.num_param
ke1s = fem.Ke1s(:,:,model.p2c);
ke2s = fem.k(kit)^2*fem.Ke2s(:,:,model.p2c);
Kes = ke1s+ke2s;
tri2node = mesh.tri2node(:,model.p2c);
N = length(model.p2c);
jac=zeros(input.num_mes,N);

pa = input.array_pair(:,1);
pb = input.array_pair(:,2);
pm = input.array_pair(:,3);
pn = input.array_pair(:,4);

for mes=1:input.num_mes  %每个实测数据
    % 这种求解雅可比矩阵的方法为伴随方程法
    src=pa(mes); 
    if src~=0
        Ua = reshape(aaa(src,tri2node),1,3,N);
        Uaa = repmat(Ua,[3,1,1]);
        c_ka = squeeze(sum(bsxfun(@times,Kes,Uaa),2));
    else
        c_ka = zeros(3,N);
    end
    rec = pb(mes);
    if rec~=0
        Ub = reshape(aaa(rec,tri2node),1,3,N);
        Ubb= repmat(Ub,[3,1,1]); 
        c_kb = squeeze(sum(bsxfun(@times,Kes,Ubb),2));
    else
        c_kb = zeros(3,N);
    end
    im = pm(mes);
    if im~=0
        Um = reshape(aaa(im,tri2node),3,N);
    else
        Um = zeros(3,model.num_param);
    end
    in = pn(mes);
    if im~=0
        Un = reshape(aaa(in,tri2node),3,N); 
    else
        Un = zeros(3,N);
    end
    jactmp = sum((c_ka-c_kb).*(Um-Un),1);  
    jac(mes,:) = input.geofac(mes)*reshape(jactmp,1,[]);    
end
%%%%%%% function end%%%%%%%%%%%%%%%%%%%%%%%%%    
end
%%
function fem = homog_field(input,model,mesh,fem)
%------------------------------------------------------------
% 求解均匀半空间的 解析解以及波数域下的解析解 
% current=1;  默认电流=1A  默认网格电阻率=1
% 将电阻率设置为1
%------------------------------------------------------------
elec_x = input.elec_xz(:,1); elec_z = input.elec_xz(:,2);
fem.Ap = zeros(input.num_elec,mesh.numNode);  
rho = model.mean_res;   %以平均电阻率作为背景电阻率
%-------------------------------------------------------------
% Analytic solution for V=I*rho/2pi*r
node_x = mesh.node(1,:); node_z = mesh.node(2,:);
for i=1:length(elec_x)
    %  时间域下的电位
    tmp11=sqrt((elec_x(i)-node_x).*(elec_x(i)-node_x)+ ...
              (elec_z(i)-node_z).*(elec_z(i)-node_z) );
    tmp12=sqrt((elec_x(i)-node_x).*(elec_x(i)-node_x)+ ...
               (-elec_z(i)-node_z).*(-elec_z(i)-node_z) );
    fem.Ap(i,:)=rho./(4*pi).*(1./tmp11 + 1./tmp12);                        
end

% singularity removal 去除电源处的 奇异值
for i=1:input.num_elec
    j=mesh.mes_node(i);
    nums=[1:j-1 j+1:mesh.numNode];
    tmp=fem.Ap(i,nums);
    fem.Ap(i,j)=1.5*max(tmp);
end
%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%
end





