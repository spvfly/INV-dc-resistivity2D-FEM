function fem=compute_femAndJ(input,model,mesh,fem,jc_flag,itr)
%-------------------------------------------
% ������������ĵ����ʲ���
for i = 1:model.num_param
    mesh.rho(mesh.c2p==i) = model.res_param_itr(i);
end
num_mes = input.num_mes;
tmp_jac=zeros(num_mes,length(model.p2c));
num_elec = length(mesh.mes_node);
aaa=zeros(num_elec,mesh.numNode);
%----------����Ԫ����+��ȫ��jacob�������----------------
for kit=1:length(fem.k)
    %-----------��������
    h2 = waitbar(kit/length(fem.k));
    drawnow;
    %------��ɸնȾ���������Ԫ��
    [tmp_aaa]=fem_calc_kit(kit,mesh,fem);
    fem.tmp_aaa{kit}=tmp_aaa;
    %-----Inverse fourier transform ʱ���������Ԫ��
    aaa = aaa+fem.g(kit).*tmp_aaa/2; %/pi;  
    %---��Jacobian����(���ǰ������ݹ�ʽ���)
    if jc_flag==2 || itr==0 
       trans_jac1=jacob_calc_kit(kit,input,model,mesh,fem,tmp_aaa);
       % Inverse fourier transform
       tmp_jac=tmp_jac+fem.g(kit)*trans_jac1/2; %/pi;
    end
end
fem.aaa = aaa;  % ��ֵģ��⣨��λ����
close(h2);
%-----------------------------------------------
app_res_itr = array_app_res(input,mesh,aaa);
% ������Ȱ�ռ�Ľ�����
if itr==0
    fem = homog_field(input,model,mesh,fem);
    app_res_ana = array_app_res(input,mesh,fem.Ap);
    % �ӵ����ʵ�������ϵ��
    fem.normalized = app_res_ana./app_res_itr; %%
end
fem.app_res_itr = app_res_itr;
% -------jacob==>log(jacob)----------------
array_jacobian = zeros(num_mes,model.num_param);
if jc_flag==2 || itr==0
    c2p = nonzeros(mesh.c2p);
    for i = 1:num_mes
        jac = accumarray(c2p,tmp_jac(i,:))';
        rhoa = fem.app_res_itr(i);    % �ӵ�����
        sig  = 1./model.res_param_itr';  % ģ�͵����ʣ��棩
        array_jacobian(i,:) = jac.*(sig./rhoa);
    end  
end
fem.array_jacobian = array_jacobian;
if(jc_flag==1 && itr>0)
    fem=quasi_newton(model,fem);
end

%%%%%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%
function fem=quasi_newton(model,fem)
    dd = model.dx1'*model.dx1;
    dm=log10(fem.app_res_itr) - log10(fem.app_res) - fem.array_jacobian*model.dx1;
    fem.array_jacobian = fem.array_jacobian + dm*model.dx1'/dd;
end


%%
% *********************************************************************
% *   fem_calc( )��������նȾ��� �� ������Է�����                    *
% *-------------------------------------------------------------------*
% * This function forms the global stiffness matrix.                  *
% *                                                                   *
% * This function calculate the stiffnes term for every element.      *
% * It is valid only for linear triangular elements                   *
% *********************************************************************/
function [tmp_aaa]=fem_calc_kit(kit,mesh,fem)
% Start to assemble all area elements in Omega 
sig = 1./mesh.rho;
%����Ԫ��һ��ϵ��
p1 = sig./(4*fem.area); 
% ����Ԫ�ڶ���ϵ��
a = sig.*(1/6)*fem.k(kit).^2.*fem.area;
b = a/2;
%
k11 = p1.*(fem.b1.^2 + fem.c1.^2) + a;
k12 = p1.*(fem.b1.*fem.b2 + fem.c1.*fem.c2)+b;
k13 = p1.*(fem.b1.*fem.b3 + fem.c1.*fem.c3)+b;
k22 = p1.*(fem.b2.^2 + fem.c2.^2)+a;
k23 = p1.*(fem.b2.*fem.b3+ fem.c2.*fem.c3)+b;
k33 = p1.*(fem.b3.^2 + fem.c3.^2)+a;  
k21 = k12;
k31 = k13;
k32 = k23;

% �ܸնȾ���
kkk = sparse([k11,k12,k13,k21,k22,k23,k31,k32,k33]);
K = accumarray(fem.kIndx,kkk);
% Start with boundary conditions. 
% homogeneous dirichlet or mixed bounary conditions
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
% ------------------------------------------------------------------------%
num_elecs = length(mesh.mes_node);
f=zeros(mesh.numNode,num_elecs);
current=1;
% ���Ӧ�ò����������������dirichlet�߽�ʱ����
for t=1:num_elecs
    source=mesh.mes_node(t);
    f(source,t)=f(source,t)+current;
end
%                 node_bc
%            k11 ... 0 ... k1n
%            ... ... 0 ... ...
% node_bc    0.. ... 1 ... ..0
%            ... ... 0 ... ...
%            kn1 ... 0 ... knn
k(mesh.bdNodeId,:)=0;
k(:,mesh.bdNodeId)=0;
Ind = sub2ind(size(k),mesh.bdNodeId,mesh.bdNodeId);
k(Ind)=1;  % k(Ind)=1e10; 
tmp_aaa=k\f;
tmp_aaa=full(tmp_aaa)';
%%%%%%%%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%
end

%%
%**************************************************************************
%*                 mixed_boundary_conditions( )                           *
%*------------------------------------------------------------------------*
function tmp_aaa = mixed_boundary_conditions(K,mesh,fem,kit)
% ��������£���ͬ�ĵ�Դ�в�ͬ�ĸնȾ������ԶԲ�ͬ��Դ�������������⡣
% For each source, I need the original k, before the bc
num_elecs = length(mesh.mes_node);
tmp_aaa=zeros(mesh.numNode,num_elecs);
% Դ��
current=1;
% ���Ӧ�ò����������������dirichlet�߽�ʱ����
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
%%
%********************************************************************
%*                          fast_JAC_CONTROL( )                     *
%*------------------------------------------------------------------*
%* This function controls the calculation of the jacobian .         *
%* it calls a series of functions in order to calculate the         *
%* transformed jacobian                                             *
%******************************************************************** 
function jac = jacob_calc_kit(kit,input,model,mesh,fem,aaa)
% �������µ�jacobian���󣬴�СΪnum_mes*model.num_param
jam=0; jan=0; jbm=0; jbn=0;
num_mes = input.num_mes;
jac=zeros(num_mes,length(model.p2c));
p2c = model.p2c;
a1=mesh.tri2node(1,p2c);    %I 
a2=mesh.tri2node(2,p2c);    %J
a3=mesh.tri2node(3,p2c);    %L

b1 = fem.b1(p2c);
b2 = fem.b2(p2c);
b3 = fem.b3(p2c);
c1 = fem.c1(p2c);
c2 = fem.c2(p2c);
c3 = fem.c3(p2c);
area = fem.area(p2c);


for mes=1:num_mes  % each data
    pa = input.array_pair(mes,1); pb = input.array_pair(mes,2);
    pm = input.array_pair(mes,3); pn = input.array_pair(mes,4);
    % Adjoint Green method
    if((pa*pm)~=0)
        src=pa; rec=pm;
        s1=aaa(src,a1)+aaa(src,a2)+aaa(src,a3);
        s2=aaa(rec,a1)+aaa(rec,a2)+aaa(rec,a3);
        sum1=s1.*s2+(aaa(src,a1).*aaa(rec,a1)+aaa(src,a2).*aaa(rec,a2)+aaa(src,a3).*aaa(rec,a3));  
        sum1 = fem.k(kit)^2.*area.*sum1/12;
        
        x_flux1=aaa(src,a1).*b1+aaa(src,a2).*b2+aaa(src,a3).*b3;
        x_flux2=aaa(rec,a1).*b1+aaa(rec,a2).*b2+aaa(rec,a3).*b3;
        z_flux1=aaa(src,a1).*c1+aaa(src,a2).*c2+aaa(src,a3).*c3;
        z_flux2=aaa(rec,a1).*c1+aaa(rec,a2).*c2+aaa(rec,a3).*c3;
        sum2 = (x_flux1.*x_flux2+z_flux1.*z_flux2)./(4.*area);
        jam=sum2 + sum1;
    end
    if((pa*pn)~=0)
        src=pa; rec=pn;
        s1=aaa(src,a1)+aaa(src,a2)+aaa(src,a3);
        s2=aaa(rec,a1)+aaa(rec,a2)+aaa(rec,a3);
        sum1=s1.*s2+(aaa(src,a1).*aaa(rec,a1)+aaa(src,a2).*aaa(rec,a2)+aaa(src,a3).*aaa(rec,a3));  
        sum1 = fem.k(kit)^2.*area.*sum1/12;
        
        x_flux1=aaa(src,a1).*b1+aaa(src,a2).*b2+aaa(src,a3).*b3;
        x_flux2=aaa(rec,a1).*b1+aaa(rec,a2).*b2+aaa(rec,a3).*b3;
        z_flux1=aaa(src,a1).*c1+aaa(src,a2).*c2+aaa(src,a3).*c3;
        z_flux2=aaa(rec,a1).*c1+aaa(rec,a2).*c2+aaa(rec,a3).*c3;
        sum2 = (x_flux1.*x_flux2+z_flux1.*z_flux2)./(4.*area);
        jan=sum2 + sum1;
    end
    if((pb*pm)~=0)
        src=pb; rec=pm;
        s1=aaa(src,a1)+aaa(src,a2)+aaa(src,a3);
        s2=aaa(rec,a1)+aaa(rec,a2)+aaa(rec,a3);
        sum1=s1.*s2+(aaa(src,a1).*aaa(rec,a1)+aaa(src,a2).*aaa(rec,a2)+aaa(src,a3).*aaa(rec,a3));  
        sum1 = fem.k(kit)^2.*area.*sum1/12;
        
        x_flux1=aaa(src,a1).*b1+aaa(src,a2).*b2+aaa(src,a3).*b3;
        x_flux2=aaa(rec,a1).*b1+aaa(rec,a2).*b2+aaa(rec,a3).*b3;
        z_flux1=aaa(src,a1).*c1+aaa(src,a2).*c2+aaa(src,a3).*c3;
        z_flux2=aaa(rec,a1).*c1+aaa(rec,a2).*c2+aaa(rec,a3).*c3;
        sum2 = (x_flux1.*x_flux2+z_flux1.*z_flux2)./(4.*area);
        jbm=sum2 + sum1;
    end
    if((pb*pn)~=0)
        src=pb; rec=pn;
        s1=aaa(src,a1)+aaa(src,a2)+aaa(src,a3);
        s2=aaa(rec,a1)+aaa(rec,a2)+aaa(rec,a3);
        sum1=s1.*s2+(aaa(src,a1).*aaa(rec,a1)+aaa(src,a2).*aaa(rec,a2)+aaa(src,a3).*aaa(rec,a3));  
        sum1 = fem.k(kit)^2.*area.*sum1/12;
        
        x_flux1=aaa(src,a1).*b1+aaa(src,a2).*b2+aaa(src,a3).*b3;
        x_flux2=aaa(rec,a1).*b1+aaa(rec,a2).*b2+aaa(rec,a3).*b3;
        z_flux1=aaa(src,a1).*c1+aaa(src,a2).*c2+aaa(src,a3).*c3;
        z_flux2=aaa(rec,a1).*c1+aaa(rec,a2).*c2+aaa(rec,a3).*c3;
        sum2 = (x_flux1.*x_flux2+z_flux1.*z_flux2)./(4.*area);
        jbn=sum2 + sum1;
    end
    jac(mes,:) = input.geofac(mes)*(jam-jbm-jan+jbn);
end

%%%%%%% function end%%%%%%%%%%%%%%%%%%%%%%%%%    
end

%%
function fem = homog_field(input,model,mesh,fem)
%------------------------------------------------------------
% �����Ȱ�ռ�� �������Լ��������µĽ����� 
% current=1;  Ĭ�ϵ���=1A  Ĭ�����������=1
% ������������Ϊ1
%------------------------------------------------------------
elec_x = input.elec_xz(:,1); elec_z = input.elec_xz(:,2);
fem.Ap = zeros(input.num_elec,mesh.numNode);  
rho = model.mean_res;   %��ƽ����������Ϊ����������
%-------------------------------------------------------------
% Analytic solution for V=I*rho/2pi*r
node_x = mesh.node(1,:); node_z = mesh.node(2,:);
for i=1:length(elec_x)
    %  ʱ�����µĵ�λ
    tmp11=sqrt((elec_x(i)-node_x).*(elec_x(i)-node_x)+ ...
              (elec_z(i)-node_z).*(elec_z(i)-node_z) );
    tmp12=sqrt((elec_x(i)-node_x).*(elec_x(i)-node_x)+ ...
               (-elec_z(i)-node_z).*(-elec_z(i)-node_z) );
    fem.Ap(i,:)=rho./(4*pi).*(1./tmp11 + 1./tmp12);                        
end

% singularity removal ȥ����Դ���� ����ֵ
for i=1:input.num_elec
    j=mesh.mes_node(i);
    nums=[1:j-1 j+1:mesh.numNode];
    tmp=fem.Ap(i,nums);
    fem.Ap(i,j)=1.5*max(tmp);
end
%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%%%%
end





