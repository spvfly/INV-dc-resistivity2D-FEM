function [input]=read_data()
% 先默认 无ip/sip数据
% 注意，当使用任意源发射-任意接收 需要把 源和接收分开
% 而高密度的阵列可以统一来索引就行--该脚本适合这种情况
input.ip_flag = 0;
input.sip_flag = 0;
[filename,pathname] = uigetfile({'*.d'},'select a data file',...
                        'MultiSelect','off');
input.mes_in = [pathname,filename];                    
% First check for Loke
tmp_d=importdata(input.mes_in);
[r,~] = size(tmp_d);
input.num_mes=r;
% This is the common part
input.ax=tmp_d(:,1);
input.az=abs(tmp_d(:,2));
input.bx=(tmp_d(:,3));
input.bz=abs(tmp_d(:,4));
input.mx=tmp_d(:,5);
input.mz=abs(tmp_d(:,6));
input.nx=tmp_d(:,7);
input.nz=abs(tmp_d(:,8));
%--------------------------------
input.real_data=tmp_d(:,9);
%%%IP DATA or SIP%%%%%%
try
    input.ip_data = tmp_d(:,10); 
    if input.sip_flag==1
        disp('SIP data found');
        input.ip_num=1;
        input.real_data=complex(input.real_data,input.ip_data);
        % if are in mag and phase
        % real_part=input.real_data.*(cos(input.ip_data./1000));
        % imag_part=input.real_data.*(sin(input.ip_data./1000));
        % input.real_data=complex(real_part,imag_part);      
    elseif input.ip_flag==1
        disp('IP data found');
        input.ip_num=2;
        input.ip_data=input.ip_data./1000; % Loke has it in msec
    else
        input.stdev_error=input.ip_data;   % Currently on fow DC
        input.ip_num=1;
    end
catch
    disp('No IP or SIP data found and No stdev info!!');
    input.ip_num=1;
    % 数据的标准差，没有怎么办？
    % input.stdev_error=ones(input.num_mes,1);
    % input.stdev_error = (input.real_data.^0.25);????
    input.stdev_error = sqrt(0.03*input.real_data + 0.1);    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here check for IP or SIP and STDEV
if input.sip_flag==1 || input.ip_flag==1
    try
        input.stdev_error=tmp_d(:,11); 
    catch
        input.stdev_error=ones(input.num_mes,1);
        % 极化率，没有测量误差时应该怎么处理呢？？
    end
end
clear tmp_d;
input.Wd=spdiags(1./input.stdev_error,0,input.num_mes,input.num_mes);
input.Rd = sparse(eye(size(input.Wd)));
%===========================================
% Try to find number of electrodes used
%--------------------------------------------------------------------------
% search for how many electrodes
error=0;
b_el=0;
n_el=0;
input.array_type=4;
ax=length(unique(input.ax)); az=length(unique(input.az));
bx=length(unique(input.bx)); bz=length(unique(input.bz)); 
mx=length(unique(input.mx)); mz=length(unique(input.mz));
nx=length(unique(input.nx)); nz=length(unique(input.nz));
% ---------check----------------------
if ax==1 && az==1
    disp('Error. A electrode must be there ALWAYS');
    error=1;
end
if mx==1 && mz==1
    disp('Error. M electrode must be there ALWAYS');
    error=1;
end

if bx==1 && bz==1  % b electrode is not exit;
    b_el=1;
    error=0;
end

if nx==1 && nz==1  % n electrode is not exit;
    n_el=1;        
    error=0;
end

if b_el==1 && n_el==0
   disp('Pole-Dipole array:AMN');
   input.array_type=3;
elseif b_el==1 && n_el==1
   disp('Pole-Pole array:AM');
   input.array_type=2;
elseif b_el==0 && n_el==1
   disp('Not tested array');
   pause
else
   disp('4 electrode array');
   input.array_type=4;
end
if error==1,pause;end  
% input.m_array_type(1:input.num_mes,1)=input.array_type;
input = calc_array_pair(input);
input = geom_factor(input);
end

function input = calc_array_pair(input)
all_probe = [input.ax, input.az; input.bx, input.bz;...
             input.mx, input.mz; input.nx, input.nz];
all_probe = unique(all_probe,'rows');
% 电极的索引
array_pair = zeros(input.num_mes,4);
elec_x = all_probe(:,1); elec_z = all_probe(:,2);
n_e = input.array_type;
for i = 1:input.num_mes
    if n_e==4  
        par1 = intersect(find(elec_x==input.ax(i)),find(elec_z==input.az(i)));
        par2 = intersect(find(elec_x==input.bx(i)),find(elec_z==input.bz(i)));
        par3 = intersect(find(elec_x==input.mx(i)),find(elec_z==input.mz(i)));
        par4 = intersect(find(elec_x==input.nx(i)),find(elec_z==input.nz(i)));
        array_pair(i,:) =[par1,par2,par3,par4];
    elseif n_e==3
        par1 = intersect(find(elec_x==input.ax(i)),find(elec_z==input.az(i)));
        par2 = 0;
        par3 = intersect(find(elec_x==input.mx(i)),find(elec_z==input.mz(i)));
        par4 = intersect(find(elec_x==input.nx(i)),find(elec_z==input.nz(i)));
        array_pair(i,:) = [par1,par2,par3,par4];
    elseif n_e==2
        par1 = intersect(find(elec_x==input.ax(i)),find(elec_z==input.az(i)));
        par2 = 0;
        par3 = intersect(find(elec_x==input.mx(i)),find(elec_z==input.mz(i)));
        par4 = 0;
        array_pair(i,:) = [par1,par2,par3,par4];
    end
end
input.array_pair = array_pair;
%  整理电极坐标，地下电极z坐标为负
input.elec_xz = all_probe;
I = (input.elec_xz(:,2)>0);
input.elec_xz(I,2) = -input.elec_xz(I,2);
input.num_elec = length(all_probe);

ind = (input.elec_xz(:,2)==0);
input.surf_elec = input.elec_xz(ind,:);
ind = (input.elec_xz(:,2)~=0);
input.bore_elec = input.elec_xz(ind,:);

end

function input = geom_factor(input)
% 只适用于水平地表时，计算 几何因子
VAM=zeros(input.num_mes,1);  VAN=zeros(input.num_mes,1);
VBM=zeros(input.num_mes,1);  VBN=zeros(input.num_mes,1);
%--------------------------------------------------------
n_e = unique(input.array_type);
if(n_e==4)
    %/* flag for 4 probe */
    aflag=1; bflag=1; mflag=1; nflag=1;
elseif(n_e==3)
    %/* flag for 3 probe */
    aflag=1; bflag=0; mflag=1; nflag=1;
elseif(n_e==2)
    %/* flag for 2 probe */
    aflag=1; mflag=1; bflag=0; nflag=0;
end
%-------------------------------------------------------
dd = [];
if aflag*mflag~=0
    r1 = sqrt( (input.ax-input.mx).^2 + (input.az-input.mz).^2  );
    r2 = sqrt((input.ax-input.mx).^2 + (-input.az-input.mz).^2 );
    VAM = 1./r1+1./r2;
    dd = [dd min(r1) max(r1)];
end
if aflag*nflag~=0
    r1 = sqrt( (input.ax-input.nx).^2 + (input.az-input.nz).^2  );
    r2 = sqrt((input.ax-input.nx).^2 + (-input.az-input.nz).^2 );
    VAN = 1./r1+1./r2;
    dd = [dd min(r1) max(r1)];
end
if bflag*mflag~=0
    r1 = sqrt( (input.bx-input.mx).^2 + (input.bz-input.mz).^2  );
    r2 = sqrt( (input.bx-input.mx).^2 + (-input.bz-input.mz).^2 );
    VBM = 1./r1+1./r2;
    dd = [dd min(r1) max(r1)];
end
if bflag*nflag~=0
    r1 = sqrt( (input.bx-input.nx).^2 + (input.bz-input.nz).^2  );
    r2 = sqrt((input.bx-input.nx).^2 + (-input.bz-input.nz).^2 );
    VBN = 1./r1+1./r2;
    dd = [dd min(r1) max(r1)];
end
input.unit_a = min(dd);
input.n = max(dd)/min(dd);
input.geofac=(4*pi)./(VAM-VBM-VAN+VBN);
end

