function [model]=update_prop(input,model,mesh,fem)
%------------更新电性参数---------------------------------
f = log(input.real_data);
WdW = input.Wd'*input.Rd*input.Wd;
f0 = log(fem.app_res_itr); % f(m0)
% --------------------------------
% m(k+1) = m(k) + alp*dm;
% --------------------------------
% m1
model.res_param_itr = model.res_param.*exp(0.5*model.dm);
fem=compute_femAndJ(input,model,mesh,fem,0,1); % 正演 但不求解jacob矩阵
f1 = log(fem.app_res_itr); %f(m1)

% m2
model.res_param_itr = model.res_param.*exp(model.dm);
fem=compute_femAndJ(input,model,mesh,fem,0,1); % 正演 但不求解jacob矩阵
f2 = log(fem.app_res_itr); %f(m2)

%=======================================
% fhi(m)
phi0 = (f-f0)'*WdW*(f-f0);
phi1 = (f-f1)'*WdW*(f-f1);
phi2 = (f-f2)'*WdW*(f-f2);

tmp1 = -0.75*phi0 +phi1 -0.25*phi2;
tmp2 = -0.5*phi0 + phi1  -0.5*phi2;
alp = min(0.5*tmp1/tmp2,1);
alp = max(0.05,alp);
if alp==0.05,disp(alp);end
model.res_param_itr = model.res_param.*exp(alp*model.dm);
%----------------------------------------------
% 检查是否有超过所设界限的 电阻率值，并限制住它
mean_res = model.mean_res;
limit_max_res = mean_res*5;
limit_min_res = mean_res/5;
ind = find(model.res_param_itr>limit_max_res);
if ~isempty(ind)
   model.res_param_itr(ind)=limit_max_res;
end
ind = find(model.res_param_itr<limit_min_res);
if ~isempty(ind)
    model.res_param_itr(ind)=limit_min_res;
end
%%%%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%%%
end

% % 黄金分割法搜索
% % m(k+1) = m(k) + alp*dm;
% m1 = model.res_param.*exp(model.dm);
% % 更新正演网格的电阻率参数
% for i=1:model.num_param
%     ind = (mesh.c2p==i);
%     mesh.prop(ind)=m1(i);
% end
% % 正演 但不求解jacob矩阵
% fem=compute_femAndJ(input,model,mesh,fem,0,1);
% f1 = (fem.app_res_itr);
% f = (input.real_data);
% % 初值
% a = 0; b=1;
% x = a+0.382*(b-a);
% y = a+0.618*(b-a);
% fx = (1-x)*f0+x*f1; phix = (f-fx)'*WdW*(f-fx);
% fy = (1-y)*f0+y*f1; phiy = (f-fy)'*WdW*(f-fy);
% while 1  
%    if phix<phiy
%        b=y; y=x; 
%        x = a+0.382*(b-a); fx = (1-x)*f0+x*f1; 
%        phix = (f-fx)'*WdW*(f-fx);
%    else
%        a = x; x = y;
%        y = a+0.618*(b-a);
%        fy = (1-y)*f0+y*f1; phiy = (f-fy)'*WdW*(f-fy); 
%    end
%    if abs(a-b)<1e-1,alp = 0.5*(a+b);break;end
% end
% disp(alp)
% model.res_param_itr = model.res_param.*exp(alp*model.dm);