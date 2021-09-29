function final = single_inversing_main(input,model,mesh,fem,inv_opts)
%----------------------------------------------------%
%             反演的一些参数设置                      %
%----------------------------------------------------%
model = model_constrains(model);
if inv_opts.refer_model_flag==1
    model.mref = '读取参考模型';
    pause;  return;
else
    model.mref = []; %model.res_param;
end
% 先不考虑激电的反演
final = [];
%----------------------------------------------------%
%             反演迭代                               
%----------------------------------------------------%
% 开始反演
jac_flag = 2;
fem=compute_femAndJ(input,model,mesh,fem,jac_flag,0);
[~,fem,model,inv_opts]=rms(input,model,fem,inv_opts,0);
for itr = 1:inv_opts.itn
    %===============================
    if itr==1
       fprintf('***********START INVERSION ITERATE*************\n');
    end
    fprintf('itr number:%d\t Lagran:%.2f\n',itr,inv_opts.lagran);
    if itr>1 && inv_opts.beta~=0
        %input = robust_constrain_data(input,fem);
        [model] = robust_constrain(model); 
    end
    %====================================================
    [model] = GN_inverse(input,model,fem,inv_opts);
    [model] = update_prop(input,model,mesh,fem);
    inv_opts.lagran = inv_opts.lagran*0.75;  % 是否有别的方法调节？
    fem=compute_femAndJ(input,model,mesh,fem,jac_flag,itr);
    [ex,fem,model,inv_opts]=rms(input,model,fem,inv_opts,itr);
    %===============================
    final=info_save(input,model,fem,inv_opts,itr,final);
    if ex==2,break;end
    if itr==inv_opts.itn 
       fprintf('***********PROGRAME TIMINAL************\n');
       break;
    end
end
% plot_inv_result(model);
%-----------------------------------------------------
% Here save outputs...
rootdir = pwd;
n = datenum(now);

tmppath = ['\results\inv_results',num2str(n),'.mat'];
filepath = [rootdir,tmppath];
save(filepath,'final');
%%%%%%%%%%function end%%%%%%%%%%%%%%%%%%%
end

function [model] = GN_inverse(input,model,fem,inv_opts)
%---------------------------------------------------------------------
ctc = model.wx'*model.Rx*model.wx+model.wz'*model.Rz*model.wz;
ctc = inv_opts.lagran*ctc;
alpha_s = model.alp_s;%0.001*(90/36)^2;  %DCIP2d
if inv_opts.refer_model_flag==1
    csc = inv_opts.lagran*alpha_s*model.ws'*model.Rs*model.cs;
    b3 = csc*(log(model.res_param(:))-log(model.mref(:)));
else
    [I,J] = size(ctc);
    csc = sparse(I,J);
    b3 = sparse(I,1);
end
% 左端项
WdW = input.Wd'*input.Rd*input.Wd;
JTJ=fem.array_jacobian'*WdW*fem.array_jacobian;
A=(JTJ + ctc+csc);
%-------------------------------------------------------------------
% 右端项
grad_fd = fem.array_jacobian'*WdW*(log(fem.app_res_itr)...
    -log(input.real_data));
grad_fm = ctc*log(model.res_param(:))+b3;
b = -(grad_fd+grad_fm);
%----------------------------------------------------------------
dx1=A\b;
model.dm = dx1;
%-------------------------------------
end

function plot_inv_result(model) 
xc = mean(model.param_xz(:,1:2:8),2);
zc = mean(model.param_xz(:,2:2:8),2);
F = TriScatteredInterp(xc,zc,model.res_param);
res = F(model.grid_x,model.grid_z);
contour(model.grid_x,model.grid_z,res)
axes1 = gca;
rhomin = min(model.res_param); 
rhomax = max(model.res_param);
set(axes1,'CLim',[rhomin, rhomax]);
colormap('jet');
axis (axes1,'equal');
min_x = min(model.grid_x(:));
max_x = max(model.grid_x(:));
xlim(axes1,[min_x max_x]);
min_z = min(model.grid_z(:));
max_z = max(model.grid_z(:));
ylim(axes1,[min_z max_z]);
% set(gca,'YTick',[-12,-8,-4]);
grid on

end











