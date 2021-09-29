function [model]=inverse_single(input,model,fem,inv_opts)
% inv_flag
% 1--Levenberg-Marquardt， 2--Occam, 3--Gauss-Newton，          
switch inv_opts.inv_flag
    case 1 
        [model] = LM_inverse(input,model,fem,inv_opts);
    case 2
        [model] = Occam(input,model,fem,inv_opts);
    case 3
        [model] = GN_inverse(input,model,fem,inv_opts);   
    case 5
        disp('Not yet support')
        % [model] = NLGC(input,model,fem,inv_opts,itr);
end
%%%%%%% function end%%%%%%%%%%%%%%%%%%%5
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
% b = fem.array_jacobian'*WdW*(log(input.real_data)...
%     -log(fem.app_res_itr(:)));
% b2 = ctc*log(model.res_param(:));
% 
% b = b-b2-b3;
b = -(grad_fd+grad_fm);
%----------------------------------------------------------------
dx1=A\b;
model.dm = dx1;
% % 更新模型
% model.res_param_itr = model.res_param.*exp(dx1);
%-------------------------------------
end


function [model] = LM_inverse(input,model,fem,inv_opts)
% GN 反演的简化版
% alpha_s = 0.001*(90/36)^2;  %DCIP2d
ctc = inv_opts.lagran*eye(model.num_param,model.num_param);
% 左端项
WdW = input.Wd'*input.Rd*input.Wd;
JTJ=fem.array_jacobian.'*WdW*fem.array_jacobian;
A = JTJ + ctc;
b = fem.array_jacobian'*WdW*(log(input.real_data)...
    -log(fem.app_res_itr(:)));
dx1 = A\b;
model.dx1 = dx1;
% % 更新模型
% model.res_param_itr = model.res_param.*exp(dx1);
end


function [model] = Occam(input,model,fem,inv_opts)
if inv_opts.refer_model_flag==1
    ctc = model.cx'*model.cx + model.cz'*model.cz + model.cs;
else
    ctc = model.cx'*model.cx + model.cz'*model.cz;
end
ctc = inv_opts.lagran*model.wt*ctc;  
% 左端项
WdW = input.Wd'*input.Rd*input.Wd;
JTJ=fem.array_jacobian.'*WdW*fem.array_jacobian;
A = JTJ + ctc;
%----------------------------------------------------------
%右端项
tmp = fem.array_jacobian*log(model.res_param);
dk = log(input.real_data)-log(fem.app_res_itr)+tmp;
b = fem.array_jacobian'*WdW*dk;
%----------------------------------------------------------
dx1 = A\b;
model.dx1 = dx1;
%---------------------------------------------------------------------
% model.res_param_itr = exp(dx1);
%-------------------------------------
end










