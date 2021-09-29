%* ******************************************************************
%*                             RMS( )                               *
%*------------------------------------------------------------------*
%* This function finds the RMS error and ends the procedure         *
%********************************************************************
function [ex,fem,model,inv_opts]=rms(input,model,fem,inv_opts,itr)
%%
% res_param_itr          试探解，    res_param 保存符合条件的解
% rms_sum_itr            试探解，    rms_sum   保存符合条件的
% app_res_itr            试探解，    app_res   保存符合条件的
% ex: 0--反演未结束， 1--发散了，结束, 2--收敛了，结束
%=========================================
ex=0;
fem.rms_sum_itr=0; % 当前的拟合误差
% Find RMS error 求拟合差
% RMS % 相对误差 
error = input.Wd*((input.real_data-fem.app_res_itr(:)));%./input.real_data; 
fem.rms_sum_itr=sqrt(sum(error.^2)/input.num_mes); %*100;

% Update for first iteration 用来对第一次迭代计算更新
if (itr==0)
    fem.rms_sum=fem.rms_sum_itr;
end
% fprintf('itr number:%d\tRMS=>%f\n',itr,fem.rms_sum_itr);
fprintf('RMS=>%f\n',fem.rms_sum_itr);
%====================================
%----截止条件1：
%----如果拟合差下降速率<预设的速率（即拟合差已经收敛），反演结束
% 默认设置，反演收敛速率为2%
tmp11=abs((fem.rms_sum-fem.rms_sum_itr)/fem.rms_sum);
if (tmp11<0.02 && itr>1)
    disp('Divergence rate is small.');
    ex = 1;
end
% 判断迭代的效果（反演是否收敛），以决定是否结束反演
%----截止条件2：
%----如果当前的拟合差>前一次的拟合差，且itr>1 则说明发散了，结束反演
if (fem.rms_sum_itr>fem.rms_sum && itr>0) 
    fprintf('Divergence occured:\n\t Old RMS=>%f New Rms=>%f\n',fem.rms_sum,fem.rms_sum_itr);
    fem.rms_sum_itr=fem.rms_sum;
    ex = 2;
end

%----截止条件3：
%----如果当前的拟合差<1，结束反演
% if (fem.rms_sum_itr<1) 
%     fprintf('Rms is too small:=>%f\n',fem.rms_sum_itr);
%     fem.rms_sum_itr=fem.rms_sum;
%     ex = 3;
% end
% 将试探解，保存

% 如果不符合上述条件，则反演过程没有结束，
% if ismember(ex,[0,1,3])
% 将试探解，保存
model.res_param=model.res_param_itr;
fem.rms_sum=fem.rms_sum_itr;
fem.app_res=fem.app_res_itr;

% 将试探解，保存

%%%%%%%%%%%%%%%% function end%%%%%%%%%%%%%%%%%%%%%
end


