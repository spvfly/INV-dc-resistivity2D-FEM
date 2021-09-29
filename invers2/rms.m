%* ******************************************************************
%*                             RMS( )                               *
%*------------------------------------------------------------------*
%* This function finds the RMS error and ends the procedure         *
%********************************************************************
function [ex,fem,model,inv_opts]=rms(input,model,fem,inv_opts,itr)
%%
% res_param_itr          ��̽�⣬    res_param ������������Ľ�
% rms_sum_itr            ��̽�⣬    rms_sum   �������������
% app_res_itr            ��̽�⣬    app_res   �������������
% ex: 0--����δ������ 1--��ɢ�ˣ�����, 2--�����ˣ�����
%=========================================
ex=0;
fem.rms_sum_itr=0; % ��ǰ��������
% Find RMS error ����ϲ�
% RMS % ������ 
error = input.Wd*((input.real_data-fem.app_res_itr(:)));%./input.real_data; 
fem.rms_sum_itr=sqrt(sum(error.^2)/input.num_mes); %*100;

% Update for first iteration �����Ե�һ�ε����������
if (itr==0)
    fem.rms_sum=fem.rms_sum_itr;
end
% fprintf('itr number:%d\tRMS=>%f\n',itr,fem.rms_sum_itr);
fprintf('RMS=>%f\n',fem.rms_sum_itr);
%====================================
%----��ֹ����1��
%----�����ϲ��½�����<Ԥ������ʣ�����ϲ��Ѿ������������ݽ���
% Ĭ�����ã�������������Ϊ2%
tmp11=abs((fem.rms_sum-fem.rms_sum_itr)/fem.rms_sum);
if (tmp11<0.02 && itr>1)
    disp('Divergence rate is small.');
    ex = 1;
end
% �жϵ�����Ч���������Ƿ����������Ծ����Ƿ��������
%----��ֹ����2��
%----�����ǰ����ϲ�>ǰһ�ε���ϲ��itr>1 ��˵����ɢ�ˣ���������
if (fem.rms_sum_itr>fem.rms_sum && itr>0) 
    fprintf('Divergence occured:\n\t Old RMS=>%f New Rms=>%f\n',fem.rms_sum,fem.rms_sum_itr);
    fem.rms_sum_itr=fem.rms_sum;
    ex = 2;
end

%----��ֹ����3��
%----�����ǰ����ϲ�<1����������
% if (fem.rms_sum_itr<1) 
%     fprintf('Rms is too small:=>%f\n',fem.rms_sum_itr);
%     fem.rms_sum_itr=fem.rms_sum;
%     ex = 3;
% end
% ����̽�⣬����

% ����������������������ݹ���û�н�����
% if ismember(ex,[0,1,3])
% ����̽�⣬����
model.res_param=model.res_param_itr;
fem.rms_sum=fem.rms_sum_itr;
fem.app_res=fem.app_res_itr;

% ����̽�⣬����

%%%%%%%%%%%%%%%% function end%%%%%%%%%%%%%%%%%%%%%
end


