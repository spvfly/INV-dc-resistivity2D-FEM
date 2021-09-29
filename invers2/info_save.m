function final=info_save(input,model,fem,inv_opts,itr,final)
if itr==1
    % input相关内容
    final.mes_in=input.mes_in;
    final.num_mes=input.num_mes;
    final.elec_x = input.elec_xz(:,1);
    final.elec_z = input.elec_xz(:,2);
    final.num_param=model.num_param;
    % 模型model 相关内容
    final.param_xz = model.param_xz;
    final.grid_x  = model.grid_x;
    final.grid_z  = model.grid_z;
    final.refer_model_flag=inv_opts.refer_model_flag;
end

final.itn=itr;
final.lagrn=inv_opts.lagran;

% Keep inversion results
final.res_param_itr(:,itr)=model.res_param;
final.array_model_data=fem.app_res_itr;
final.RMS(itr)=fem.rms_sum;
end

    


