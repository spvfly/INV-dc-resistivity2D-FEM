clear all; clc;
addpath forward0 invers0
input = read_data();
model = create_model_grid(input);
[model,mesh] = create_mesh_mix(input,model);
%========================================================
fem = pdetrgm(mesh); % 单元分析
fem.k = [0.0031677 0.0301330 0.1285886 0.4599185 1.5842125]; % 波数
fem.g = [0.0067253 0.0314373 0.1090454 0.3609340 1.3039204];
fem.bc_option = 1;  %采用第三类边界
jacobian_flag = 2;
%======================================================
inv_opts.lagran = 0.15;
inv_opts.itn = 5;
inv_opts.beta = 0;
inv_opts.refer_model_flag=0;
final = single_inversing_main(input,model,mesh,fem,inv_opts);
plot_res_inv(final);
rmpath forward0 invers0