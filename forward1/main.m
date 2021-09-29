clear all; clc;
input = read_data();

model = create_model_grid(input);
[model,mesh] = create_mesh_mix(input,model);
% ��Ԫ����
fem = pdetrgm(mesh);

% ����
fem.k = [0.0031677 0.0301330 0.1285886 0.4599185 1.5842125];
fem.g = [0.0067253 0.0314373 0.1090454 0.3609340 1.3039204];

% ----------------------------------
%          ����Ԫ
% ----------------------------------
fem.bc_option = 2;  %���õ�����߽�
jacobian_flag = 2;
fem = compute_femAndJ(input,model,mesh,fem,jacobian_flag,0);