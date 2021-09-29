function [model] = robust_constrain(model)
% cut_off = 1e-7;
% ��Ϊfhi_dΪL2norm,��Ҫ��fhi_mΪLp��p=(0,1)
% ���� �˴�Ӧ����һ��ϵ��0.5  �� L1/L2 = 1/2 
m = log10(model.res_param_itr);
rx = abs(model.wx*m);
cut_off = mean(rx);
Rx = 2*ones(size(rx));
indx = (rx>cut_off);
Rx(indx) = 2*cut_off./rx(indx);
Rx = 0.5*Rx;
model.Rx = diag(Rx);
% 
rz = abs(model.wz*m);
cut_off = mean(rz);
Rz = 2*ones(size(rz));
indz= (rz>cut_off);
Rz(indz) = 2*cut_off./rz(indz);
Rz = 0.5*Rz;
model.Rz = diag(Rz);

try
   mref = log10(model.mref);
   rs = abs(model.ws*(m-mref));
   Rs = 2*ones(size(rs));
   inds = find(rs>cut_off);
   Rs(inds) = 2*cut_off./rs(inds);
   model.Rs = diag(Rs); 
end

end

