function input = robust_constrain_data(input,fem)
rd = input.Wd*(log(input.real_data)-log(fem.app_res_itr(:)));
cut_off = mean(rd);
Rd = 2*ones(size(rd));
indx = (rd>cut_off);
Rd(indx) = 2*cut_off./rd(indx);
% Rx = 2*cut_off^2./(rx.^2+cut_off^2);
input.Rd = diag(Rd);
end