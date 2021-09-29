function [model] = robust_constrain(model)
% cut_off = 1e-7;
m = log(model.res_param);
rc = abs(model.wc*m);
cut_off = mean(rc);
Rc = 2*ones(size(rc));
indx = (rc>cut_off);
Rc(indx) = 2*cut_off./rc(indx);
model.Rc = sparse(diag(0.5*Rc));

if ~isempty(model.mref)
    mref = log(model.mref);
    rs = abs(model.ws*(m-mref));
    Rs = 2*ones(size(rs));
    inds = find(rs>cut_off);
    Rs(inds) = 2*cut_off./rs(inds);
    model.Rs = diag(0.5*Rs); 
end

end

