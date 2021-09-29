function plot_res_inv(final)
xc = final.param_xz(:,1);
zc = final.param_xz(:,2);
rho = log10(final.res_param_itr(:,end));
rhomin = min(rho); rhomax = max(rho);
F = TriScatteredInterp(xc,zc,rho,'natural');
rr = F(final.grid_x,final.grid_z);
contourf(final.grid_x,final.grid_z,rr,10,'EdgeColor','none');
h1 = gca;
set(h1,'CLim',[rhomin, rhomax]);
colormap(h1,'jet');
axis(h1,'equal');
min_x = min(final.grid_x(:));  min_z = min(final.grid_z(:));
max_x = max(final.grid_x(:));  max_z = max(final.grid_z(:));
xlim(h1,[min_x max_x]);
ylim(h1,[min_z max_z]);

% set(gca,'YTick',[-12,-8,-4]);
grid on
end



