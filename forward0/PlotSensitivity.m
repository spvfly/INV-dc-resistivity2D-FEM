function PlotSensitivity(input,model,fem,pds)
% 不推荐使用 TriScatteredInterp。请改用 scatteredInterpolant。
xc = mean(model.params(:,1:2:8),2); zc = mean(model.params(:,2:2:8),2);
F = TriScatteredInterp(xc,zc,fem.array_jacobian(pds,:)');
qz=F(model.grid_x,model.grid_z);
contourf(model.grid_x, model.grid_z, qz,'LineStyle','none');
%     shading interp
colorbar
axis equal;% off; 
ext = 0.5;
left = min(min(model.grid_x))-0.5; right = max(max(model.grid_x))+0.5;
bott = min(min(model.grid_z))-0.5; top = max(max(model.grid_z))+0.5;
xlim([left,right]); 
ylim([bott,top]);


hold on
plot(input.ax(pds),input.az(pds),'ko','LineWidth',4);
text(input.ax(pds),input.az(pds)+input.unit_a,'A','FontSize',22);

if input.n_e == 4
    plot(input.bx(pds),input.bz(pds),'ko','LineWidth',4);
    text(input.bx(pds),input.bz(pds)+input.unit_a,'B','FontSize',22);
end

plot(input.mx(pds),input.mz(pds),'ko','LineWidth',4);
text(input.mx(pds),input.mz(pds)+input.unit_a,'M','FontSize',22);

if input.n_e==3 || input.n_e==4
    plot(input.nx(pds),input.nz(pds),'ko','LineWidth',4);
    text(input.nx(pds),input.nz(pds)+input.unit_a,'N','FontSize',22);
end

fn = 'Microsoft Sans Serif'; fs = 11;
xlabel('Distance (m)','FontName',fn,'FontSize',fs)
ylabel('Pseudo-Depth (m)','FontName',fn,'FontSize',fs)
title(['Calculated ',' Sensitivity'],'FontName',fn,'FontSize',fs);

end