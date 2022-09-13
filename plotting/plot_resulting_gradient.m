function  plot_resulting_gradient(coil_layouts,single_ind_to_plot,plot_title)

pos_data=coil_layouts(single_ind_to_plot).out.target_field.coords;
std_plot_factor=5;

%remove noise points
coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx(abs(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx)>5)=0;
coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy(abs(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy)>5)=0;
coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz(abs(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz)>5)=0;

%[stongest_gradient_val,stongest_gradient]=max([mean(abs(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx),'omitnan'),[mean(abs(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy),'omitnan'),[mean(abs(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz),'omitnan')]]]);
[~,stongest_gradient]=max(abs([coil_layouts(single_ind_to_plot).out.layout_gradient.mean_local_gx_loops/coil_layouts(single_ind_to_plot).out.layout_gradient.std_local_gx_loops ...
    coil_layouts(single_ind_to_plot).out.layout_gradient.mean_local_gy_loops/coil_layouts(single_ind_to_plot).out.layout_gradient.std_local_gy_loops ...
    coil_layouts(single_ind_to_plot).out.layout_gradient.mean_local_gz_loops/coil_layouts(single_ind_to_plot).out.layout_gradient.std_local_gz_loops]));
switch stongest_gradient
    case 1
stongest_gradient_val=coil_layouts(single_ind_to_plot).out.layout_gradient.mean_local_gx_loops;
    case 2
stongest_gradient_val=coil_layouts(single_ind_to_plot).out.layout_gradient.mean_local_gy_loops;
    case 3
stongest_gradient_val=coil_layouts(single_ind_to_plot).out.layout_gradient.mean_local_gz_loops;
end

std_gradients=[std(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx_loops,'omitnan') std(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy_loops,'omitnan') std(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz_loops,'omitnan')];
plot_range=[stongest_gradient_val-std_plot_factor*std_gradients(stongest_gradient),stongest_gradient_val+std_plot_factor*std_gradients(stongest_gradient)];

if std_gradients(stongest_gradient)==0
plot_range(1)=stongest_gradient_val*0.5;
plot_range(2)=stongest_gradient_val*1.5;
end

dot_size=100;

%Plot Histograms of the gradient
figure('name',plot_title);
tiledlayout('flow');


nexttile;
hold on;
title("Gmag"+" "+plot_title, 'interpreter', 'none');
g_mag=[coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx_loops; coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy_loops; coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz_loops].^2;
g_mag=sum(g_mag,1,'omitnan');
histogram(g_mag);
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx./coil_layouts.out.potential_step);
xlabel('[mT/m/A]');
hold off;
nexttile;
hold on
axis equal tight;
title('Gmag [mT/m/A]');
plot_colors=g_mag;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
%caxis([mean(g_mag).*0.9 mean(g_mag).*1.1]);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');

switch stongest_gradient
    case 1
nexttile;
hold on;
title("Gx"+" "+plot_title, 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx_loops(~isnan(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx_loops)));
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx./coil_layouts.out.potential_step);
xlabel('[mT/m/A]');
hold off;
nexttile;
hold on
axis equal tight;
title('Gx [mT/m/A]', 'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.local_gx_loops;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
caxis(plot_range);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
%Plot also the target gradient defined by the field shape
nexttile;
hold on;
title("Target Gx", 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx(~isnan(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx)));
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx./coil_layouts.out.potential_step);
xlabel('[mT/m]');
hold off;
nexttile;
hold on
axis equal tight;
title('Target Gx [mT/m]', 'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
caxis(plot_range);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');

    case 2
nexttile;
hold on;
title("Gy"+" "+plot_title, 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy_loops(~isnan(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy_loops)));
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gy./coil_layouts.out.potential_step);
xlabel('[mT/m/A]');
hold off;
hold off
nexttile;
hold on
axis equal tight;
title('Gy [mT/m/A]', 'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.local_gy_loops;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
caxis(plot_range);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off
%Plot also the target gradient defined by the field shape
nexttile;
hold on;
title("Target Gy", 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gy(~isnan(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gy)));
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx./coil_layouts.out.potential_step);
xlabel('[mT/m]');
hold off;
nexttile;
hold on
axis equal tight;
title('Target Gy [mT/m]', 'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gy;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
caxis(plot_range);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');

    case 3
nexttile;
hold on;
title("Gz"+" "+plot_title, 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz_loops(~isnan(coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz_loops)));
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gz./coil_layouts.out.potential_step);
xlabel('[mT/m/A]');
hold off;
nexttile;
hold on
axis equal tight;
title('Gz [mT/m/A]', 'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.local_gz_loops;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
caxis(plot_range);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
hold off
%Plot also the target gradient defined by the field shape
nexttile;
hold on;
title("Target Gz", 'interpreter', 'none');
histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gz(~isnan(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gz)));
%histogram(coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gx./coil_layouts.out.potential_step);
xlabel('[mT/m]');
hold off;
nexttile;
hold on
axis equal tight;
title('Target Gz [mT/m]', 'interpreter', 'none');
plot_colors=coil_layouts(single_ind_to_plot).out.layout_gradient.local_target_gz;
view(45,45);
colorbar;
colormap(parula);
scatter3(pos_data(1,:),pos_data(2,:),pos_data(3,:),dot_size*ones(1,numel(pos_data,2)),plot_colors,'filled');
caxis(plot_range);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
end

set(gcf,'color','w');

end

