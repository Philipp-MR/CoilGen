clc; clear all; close all;

%Select result parameters and plot them

if ispc
cd('..\');
else
cd('../');
end


load Results\results_parameter_study.mat





%% Plot results


plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{'all' [0.2] 'x' [1 0 0 0] 'all'},'Plate Size [m]',5);
plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{'all' 'all' 'x' [1 0 0 0] [0.5]},'Plate Distance [m]',2);




function plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,target_parameters,x_label,para_plot_ind)

inds(numel(target_parameters)).match=[];
plot_parameter.max_error=[];
plot_parameter.mean_error=[];
plot_parameter.sensitivity=[];
plot_parameter.wire_length=[];
plot_parameter.tikonov_value=[];
plot_parameter.plate_distances=[];
plot_parameter.plate_size=[];


if ~strcmp(target_parameters{1},'all')
inds(1).match=find(para_ind_grid{1}==find(sampling_parameters.tikonov_values==target_parameters{1}));
else
inds(1).match=find(para_ind_grid{1});
end
if ~strcmp(target_parameters{2},'all')
inds(2).match=find(para_ind_grid{2}==find(sampling_parameters.plate_distances==target_parameters{2}));
else
inds(2).match=find(para_ind_grid{2});
end
if ~strcmp(target_parameters{3},'all')
inds(3).match=find(para_ind_grid{3}==find(strcmp(sampling_parameters.field_shape_functions,target_parameters{3})));
else
inds(3).match=find(para_ind_grid{3});
end
if ~strcmp(target_parameters{4},'all')
inds(4).match=find(para_ind_grid{4}==find(arrayfun(@(x) isequal(sampling_parameters.plate_alignements{x},target_parameters{4}),1:numel(sampling_parameters.plate_alignements))));
else
inds(4).match=find(para_ind_grid{4});
end
if ~strcmp(target_parameters{5},'all')
inds(5).match=find(para_ind_grid{5}==find(sampling_parameters.plate_size==target_parameters{5}));
else
inds(5).match=find(para_ind_grid{5});
end


%find the case indices which match all the criteria
inds_out=inds(1).match;
if numel(inds)>1
for array_ind=2:numel(inds)
        inds_out= intersect(inds_out,inds(array_ind).match);
end
end

%write the paramters for the matching cases
if ~isempty(inds_out)
plot_parameter.max_error=result_parameters.max_error(inds_out);
plot_parameter.mean_error=result_parameters.mean_error(inds_out);
plot_parameter.sensitivity=result_parameters.sensitivity(inds_out);
plot_parameter.wire_length=result_parameters.wire_length(inds_out);
plot_parameter.tikonov_value=sampling_parameters.tikonov_values(para_ind_grid{1}(inds_out));
plot_parameter.plate_distances=sampling_parameters.plate_distances(para_ind_grid{2}(inds_out));
plot_parameter.plate_size=sampling_parameters.plate_size(para_ind_grid{5}(inds_out));
end

switch para_plot_ind
    case 1
x_plot_vals=plot_parameter.tikonov_values;
    case 2
x_plot_vals=plot_parameter.plate_distances;
    case 3
x_plot_vals=plot_parameter.field_shape_functions;
    case 4
x_plot_vals=plot_parameter.plate_alignements;
    case 5
x_plot_vals=plot_parameter.plate_size;
end

figure;
hold on;
plot(x_plot_vals,plot_parameter.sensitivity,'LineWidth',2);
grid on;
%axis equal;
xlabel(x_label, 'Interpreter', 'none');
ylabel('Sensitivity [mT\m\A]', 'Interpreter', 'none');
axis([min(x_plot_vals) max(x_plot_vals) 0 max(plot_parameter.sensitivity)]);
set(gcf,'color','w');
set(gcf, 'Position',  [100, 100, 400, 400])
hold off;


end




% % % best_orientation_x=plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{7.943282347242805 0.2 'x' 0.4 'all'},'Alignment [m]',4);
% % % best_orientation_y=plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{7.943282347242805 0.2 'y' 0.4 'all'},'Alignment [m]',4);
% % % best_orientation_z=plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,{7.943282347242805 0.2 'z' 0.4 'all'},'Alignment [m]',4);
% % % 
% % % 
% % % 
% % % 
% % % 
% % % function best_orientation=plot_parameter_set(para_ind_grid,sampling_parameters,result_parameters,target_parameters,x_label,para_plot_ind)
% % % 
% % % inds(numel(target_parameters)).match=[];
% % % plot_parameter.max_error=[];
% % % plot_parameter.mean_error=[];
% % % plot_parameter.sensitivity=[];
% % % plot_parameter.wire_length=[];
% % % plot_parameter.tikonov_value=[];
% % % plot_parameter.plate_distances=[];
% % % plot_parameter.plate_size=[];
% % % 
% % % 
% % % if ~strcmp(target_parameters{1},'all')
% % % inds(1).match=find(para_ind_grid{1}==find(sampling_parameters.tikonov_values==target_parameters{1}));
% % % else
% % % inds(1).match=find(para_ind_grid{1});
% % % end
% % % if ~strcmp(target_parameters{2},'all')
% % % inds(2).match=find(para_ind_grid{2}==find(sampling_parameters.plate_distances==target_parameters{2}));
% % % else
% % % inds(2).match=find(para_ind_grid{2});
% % % end
% % % if ~strcmp(target_parameters{3},'all')
% % % inds(3).match=find(para_ind_grid{3}==find(strcmp(sampling_parameters.field_shape_functions,target_parameters{3})));
% % % else
% % % inds(3).match=find(para_ind_grid{3},'all');
% % % end
% % % if ~strcmp(target_parameters{4},'all')
% % % inds(4).match=find(para_ind_grid{4}==find(sampling_parameters.plate_size==target_parameters{4}));
% % % else
% % % inds(4).match=find(para_ind_grid{4});
% % % end
% % % if ~strcmp(target_parameters{5},'all')
% % % inds(5).match=find(para_ind_grid{5}==find(arrayfun(@(x) isequal(sampling_parameters.plate_alignements{x},target_parameters{5}),1:numel(sampling_parameters.plate_alignements))));
% % % else
% % % inds(5).match=find(para_ind_grid{5});
% % % end
% % % 
% % % %find the case indices which match all the criteria
% % % inds_out=inds(1).match;
% % % if numel(inds)>1
% % % for array_ind=2:numel(inds)
% % %         inds_out= intersect(inds_out,inds(array_ind).match);
% % % end
% % % end
% % % 
% % % %write the paramters for the matching cases
% % % if ~isempty(inds_out)
% % % plot_parameter.max_error=result_parameters.max_error(inds_out);
% % % plot_parameter.mean_error=result_parameters.mean_error(inds_out);
% % % plot_parameter.sensitivity=result_parameters.sensitivity(inds_out);
% % % plot_parameter.wire_length=result_parameters.wire_length(inds_out);
% % % plot_parameter.tikonov_value=sampling_parameters.tikonov_values(para_ind_grid{1}(inds_out));
% % % plot_parameter.plate_distances=sampling_parameters.plate_distances(para_ind_grid{2}(inds_out));
% % % plot_parameter.plate_size=sampling_parameters.plate_size(para_ind_grid{4}(inds_out));
% % % plot_parameter.plate_alignements=sampling_parameters.plate_alignements(para_ind_grid{5}(inds_out));
% % % end
% % % 
% % % switch para_plot_ind
% % %     case 1
% % % x_plot_vals=plot_parameter.tikonov_values;
% % %     case 2
% % % x_plot_vals=plot_parameter.plate_distances;
% % %     case 3
% % % x_plot_vals=plot_parameter.field_shape_functions;
% % %     case 4
% % % x_plot_vals=plot_parameter.plate_size;
% % %     case 5
% % % x_plot_vals=plot_parameter.plate_alignements;
% % % end
% % % 
% % % X=arrayfun(@(x) plot_parameter.plate_alignements{x}(1),1:numel(plot_parameter.plate_alignements));
% % % Y=arrayfun(@(x) plot_parameter.plate_alignements{x}(2),1:numel(plot_parameter.plate_alignements));
% % % Z=arrayfun(@(x) plot_parameter.plate_alignements{x}(3),1:numel(plot_parameter.plate_alignements));
% % % thetas=acos(Z./((X.^2+Y.^2+Z.^2).^(1/2))).*(360/2/pi);
% % % phis=atan2(Y,X).*(360/2/pi);
% % % delaunay_tris = delaunay(thetas,phis);
% % % delaunay_tris_3d= [convhull(X,Y,Z)]';
% % % 
% % % [~,max_ind_sensitivity]=max(plot_parameter.sensitivity);
% % % best_orientation=plot_parameter.plate_alignements{max_ind_sensitivity};
% % % 
% % % 
% % % 
% % % figure;
% % % hold on;
% % % title("Alignment vs Sensitivity, Channel "+" "+ target_parameters{3});
% % % normed_pot=plot_parameter.sensitivity;
% % % normed_pot=normed_pot./max(abs(normed_pot))*100;
% % % x_vals_1=X(delaunay_tris_3d(1,:));
% % % x_vals_2=X(delaunay_tris_3d(2,:));
% % % x_vals_3=X(delaunay_tris_3d(3,:));
% % % y_vals_1=Y(delaunay_tris_3d(1,:));
% % % y_vals_2=Y(delaunay_tris_3d(2,:));
% % % y_vals_3=Y(delaunay_tris_3d(3,:));
% % % z_vals_1=Z(delaunay_tris_3d(1,:));
% % % z_vals_2=Z(delaunay_tris_3d(2,:));
% % % z_vals_3=Z(delaunay_tris_3d(3,:));
% % % normed_pot_1=normed_pot(delaunay_tris_3d(1,:));
% % % normed_pot_2=normed_pot(delaunay_tris_3d(2,:));
% % % normed_pot_3=normed_pot(delaunay_tris_3d(3,:));
% % % fill3([x_vals_1' x_vals_2' x_vals_3']',[y_vals_1' y_vals_2' y_vals_3']',[z_vals_1' z_vals_2' z_vals_3']',[normed_pot_1 normed_pot_2 normed_pot_3]','EdgeAlpha',0);
% % % %scatter3(X,Y,Z,1/2.*ones(size(X)),"black",'filled');
% % % scatter3(best_orientation(1),best_orientation(2),best_orientation(3),80,'red','filled');
% % % %Plot also a spherical coordiante system
% % % azimuth_resol=10;
% % % polar_resol=10;
% % % z_val=sin((0:(1/(azimuth_resol-1)):1).*(pi/2));
% % % cut_radia=(1-(z_val).^2).^(1/2);
% % % lin_resol=100;
% % % lin_points_num=floor(cut_radia.*lin_resol);
% % % for az_ind=1:azimuth_resol
% % % lin_incr=0:1/(lin_points_num(az_ind)-1):1;
% % % x_vals=sin(lin_incr.*(2*pi)).*cut_radia(az_ind).*1.001;
% % % y_vals=cos(lin_incr.*(2*pi)).*cut_radia(az_ind).*1.001;
% % % z_vals=z_val(az_ind).*ones(size(lin_incr)).*1.001;
% % % plot3(x_vals,y_vals,z_vals,'-k','LineWidth',1);
% % % end
% % % lin_incr=0:1/(lin_resol-1):1;
% % % x_vals=sin(lin_incr.*pi./2);
% % % y_vals=zeros(size(x_vals));
% % % z_vals=cos(lin_incr.*pi./2);
% % % for pol_ind=0:polar_resol
% % % rot_matrix=[cos(pol_ind/polar_resol*2*pi) (-1).*sin(pol_ind/polar_resol*2*pi); sin(pol_ind/polar_resol*2*pi) cos(pol_ind/polar_resol*2*pi)];
% % % rotated_mat=rot_matrix*[x_vals;y_vals];
% % % plot3(rotated_mat(1,:).*1.001,rotated_mat(2,:).*1.001,z_vals.*1.001,'k','LineWidth',1);
% % % end
% % % xlabel('X');
% % % ylabel('Y');
% % % zlabel('Z');
% % % axis equal;
% % % grid on;
% % % c = colorbar;
% % % c.Label.String = 'Sensitivity (Relativ to Best) [%]';
% % % view(45,45);
% % % set(gcf,'color','w');
% % % hold off;
% % % 
% % % 
% % % figure;
% % % hold on;
% % % title("Alignment vs Field Error, Channel "+" "+ target_parameters{3});
% % % normed_pot=plot_parameter.mean_error;
% % % normed_pot=normed_pot./max(abs(normed_pot))*100;
% % % x_vals_1=X(delaunay_tris_3d(1,:));
% % % x_vals_2=X(delaunay_tris_3d(2,:));
% % % x_vals_3=X(delaunay_tris_3d(3,:));
% % % y_vals_1=Y(delaunay_tris_3d(1,:));
% % % y_vals_2=Y(delaunay_tris_3d(2,:));
% % % y_vals_3=Y(delaunay_tris_3d(3,:));
% % % z_vals_1=Z(delaunay_tris_3d(1,:));
% % % z_vals_2=Z(delaunay_tris_3d(2,:));
% % % z_vals_3=Z(delaunay_tris_3d(3,:));
% % % normed_pot_1=normed_pot(delaunay_tris_3d(1,:));
% % % normed_pot_2=normed_pot(delaunay_tris_3d(2,:));
% % % normed_pot_3=normed_pot(delaunay_tris_3d(3,:));
% % % fill3([x_vals_1' x_vals_2' x_vals_3']',[y_vals_1' y_vals_2' y_vals_3']',[z_vals_1' z_vals_2' z_vals_3']',[normed_pot_1 normed_pot_2 normed_pot_3]','EdgeAlpha',0);
% % % %scatter3(X,Y,Z,1/2.*ones(size(X)),"black",'filled');
% % % scatter3(best_orientation(1),best_orientation(2),best_orientation(3),80,'red','filled');
% % % %Plot also a spherical coordiante system
% % % azimuth_resol=10;
% % % polar_resol=10;
% % % z_val=sin((0:(1/(azimuth_resol-1)):1).*(pi/2));
% % % cut_radia=(1-(z_val).^2).^(1/2);
% % % lin_resol=100;
% % % lin_points_num=floor(cut_radia.*lin_resol);
% % % for az_ind=1:azimuth_resol
% % % lin_incr=0:1/(lin_points_num(az_ind)-1):1;
% % % x_vals=sin(lin_incr.*(2*pi)).*cut_radia(az_ind).*1.001;
% % % y_vals=cos(lin_incr.*(2*pi)).*cut_radia(az_ind).*1.001;
% % % z_vals=z_val(az_ind).*ones(size(lin_incr)).*1.001;
% % % plot3(x_vals,y_vals,z_vals,'-k','LineWidth',1);
% % % end
% % % lin_incr=0:1/(lin_resol-1):1;
% % % x_vals=sin(lin_incr.*pi./2);
% % % y_vals=zeros(size(x_vals));
% % % z_vals=cos(lin_incr.*pi./2);
% % % for pol_ind=0:polar_resol
% % % rot_matrix=[cos(pol_ind/polar_resol*2*pi) (-1).*sin(pol_ind/polar_resol*2*pi); sin(pol_ind/polar_resol*2*pi) cos(pol_ind/polar_resol*2*pi)];
% % % rotated_mat=rot_matrix*[x_vals;y_vals];
% % % plot3(rotated_mat(1,:).*1.001,rotated_mat(2,:).*1.001,z_vals.*1.001,'k','LineWidth',1);
% % % end
% % % xlabel('X');
% % % ylabel('Y');
% % % zlabel('Z');
% % % axis equal;
% % % grid on;
% % % c = colorbar;
% % % c.Label.String = 'Sensitivity (Relativ to Best) [%]';
% % % view(45,45);
% % % set(gcf,'color','w');
% % % hold off;
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % figure;
% % % hold on;
% % % title("Alignment vs Sensitivity, Channel "+" "+ target_parameters{3});
% % % patch('Faces',delaunay_tris,'Vertices',[phis; thetas]','FaceVertexCData',plot_parameter.sensitivity,'FaceColor','interp','edgealpha',0.0);
% % % scatter(phis,thetas,10*ones(size(thetas)),'black','filled');
% % % ylim([0 90]);
% % % xlim([-181 181]);
% % % xlabel('Phi [Degree]');
% % % ylabel('Theta [Degree]');
% % % axis equal;
% % % grid on;
% % % c = colorbar;
% % % c.Label.String = 'Sensitivity [mT/m/A]';
% % % hold off;
% % % 
% % % 
% % % 
% % % % figure;
% % % % hold on;
% % % % plot(x_plot_vals,plot_parameter.sensitivity,'LineWidth',2);
% % % % grid on;
% % % % %axis equal;
% % % % xlabel(x_label, 'Interpreter', 'none');
% % % % ylabel('Sensitivity [mT\m\A]', 'Interpreter', 'none');
% % % % axis([min(x_plot_vals) max(x_plot_vals) 0 max(plot_parameter.sensitivity)]);
% % % % set(gcf,'color','w');
% % % % set(gcf, 'Position',  [100, 100, 400, 400])
% % % % hold off;
% % % 
% % % 
% % % end




