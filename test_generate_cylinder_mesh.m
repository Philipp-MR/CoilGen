clc; clear all; 
close all;

cylinder_height=1; % in m

cylinder_radius=0.3;

num_arch_segments=20;
num_cylinder_rings=20;

%build circular cut shapes
x_positions=[sin(0:(2*pi)/(num_arch_segments):2*pi)].*cylinder_radius;
y_positions=[cos(0:(2*pi)/(num_arch_segments):2*pi)].*cylinder_radius;
x_positions(end)=[]; y_positions(end)=[]; %remove repetions at the ends
z_positions=((-1)*cylinder_height/2):(cylinder_height/(num_cylinder_rings+1)):(cylinder_height/2);

vertices_x=repmat(x_positions,[1 size(z_positions,2)-1]);
vertices_y=repmat(y_positions,[1 size(z_positions,2)-1]);
vertices_z=repelem(z_positions(1:end-1),size(x_positions,2));

vertices=[vertices_x; vertices_y; vertices_z];

tri_1_vert_inds_1=1:((num_arch_segments)*(num_cylinder_rings));
tri_1_vert_inds_2=tri_1_vert_inds_1+1;
%take care of indice overflow at the end of the rings
tri_1_vert_inds_2(mod(tri_1_vert_inds_2-1,num_arch_segments)==0)=tri_1_vert_inds_2(mod(tri_1_vert_inds_2-1,num_arch_segments)==0)-num_arch_segments;
tri_1_vert_inds_3=tri_1_vert_inds_2+num_arch_segments;

tri_2_vert_inds_1=tri_1_vert_inds_1;
tri_2_vert_inds_2=tri_1_vert_inds_3;
tri_2_vert_inds_3=[1:((num_arch_segments)*(num_cylinder_rings))]+num_arch_segments;

faces_1=[tri_1_vert_inds_2; tri_1_vert_inds_1; tri_1_vert_inds_3];
faces_2=[tri_2_vert_inds_2; tri_2_vert_inds_1; tri_2_vert_inds_3];

tri_mesh=triangulation([faces_1';faces_2'],vertices');

stlwrite(tri_mesh,cd+"\test.stl",'text') 

figure; 
hold on;
axis equal;
trisurf(tri_mesh);
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
view(45,45)
hold off;
