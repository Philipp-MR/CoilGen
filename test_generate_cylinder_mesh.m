clc; clear all; close all;

cylinder_length=1; % in m

cylinder_radius=0.3;

num_arch_segments=20;

num_cylinder_rings=10;

%build circular cut shapes
opening_circle=[sin(0:(2*pi)/(num_arch_segments-1):2*pi); cos(0:(2*pi)/(num_arch_segments-1):2*pi)];
opening_circle(:,end)=[];
z_positions=((-1)*cylinder_length/2):(cylinder_length)/(num_cylinder_rings-1):(cylinder_length/2);

vertices_x=repmat(opening_circle(1,:),[1 size(z_positions,2)]);
vertices_y=repmat(opening_circle(2,:),[1 size(z_positions,2)]);
vertices_z=repelem(z_positions,size(opening_circle,2));

vertices=[vertices_x; vertices_y; vertices_z];

num_vertices=size(vertices,2);

circle_inds=repmat(1:size(vertices,2),[1 size(z_positions,2)]);

tri_1_vert_inds_1=1:(num_vertices-(num_arch_segments));
tri_1_vert_inds_2=circshift(1:(num_vertices-(num_arch_segments)),-1);
%take care of indice overflow at the end of the rings
%tri_1_vert_inds_2(mod(tri_1_vert_inds_2,num_arch_segments+1)==0)=tri_1_vert_inds_2(mod(tri_1_vert_inds_2,num_arch_segments+1)==0)-(num_arch_segments-1);
%tri_1_vert_inds_3=circshift(1:(num_vertices-numel(z_positions)),(-1)*numel(z_positions));
tri_1_vert_inds_3=(num_arch_segments):(num_vertices-1);


% tri_2_vert_inds_1=
% tri_2_vert_inds_2=
% tri_2_vert_inds_3=

faces=[tri_1_vert_inds_1; tri_1_vert_inds_2; tri_1_vert_inds_3];


tri_mesh=triangulation(faces',vertices');

figure; 
hold on;
axis equal;
trisurf(tri_mesh);
hold off;
