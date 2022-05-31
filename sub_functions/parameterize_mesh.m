function coil_parts=parameterize_mesh(coil_parts,input)
%create the paramterized 2D mesh

%The non-cylidnrical parameterization is taken from "matlabmesh @ Ryan
%Schmidt  rms@dgp.toronto.edu" based on desbrun et al (2002), "Intrinsic Parameterizations of {Surface} Meshes",

surface_is_cylinder=input.surface_is_cylinder_flag;


for part_ind=1:numel(coil_parts)

face_normals=faceNormal(triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'));
max_face_normal_std=max([std(face_normals(:,1)) std(face_normals(:,2)) std(face_normals(:,3))]);

coil_parts(part_ind).coil_mesh.v=coil_parts(part_ind).coil_mesh.vertices';
coil_parts(part_ind).coil_mesh.fn=faceNormal(triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.v));

%Check if vertex coordinates are rather constant in one the three
%dimensions
if ~(max_face_normal_std<10^(-6)) %Check if vertex coordinates are rather constant in one the three dimensions


%go for the parameterization; distuingish between cylinder and non-cylinder

if ~surface_is_cylinder
coil_parts(part_ind).coil_mesh = mesh_parameterization_iterative( coil_parts(part_ind).coil_mesh );
%create a the 2D data set for fit
else
%Planarization of cylinder
boundary_edges=freeBoundary(triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'));
    %Build the boundary loops form the boundary edges
    is_new_node=[boundary_edges(:,1)' 0]==[0 boundary_edges(:,2)'];
    is_new_node(1)=1;
    is_new_node(end)=1;
    is_new_node=~is_new_node;
    num_boundaries=numel(find(is_new_node))+1;
    boundary_start=[1 find(is_new_node)];
    boundary_end=[find(is_new_node)-1 size(boundary_edges,1)];
    boundary_loop_nodes=cell(1,num_boundaries);
    for boundary_ind=1:num_boundaries
        boundary_loop_nodes{boundary_ind}=[boundary_edges(boundary_start(boundary_ind):boundary_end(boundary_ind),1); boundary_edges(boundary_start(boundary_ind),1)];
    end

%check if the cylinder is oriented along the z-axis
%if make a rotated copy for the parameterization
opening_mean=mean(coil_parts(part_ind).coil_mesh.vertices(:,boundary_loop_nodes{1}),2);
overall_mean=mean(coil_parts(part_ind).coil_mesh.vertices,2);
old_orientation_vector=(opening_mean-overall_mean)./repmat(vecnorm(opening_mean-overall_mean),[3 1]);
z_vec=[0 0 1]';
sina = vecnorm(cross(old_orientation_vector,z_vec)) / ( vecnorm(old_orientation_vector) * vecnorm(z_vec) );
cosa = vecnorm(dot(old_orientation_vector,z_vec)) / ( vecnorm(old_orientation_vector) * vecnorm(z_vec) );
angle = atan2( sina, cosa );
rotation_vector=cross(old_orientation_vector,[0 0 1]')./repmat(vecnorm(cross(old_orientation_vector,[0 0 1]')),[3 1]);
rot_mat= calc_3d_rotation_matrix_by_vector(rotation_vector,angle);
rotated_vectices=rot_mat*coil_parts(part_ind).coil_mesh.vertices;
point_coords=rotated_vectices;
min_z_cylinder=min(point_coords(3,:));
point_coords(3,:)=point_coords(3,:)+min_z_cylinder;
phi_coord=atan2(point_coords(2,:),point_coords(1,:));
r_coord=(point_coords(1,:).^2+point_coords(2,:).^2).^(1/2);
u_coord=(r_coord - point_coords(3,:)).*sin(phi_coord);
v_coord=(r_coord - point_coords(3,:)).*cos(phi_coord);
coil_parts(part_ind).coil_mesh.uv=[u_coord;v_coord];
coil_parts(part_ind).coil_mesh.n=vertexNormal(triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'))';
coil_parts(part_ind).coil_mesh.boundary=boundary_loop_nodes;

end

else  % the 3D mesh is already planar but the normals must be aligned to the z-axis
    
%Rotate the planar mesh in the xy plane    
mean_norm=mean(face_normals,1);
new_norm=[0 0 1];
v_c=cross(mean_norm,new_norm);
if norm(v_c)>10^(-8) %check wether the normals are already aligned to z
v_d=dot(mean_norm,new_norm);
mat_v=[0 (-1)*v_c(3) v_c(2); v_c(3) 0 (-1)*v_c(1); (-1)*v_c(2) v_c(1) 0];
rot_mat=eye(3)+mat_v+mat_v*mat_v*(1/(1+v_d));
out_a=sum(repmat(rot_mat(1,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1]).*coil_parts(part_ind).coil_mesh.vertices',2);
out_b=sum(repmat(rot_mat(2,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1]).*coil_parts(part_ind).coil_mesh.vertices',2);
out_c=sum(repmat(rot_mat(3,:),[size(coil_parts(part_ind).coil_mesh.vertices,2) 1]).*coil_parts(part_ind).coil_mesh.vertices',2);
coil_parts(part_ind).coil_mesh.uv=[out_a'; out_b'];
else
coil_parts(part_ind).coil_mesh.uv=[coil_parts(part_ind).coil_mesh.vertices(1,:); coil_parts(part_ind).coil_mesh.vertices(2,:)];
end

boundary_edges=freeBoundary(triangulation(coil_parts(part_ind).coil_mesh.faces',[coil_parts(part_ind).coil_mesh.uv; zeros(1,size(coil_parts(part_ind).coil_mesh.uv(1,:),2))]'));
%Build the boundary loops form the boundary edges
is_new_node=[boundary_edges(:,1)' 0]==[0 boundary_edges(:,2)'];
is_new_node(1)=1;
is_new_node(end)=1;
is_new_node=~is_new_node;
num_boundaries=numel(find(is_new_node))+1;
boundary_start=[1 find(is_new_node)];
boundary_end=[find(is_new_node)-1 size(boundary_edges,1)];
boundary_loop_nodes=cell(1,num_boundaries);
for boundary_ind=1:num_boundaries
boundary_loop_nodes{boundary_ind}=[boundary_edges(boundary_start(boundary_ind):boundary_end(boundary_ind),1); boundary_edges(boundary_start(boundary_ind),1)];
end
coil_parts(part_ind).coil_mesh.n=vertexNormal(triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'))';
coil_parts(part_ind).coil_mesh.boundary=boundary_loop_nodes;
    
    
end

end

 
end