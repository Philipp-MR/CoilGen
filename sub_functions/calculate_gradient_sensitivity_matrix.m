function coil_parts=calculate_gradient_sensitivity_matrix(coil_parts,target_field,input)


target_points=target_field.coords;
gauss_order=input.gauss_order;
[u_coord,v_coord,gauss_weight] = gauss_legendre_integration_points_triangle(gauss_order); % find uv-coords and weights for the gaus-legendre integration


% calculate the sensitivity matrix
coil_parts(numel(coil_parts)).gradient_sensitivity_matrix=[];

for part_ind=1:numel(coil_parts)

if ~input.temp_evalution.use_preoptimization_temp



% calculate the weights and the test point for the gauss legendre
% integration on each triangle is done by gauss legendre integration
biot_savart_coeff=10^(-7); % biot savart coefficient
plate_thickness=0.001; % in meter
num_nodes=numel(coil_parts(part_ind).basis_elements);
num_target_points=size(target_points,2);
gradient_sensitivity_matrix=zeros(3,num_target_points,num_nodes); %initialize the sensitivity matrix
num_triangles_per_node=arrayfun(@(x) numel(coil_parts(part_ind).basis_elements(x).area),1:numel(coil_parts(part_ind).basis_elements));
verts1=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(1,:));
verts2=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(2,:));
verts3=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).coil_mesh.faces(3,:));

x=target_points(1,:);
y=target_points(2,:);
z=target_points(3,:);


for node_ind=1:num_nodes % iterate of number of nodes
    
    
DBzdx=zeros(1,num_target_points);
DBzdy=zeros(1,num_target_points);
DBzdz=zeros(1,num_target_points);



for tri_ind=1:num_triangles_per_node(node_ind) % iterate of number of triangles of that node

    
    
% node_point=coil_parts(part_ind).coil_mesh.vertices(:,node_ind); 
% point_b=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).basis_elements(node_ind).one_ring(tri_ind,1));
% point_c=coil_parts(part_ind).coil_mesh.vertices(:,coil_parts(part_ind).basis_elements(node_ind).one_ring(tri_ind,2));
node_point=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,:,1);
point_b=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,:,2);
point_c=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,:,3);

    
%extract the coords of the triangle corners
x1=node_point(1);
y1=node_point(2);
z1=node_point(3);
x2=point_b(1);
y2=point_b(2);
z2=point_b(3);
x3=point_c(1);
y3=point_c(2);
z3=point_c(3);
    
    
% % %extract the coords of the triangle corners
% % x1=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,1,1); %1st ind: num triangle, 2nd ind: coords, 3th ind: corner ind
% % y1=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,2,1);
% % z1=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,3,1);
% % x2=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,1,2);
% % y2=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,2,2);
% % z2=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,3,2);
% % x3=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,1,3);
% % y3=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,2,3);
% % z3=coil_parts(part_ind).basis_elements(node_ind).triangle_points_ABC(tri_ind,3,3);

d_l_x=coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,1);
d_l_y=coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,2);
d_l_z=coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,3);

% current_vec=(point_c-point_b)./norm((point_c-point_b));
% dist_norm_factor=norm(cross(node_point-point_b,node_point-point_c))/norm(point_b-point_c);
% current_vec=current_vec./dist_norm_factor;
% vx=current_vec(1);
% vy=current_vec(2);
% vz=current_vec(3);

for gauss_ind=1:numel(gauss_weight) %build the sum over the gauss-legendre test points

l_x=x1*u_coord(gauss_ind)+x2*v_coord(gauss_ind)+x3*(1-u_coord(gauss_ind)-v_coord(gauss_ind)); %express the test points coords with uv coords
l_y=y1*u_coord(gauss_ind)+y2*v_coord(gauss_ind)+y3*(1-u_coord(gauss_ind)-v_coord(gauss_ind));
l_z=z1*u_coord(gauss_ind)+z2*v_coord(gauss_ind)+z3*(1-u_coord(gauss_ind)-v_coord(gauss_ind));

% l_x=(x1+x2+x3)./3;
% l_y=(y1+y2+y3)./3;
% l_z=(z1+z2+z3)./3;

E=((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(3/2);
dEdx=3.*(x-l_x).*((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(1/2);
dEdy=3.*(y-l_y).*((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(1/2);
dEdz=3.*(z-l_z).*((x-l_x).^(2)+(y-l_y).^(2)+(z-l_z).^(2)).^(1/2); 
C=10^(-7);
phi_x=(d_l_y.*(z-l_z)-d_l_z.*(y-l_y));
phi_y=(d_l_z.*(x-l_x)-d_l_x.*(z-l_z));
phi_z=(d_l_x.*(y-l_y)-d_l_y.*(x-l_x));
theta_factor=((E.*E).^(-1)).*(-1).*C;
%gradient of bz compnent
dBz_dx=theta_factor.*dEdx.*phi_z+(E.^(-1)).*C.*(-1).*d_l_y;
%gradient in y of bz compnent
dBz_dy=theta_factor.*dEdy.*phi_z+(E.^(-1)).*C.*d_l_x;
%gradient in z of bz compnent
dBz_dz=theta_factor.*dEdz.*phi_z;

DBzdx=DBzdx+dBz_dx.*gauss_weight(gauss_ind);
DBzdy=DBzdy+dBz_dy.*gauss_weight(gauss_ind);
DBzdz=DBzdz+dBz_dz.*gauss_weight(gauss_ind);


% %calcute the contribution to the sensitivity matrix and add them up
% dCx=dCx+((zgauss_in_uv-z_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,2)-(ygauss_in_uv-y_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,3)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
% dCy=dCy+((xgauss_in_uv-x_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,3)-(zgauss_in_uv-z_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,1)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
% dCz=dCz+((ygauss_in_uv-y_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,1)-(xgauss_in_uv-x_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,2)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
end
end
gradient_sensitivity_matrix(:,:,node_ind)=[DBzdx' DBzdy' DBzdz']'.*biot_savart_coeff.*plate_thickness;
end

coil_parts(part_ind).gradient_sensitivity_matrix=gradient_sensitivity_matrix;

else

coil_parts(part_ind).gradient_sensitivity_matrix=input.temp.coil_parts(part_ind).gradient_sensitivity_matrix;


end

end


end