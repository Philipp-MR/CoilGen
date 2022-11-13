function coil_parts=calculate_sensitivity_matrix(coil_parts,target_field,input)
% calculate the sensitivity matrix



coil_parts(numel(coil_parts)).sensitivity_matrix=[];

for part_ind=1:numel(coil_parts)


if ~input.temp_evalution.use_preoptimization_temp

target_points=target_field.coords;
gauss_order=input.gauss_order;

% calculate the weights and the test point for the gauss legendre
% integration on each triangle is done by gauss legendre integration
[u_coord,v_coord,gauss_weight] = gauss_legendre_integration_points_triangle(gauss_order); % find uv-coords and weights for the gaus-legendre integration
num_gauss_points=numel(gauss_weight); % the number of gaus points
biot_savart_coeff=10^(-7); % biot savart coefficient
num_nodes=numel(coil_parts(part_ind).basis_elements);
num_target_points=size(target_points,2);
sensitivity_matrix=zeros(3,num_target_points,num_nodes); %initialize the sensitivity matrix
num_triangles_per_node=arrayfun(@(x) numel(coil_parts(part_ind).basis_elements(x).area),1:numel(coil_parts(part_ind).basis_elements));

x_target=target_points(1,:);
y_target=target_points(2,:);
z_target=target_points(3,:);


for node_ind=1:num_nodes % iterate of number of nodes
    
    
dCx=zeros(1,num_target_points);
dCy=zeros(1,num_target_points);
dCz=zeros(1,num_target_points);
    
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

vx=coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,1);
vy=coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,2);
vz=coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,3);

% current_vec=(point_c-point_b)./norm((point_c-point_b));
% dist_norm_factor=norm(cross(node_point-point_b,node_point-point_c))/norm(point_b-point_c);
% current_vec=current_vec./dist_norm_factor;
% vx=current_vec(1);
% vy=current_vec(2);
% vz=current_vec(3);

for gauss_ind=1:num_gauss_points %build the sum over the gauss-legendre test points
xgauss_in_uv=x1*u_coord(gauss_ind)+x2*v_coord(gauss_ind)+x3*(1-u_coord(gauss_ind)-v_coord(gauss_ind)); %express the test points coords with uv coords
ygauss_in_uv=y1*u_coord(gauss_ind)+y2*v_coord(gauss_ind)+y3*(1-u_coord(gauss_ind)-v_coord(gauss_ind));
zgauss_in_uv=z1*u_coord(gauss_ind)+z2*v_coord(gauss_ind)+z3*(1-u_coord(gauss_ind)-v_coord(gauss_ind));


% calculate the distance norm for the biot savart law
distance_norm=((xgauss_in_uv-x_target).^2+(ygauss_in_uv-y_target).^2+(zgauss_in_uv-z_target).^2).^(-3/2);


%calcute the contribution to the sensitivity matrix and add them up
dCx=dCx+((-1).*vz.*(y_target-ygauss_in_uv)+vy.*(z_target-zgauss_in_uv)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
dCy=dCy+((-1).*vx.*(z_target-zgauss_in_uv)+vz.*(x_target-xgauss_in_uv)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
dCz=dCz+((-1).*vy.*(x_target-xgauss_in_uv)+vx.*(y_target-ygauss_in_uv)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);

% %calcute the contribution to the sensitivity matrix and add them up
% dCx=dCx+((zgauss_in_uv-z_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,2)-(ygauss_in_uv-y_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,3)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
% dCy=dCy+((xgauss_in_uv-x_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,3)-(zgauss_in_uv-z_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,1)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
% dCz=dCz+((ygauss_in_uv-y_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,1)-(xgauss_in_uv-x_target).*coil_parts(part_ind).basis_elements(node_ind).current(tri_ind,2)).*distance_norm.*2.*coil_parts(part_ind).basis_elements(node_ind).area(tri_ind).*gauss_weight(gauss_ind);
end
end
sensitivity_matrix(:,:,node_ind)=[dCx' dCy' dCz']'.*biot_savart_coeff;
end

coil_parts(part_ind).sensitivity_matrix=sensitivity_matrix;

else

coil_parts(part_ind).sensitivity_matrix=input.temp.coil_parts(part_ind).sensitivity_matrix;

end
end
end