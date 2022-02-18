function inductance_matrix=calculate_inductance_matrix(basis_elements,gauss_order,is_real_triangle_mat,triangle_corner_coord_mat,current_mat,area_mat)
% calculate the inductance matrix Lmn


% % % % % calculate the weights and the test point for the gauss legendre
% % % % % integration on each triangle is done by gauss legendre integration
% % % % [u_coord,v_coord,gauss_weight] = gauss_legendre_integration_points_triangle(gauss_order); % find uv-coords and weights for the gaus-legendre integration
% % % % num_gauss_points=numel(gauss_weight); % the number of gaus points
% % % % biot_savart_coeff=10^(-7); % biot savart coefficient
% % % % num_nodes=numel(basis_elements);
% % % % inductance_matrix=zeros(num_nodes,num_nodes); %initialize the sensitivity matrix
% % % % 
% % % % triangle_num=size(current_mat,2);
% % % % 
% % % % %extract the coords of the triangle corners
% % % % x11=triangle_corner_coord_mat(:,:,1,1);
% % % % y11=triangle_corner_coord_mat(:,:,2,1);
% % % % z11=triangle_corner_coord_mat(:,:,3,1);
% % % % x12=triangle_corner_coord_mat(:,:,1,2);
% % % % y12=triangle_corner_coord_mat(:,:,2,2);
% % % % z12=triangle_corner_coord_mat(:,:,3,2);
% % % % x13=triangle_corner_coord_mat(:,:,1,3);
% % % % y13=triangle_corner_coord_mat(:,:,2,3);
% % % % z13=triangle_corner_coord_mat(:,:,3,3);
% % % % 
% % % % x21=triangle_corner_coord_mat(:,:,1,1);
% % % % y21=triangle_corner_coord_mat(:,:,2,1);
% % % % z21=triangle_corner_coord_mat(:,:,3,1);
% % % % x22=triangle_corner_coord_mat(:,:,1,2);
% % % % y22=triangle_corner_coord_mat(:,:,2,2);
% % % % z22=triangle_corner_coord_mat(:,:,3,2);
% % % % x23=triangle_corner_coord_mat(:,:,1,3);
% % % % y23=triangle_corner_coord_mat(:,:,2,3);
% % % % z23=triangle_corner_coord_mat(:,:,3,3);
% % % % 
% % % % 
% % % % % for aaaa=1:num_nodes
% % % % % for bbbb=1:num_nodes 
% % % % % for cccc=1:num_triangles_per_node(aaaa)
% % % % % for dddd=1:num_triangles_per_node(bbbb)
% % % % 
% % % % dot_product=zeros(num_nodes,num_nodes,triangle_num,triangle_num);
% % % % for hhhh=1:triangle_num
% % % % for gggg=1:triangle_num
% % % % dot_product(:,:,hhhh,gggg)=squeeze(current_mat(:,hhhh,:))*squeeze(current_mat(:,gggg,:))';
% % % % end
% % % % end
% % % % 
% % % % for hhhh=1:triangle_num
% % % % for gggg=1:triangle_num
% % % % for eeee=1:num_gauss_points
% % % % for ffff=1:num_gauss_points
% % % % x1gauss_in_uv=x11(:,hhhh)*u_coord(eeee)+x12(:,hhhh)*v_coord(eeee)+x13(:,hhhh)*(1-u_coord(eeee)-v_coord(eeee)); %express the test points coords with uv coords
% % % % y1gauss_in_uv=y11(:,hhhh)*u_coord(eeee)+y12(:,hhhh)*v_coord(eeee)+y13(:,hhhh)*(1-u_coord(eeee)-v_coord(eeee));
% % % % z1gauss_in_uv=z11(:,hhhh)*u_coord(eeee)+z12(:,hhhh)*v_coord(eeee)+z13(:,hhhh)*(1-u_coord(eeee)-v_coord(eeee));
% % % % x2gauss_in_uv=x21(:,gggg)*u_coord(ffff)+x22(:,gggg)*v_coord(ffff)+x23(:,gggg)*(1-u_coord(ffff)-v_coord(ffff)); %express the test points coords with uv coords
% % % % y2gauss_in_uv=y21(:,gggg)*u_coord(ffff)+y22(:,gggg)*v_coord(ffff)+y23(:,gggg)*(1-u_coord(ffff)-v_coord(ffff));
% % % % z2gauss_in_uv=z21(:,gggg)*u_coord(ffff)+z22(:,gggg)*v_coord(ffff)+z23(:,gggg)*(1-u_coord(ffff)-v_coord(ffff));
% % % % % calculate the distance norm for the biot savart law
% % % % distance_norm=((x1gauss_in_uv-x2gauss_in_uv').^2+(y1gauss_in_uv-y2gauss_in_uv').^2+(z1gauss_in_uv-z2gauss_in_uv').^2).^(-3/2);
% % % % 
% % % % distance_norm(isinf(distance_norm)|isnan(distance_norm)) = 0;
% % % % 
% % % % fact_a=dot_product(:,:,hhhh,gggg);
% % % % fact_b=2*[area_mat(:,hhhh)*area_mat(:,gggg)'];
% % % % fact_c=gauss_weight(eeee)*gauss_weight(ffff);
% % % % fact_d=[is_real_triangle_mat(:,hhhh)*is_real_triangle_mat(:,gggg)'];
% % % % 
% % % % 
% % % % % isfinite(distance_norm)
% % % % % 
% % % % % distance_norm(isfinite(distance_norm))
% % % % 
% % % % inductance_matrix=inductance_matrix+fact_a.*distance_norm.*fact_b.*fact_c.*fact_d;
% % % % 
% % % % end
% % % % end
% % % % end
% % % % end
% % % % % end
% % % % % end
% % % % % end
% % % % % end
% % % % 
% % % % 
% % % % %take care of the singularities of triangle self inductance
% % % % for aaaa=1:num_nodes
% % % % triangual_inductance_self_contribution=0;
% % % % for bbbb=num_triangles_per_node(aaaa)
% % % % for cccc=num_triangles_per_node(aaaa)
% % % % r_1=squeeze(basis_elements(aaaa).triangle_points_ABC(bbbb,:,1));
% % % % r_2=squeeze(basis_elements(aaaa).triangle_points_ABC(bbbb,:,2));
% % % % r_3=squeeze(basis_elements(aaaa).triangle_points_ABC(bbbb,:,3));
% % % % term_a=dot((r_3-r_1),(r_3-r_1));
% % % % term_b=dot((r_3-r_1),(r_3-r_2));
% % % % term_c=dot((r_3-r_2),(r_3-r_2));
% % % % dot_product=basis_elements(aaaa).current(bbbb,:)*basis_elements(aaaa).current(cccc,:)';
% % % % pre_factor=4*basis_elements(aaaa).area(bbbb)*basis_elements(aaaa).area(cccc);
% % % % 
% % % % triangual_inductance_self_contribution=triangual_inductance_self_contribution+pre_factor*dot_product*( ...
% % % % 1/(6*sqrt(term_a))*log( ( term_a-term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c)*(term_b+sqrt(term_a*term_c)) ) / ...
% % % %                                     ( ((-1)*term_b+sqrt(term_a*term_c))*((-1)*term_a+term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c)))) ...
% % % % +1/(6*sqrt(term_c))*log( ( (term_b+sqrt(term_a*term_c))* ((-1)*term_b+term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))) /...
% % % %                                         (( term_b-term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))*((-1)*term_b+sqrt(term_a*term_c)))) ...
% % % % +1/(6*sqrt(term_a-2*term_b+term_c))*log( ( (term_a-term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c))*((-1)*term_b+term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))) /...
% % % %                                                                         ( (term_b-term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))*((-1)*term_a+term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c)))));
% % % % end
% % % % end
% % % % inductance_matrix(aaaa,aaaa)=triangual_inductance_self_contribution;
% % % % end
% % % % 
% % % % inductance_matrix=inductance_matrix.*biot_savart_coeff;

% calculate the weights and the test point for the gauss legendre
% integration on each triangle is done by gauss legendre integration
[u_coord,v_coord,gauss_weight] = gauss_legendre_integration_points_triangle(gauss_order); % find uv-coords and weights for the gaus-legendre integration
num_gauss_points=numel(gauss_weight); % the number of gaus points
biot_savart_coeff=10^(-7); % biot savart coefficient
num_nodes=numel(basis_elements);
inductance_matrix=zeros(num_nodes,num_nodes); %initialize the sensitivity matrix
num_triangles_per_node=[arrayfun(@(x) numel(basis_elements(x).area),1:num_nodes)]';

for node_ind1=1:num_nodes
for node_ind2=1:num_nodes


for aaaa=1:num_nodes
for bbbb=1:num_nodes 
for cccc=1:num_triangles_per_node(aaaa)
for dddd=1:num_triangles_per_node(bbbb)
    

%extract the coords of the triangle corners
x11=basis_elements(aaaa).triangle_points_ABC(cccc,1,1);
y11=basis_elements(aaaa).triangle_points_ABC(cccc,2,1);
z11=basis_elements(aaaa).triangle_points_ABC(cccc,3,1);
x12=basis_elements(aaaa).triangle_points_ABC(cccc,1,2);
y12=basis_elements(aaaa).triangle_points_ABC(cccc,2,2);
z12=basis_elements(aaaa).triangle_points_ABC(cccc,3,2);
x13=basis_elements(aaaa).triangle_points_ABC(cccc,1,3);
y13=basis_elements(aaaa).triangle_points_ABC(cccc,2,3);
z13=basis_elements(aaaa).triangle_points_ABC(cccc,3,3);

x21=basis_elements(bbbb).triangle_points_ABC(dddd,1,1);
y21=basis_elements(bbbb).triangle_points_ABC(dddd,2,1);
z21=basis_elements(bbbb).triangle_points_ABC(dddd,3,1);
x22=basis_elements(bbbb).triangle_points_ABC(dddd,1,2);
y22=basis_elements(bbbb).triangle_points_ABC(dddd,2,2);
z22=basis_elements(bbbb).triangle_points_ABC(dddd,3,2);
x23=basis_elements(bbbb).triangle_points_ABC(dddd,1,3);
y23=basis_elements(bbbb).triangle_points_ABC(dddd,2,3);
z23=basis_elements(bbbb).triangle_points_ABC(dddd,3,3);

dot_product=basis_elements(aaaa).current(cccc,:)*basis_elements(bbbb).current(dddd,:)';

for eeee=1:num_gauss_points
for ffff=1:num_gauss_points
x1gauss_in_uv=x11*u_coord(eeee)+x12*v_coord(eeee)+x13*(1-u_coord(eeee)-v_coord(eeee)); %express the test points coords with uv coords
y1gauss_in_uv=y11*u_coord(eeee)+y12*v_coord(eeee)+y13*(1-u_coord(eeee)-v_coord(eeee));
z1gauss_in_uv=z11*u_coord(eeee)+z12*v_coord(eeee)+z13*(1-u_coord(eeee)-v_coord(eeee));
x2gauss_in_uv=x21*u_coord(ffff)+x22*v_coord(ffff)+x23*(1-u_coord(ffff)-v_coord(ffff)); %express the test points coords with uv coords
y2gauss_in_uv=y21*u_coord(ffff)+y22*v_coord(ffff)+y23*(1-u_coord(ffff)-v_coord(ffff));
z2gauss_in_uv=z21*u_coord(ffff)+z22*v_coord(ffff)+z23*(1-u_coord(ffff)-v_coord(ffff));
% calculate the distance norm for the biot savart law
distance_norm=((x1gauss_in_uv-x2gauss_in_uv)^2+(y1gauss_in_uv-y2gauss_in_uv)^2+(z1gauss_in_uv-z2gauss_in_uv)^2)^(-3/2);
inductance_matrix(aaaa,bbbb)=inductance_matrix(aaaa,bbbb)+dot_product*distance_norm*biot_savart_coeff*2*basis_elements(node_ind1).area(node_ind2)*gauss_weight(eeee)*gauss_weight(ffff);
end
end

end
end
end
end

end
end


%take care of the singularities of triangle self inductance
for aaaa=1:num_nodes
triangual_inductance_self_contribution=0;
for bbbb=num_triangles_per_node(aaaa)
for cccc=num_triangles_per_node(aaaa)
r_1=squeeze(basis_elements(aaaa).triangle_points_ABC(bbbb,:,1));
r_2=squeeze(basis_elements(aaaa).triangle_points_ABC(bbbb,:,2));
r_3=squeeze(basis_elements(aaaa).triangle_points_ABC(bbbb,:,3));
term_a=dot((r_3-r_1),(r_3-r_1));
term_b=dot((r_3-r_1),(r_3-r_2));
term_c=dot((r_3-r_2),(r_3-r_2));
dot_product=basis_elements(aaaa).current(bbbb,:)*basis_elements(aaaa).current(cccc,:)';
pre_factor=4*basis_elements(aaaa).area(bbbb)*basis_elements(aaaa).area(cccc)*biot_savart_coeff;

triangual_inductance_self_contribution=triangual_inductance_self_contribution+pre_factor*dot_product*( ...
1/(6*sqrt(term_a))*log( ( term_a-term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c)*(term_b+sqrt(term_a*term_c)) ) / ...
                                    ( ((-1)*term_b+sqrt(term_a*term_c))*((-1)*term_a+term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c)))) ...
+1/(6*sqrt(term_c))*log( ( (term_b+sqrt(term_a*term_c))* ((-1)*term_b+term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))) /...
                                        (( term_b-term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))*((-1)*term_b+sqrt(term_a*term_c)))) ...
+1/(6*sqrt(term_a-2*term_b+term_c))*log( ( (term_a-term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c))*((-1)*term_b+term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))) /...
                                                                        ( (term_b-term_c+sqrt(term_c)*sqrt(term_a-2*term_b+term_c))*((-1)*term_a+term_b+sqrt(term_a)*sqrt(term_a-2*term_b+term_c)))));
end
end
inductance_matrix(aaaa,aaaa)=triangual_inductance_self_contribution;
end


end
