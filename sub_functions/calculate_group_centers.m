function coil_parts=calculate_group_centers(coil_parts)
%calculate group centers


coil_parts(numel(coil_parts)).group_centers=[];

for part_ind=1:numel(coil_parts)




total_center=mean(coil_parts(part_ind).coil_mesh.uv,2);

group_centers_2d=zeros(2,numel(coil_parts(part_ind).groups));
total_group_center.uv=zeros(2,numel(coil_parts(part_ind).groups));
total_group_center.v=zeros(3,numel(coil_parts(part_ind).groups));


for group_ind=1:numel(coil_parts(part_ind).groups)
point_sum_uv=[];
point_sum_v=[];
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
point_sum_uv=[point_sum_uv coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv];
point_sum_v=[point_sum_v coil_parts(part_ind).groups(group_ind).loops(loop_ind).v];  
end
total_group_center.uv(:,group_ind)=mean(point_sum_uv,2);
total_group_center.v(:,group_ind)=mean(point_sum_v,2);
inner_center=mean(coil_parts(part_ind).groups(group_ind).loops(end).uv,2);
%check if the total group center is within the most inner loop of the group
inner_test_loop=[coil_parts(part_ind).groups(group_ind).loops(end).uv-mean(coil_parts(part_ind).groups(group_ind).loops(end).uv,2)].*0.9+mean(coil_parts(part_ind).groups(group_ind).loops(end).uv,2);
%total_group_center_is_in=inpolygon(total_group_center.uv(1,group_ind),total_group_center.uv(2,group_ind),inner_test_loop(1,:),inner_test_loop(2,:));
total_group_center_is_in=check_mutual_loop_inclusion(total_group_center.uv(:,group_ind),inner_test_loop);
if total_group_center_is_in==1
% total group center and inner center is within inner loop
group_centers_2d(:,group_ind)=total_group_center.uv(:,group_ind);
else
% [~,min_ind]=min(inner_test_loop-total_center);
% inner_line_point=inner_test_loop(:,min_ind);    
scale_ind=1000;
cut_line_x=[inner_center(1)+(total_center(1)-inner_center(1))*scale_ind inner_center(1)-(total_center(1)-inner_center(1))*scale_ind];
cut_line_y=[inner_center(2)+(total_center(2)-inner_center(2))*scale_ind inner_center(2)-(total_center(2)-inner_center(2))*scale_ind];
intersection_points=find_segment_intersections(coil_parts(part_ind).groups(group_ind).loops(end).uv,[cut_line_x; cut_line_y]);
line_cut_inner_total_x=intersection_points.uv(1,:)';
line_cut_inner_total_y=intersection_points.uv(2,:)';
if isempty(line_cut_inner_total_x)
group_centers_2d(:,group_ind)=inner_center;
else 
%sort the cut points for their distance to the inner center
[~,min_ind]=sort(vecnorm([line_cut_inner_total_x line_cut_inner_total_y]'-total_center));
inner_cut_point=[mean(line_cut_inner_total_x(min_ind([1 2]))) mean(line_cut_inner_total_y(min_ind([1 2])))]';
group_centers_2d(:,group_ind)=inner_cut_point;
% group_centers_2d(:,group_ind)=total_center+([line_cut_inner_total_x(min_ind(1)) line_cut_inner_total_y(min_ind(1))]'-total_center).*(1);
% group_centers_2d(:,group_ind)=[mean([line_cut_inner_total_x(min_ind(1)) line_cut_inner_total_x(min_ind(2))]);mean([line_cut_inner_total_y(min_ind(1)) line_cut_inner_total_y(min_ind(2))])];
end
end
end 
%     if inner_center_is_in==1
%     [~,min_ind]=min(vecnorm(coil_parts(part_ind).groups(group_ind).loops(end).uv-total_group_center.uv));
%     inner_line_point=coil_parts(part_ind).groups(group_ind).loops(end).uv(:,min_ind);    
%     scale_ind=1000;
%     cut_line_x=[inner_line_point(1)+(total_group_center.uv(1)-inner_line_point(1))*scale_ind inner_line_point(1)-(total_group_center.uv(1)-inner_line_point(1))*scale_ind];
%     cut_line_y=[inner_line_point(2)+(total_group_center.uv(2)-inner_line_point(2))*scale_ind inner_line_point(2)-(total_group_center.uv(2)-inner_line_point(2))*scale_ind];
%     [line_cut_inner_total_x,line_cut_inner_total_y] = polyxpoly(coil_parts(part_ind).groups(group_ind).loops(end).uv(1,:),coil_parts(part_ind).groups(group_ind).loops(end).uv(2,:),cut_line_x,cut_line_y);
%     if ~isempty(line_cut_inner_total_x)
%     group_centers_2d(:,group_ind)=[mean([line_cut_inner_total_x(1) line_cut_inner_total_x(2)]);mean([line_cut_inner_total_y(1) line_cut_inner_total_y(2)])];
%     else
%     group_centers_2d(:,group_ind)=inner_center;
%     end
%     else %try to fix the inner center so that it is within the inner loop
%     [~,min_ind]=min(vecnorm(coil_parts(part_ind).groups(group_ind).loops(end).uv-inner_center));
%     inner_line_point=coil_parts(part_ind).groups(group_ind).loops(end).uv(:,min_ind);    
%     scale_ind=1000;
%     cut_line_x=[inner_line_point(1)+(inner_center(1)-inner_line_point(1))*scale_ind inner_line_point(1)-(inner_center(1)-inner_line_point(1))*scale_ind];
%     cut_line_y=[inner_line_point(2)+(inner_center(2)-inner_line_point(2))*scale_ind inner_line_point(2)-(inner_center(2)-inner_line_point(2))*scale_ind];
%     [line_cut_inner_total_x,line_cut_inner_total_y] = polyxpoly(coil_parts(part_ind).groups(group_ind).loops(end).uv(1,:),coil_parts(part_ind).groups(group_ind).loops(end).uv(2,:),cut_line_x,cut_line_y);
%     if ~isempty(line_cut_inner_total_x)
%     group_centers_2d(:,group_ind)=[mean([line_cut_inner_total_x(1) line_cut_inner_total_x(2)]);mean([line_cut_inner_total_y(1) line_cut_inner_total_y(2)])];
%     else
%     group_centers_2d(:,group_ind)=inner_center;
%     end
%     end

% if ~isempty(inner_center_is_in) & isempty(total_group_center_is_in)
% %  inner center is within inner loop but not the total group center
% else
% % both center are not within the inner loop
% group_centers_2d(:,group_ind)=inner_center;
% end
%end

% group_size=[max(coil_parts(part_ind).groups(3).loops(1).uv(1,:)) min(coil_parts(part_ind).groups(3).loops(1).uv(1,:)); max(coil_parts(part_ind).groups(3).loops(1).uv(2,:)) min(coil_parts(part_ind).groups(3).loops(1).uv(2,:))];
% group_diameter=vecnorm(group_size(:,2)-group_size(:,1));
% [~,min_ind]=min(vecnorm(coil_parts(part_ind).groups(tttt).loops(end).uv-inner_center));
% point_b=coil_parts(part_ind).groups(tttt).loops(end).uv(:,min_ind);
% point_c=inner_center+(point_b-inner_center)./norm((point_b-inner_center))*group_diameter;
% point_d=inner_center-(point_b-inner_center)./norm((point_b-inner_center))*group_diameter;
% [xi,yi] = polyxpoly(coil_parts(part_ind).groups(tttt).loops(end).uv(1,:),coil_parts(part_ind).groups(tttt).loops(end).uv(2,:),[point_c(1) point_d(1)],[point_c(2) point_d(2)]);
% group_centers_2d(:,tttt)=[mean(xi);mean(yi)];
%end


%Set the centers, consider that the possibility of non-mesh points
group_centers_3d=zeros(3,size(group_centers_2d,2));
planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices');
for rrrr=1:numel(coil_parts(part_ind).groups)
%set centers outside the 2D mesh in the center of the 3D volume
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,group_centers_2d(1,rrrr),group_centers_2d(2,rrrr));
if ~isnan(target_triangle)
group_centers_3d(:,rrrr) = barycentricToCartesian(curved_mesh,target_triangle,bary_centric_coord)';
else
group_centers_3d(:,rrrr)=total_group_center.v(:,rrrr);                     
end
end

coil_parts(part_ind).group_centers.uv=group_centers_2d;
coil_parts(part_ind).group_centers.v=group_centers_3d;


end

end  
