
function [min_dist,near_points_a,min_ind_a,near_points_b,min_ind_b]=find_min_mutual_loop_distance(loop_a,loop_b,only_point_flag)
%Calculate the mutal nearest positions and segment indices between two
%loops 
% Copyright: Philipp Amrein, Uniklinik Freiburg Mai 2022

if ~only_point_flag

near_dists=zeros(1,size(loop_b.v,2)-1);
near_t_b=zeros(1,size(loop_b.v,2)-1);
near_points_b.v=zeros(3,size(loop_b.v,2)-1);
for test_point_ind=1:size(loop_b.v,2)-1
x1=repmat(loop_b.v(:,test_point_ind),[1 size(loop_a.v,2)]);
x2=repmat(loop_b.v(:,test_point_ind+1),[1 size(loop_a.v,2)]);
t=dot(loop_a.v-x1,x2-x1)./sum((x2-x1).^2,1);
t(t<0)=0;
t(t>1)=1;
all_near_points_b=x1+(x2-x1).*repmat(t,[3 1]);
all_dists=vecnorm(all_near_points_b-loop_a.v);
[~,min_ind_b]=min(all_dists);
near_dists(test_point_ind)=all_dists(min_ind_b);
near_t_b(test_point_ind)=t(min_ind_b);
near_points_b.v(:,test_point_ind)=all_near_points_b(:,min_ind_b);
end
 %find the nearest point within loop b
[min_dist,min_ind_b]=min(near_dists);
near_points_b.v=near_points_b.v(:,min_ind_b);
near_t_b=near_t_b(min_ind_b);
near_points_b.uv=loop_b.uv(:,min_ind_b)+(loop_b.uv(:,min_ind_b+1)-loop_b.uv(:,min_ind_b)).*repmat(near_t_b,[2 1]);


%find the nearest respective point within loop a
x1=loop_a.v(:,1:end-1);
x2=loop_a.v(:,2:end);
t=dot(near_points_b.v-x1,x2-x1)./sum((x2-x1).^2,1);
t(t<0)=0;
t(t>1)=1;
all_near_points_a=x1+(x2-x1).*repmat(t,[3 1]);
all_dists=vecnorm(all_near_points_a-near_points_b.v);
[~,min_ind_a]=min(all_dists);
near_points_a.v=all_near_points_a(:,min_ind_a);
near_points_a.uv=loop_a.uv(:,min_ind_a)+(loop_a.uv(:,min_ind_a+1)-loop_a.uv(:,min_ind_a)).*repmat(t(min_ind_a),[2 1]);


else
%find mutual nearest points; neglecting positions between points

min_test_ind=zeros(1,size(loop_b.v,2));
min_dist_ind=zeros(1,size(loop_b.v,2));
for test_point_ind=1:size(loop_b.v,2)
[min_test_ind(test_point_ind),min_dist_ind(test_point_ind)]=min(vecnorm(loop_a.v-repmat(loop_b.v(:,test_point_ind),[1 size(loop_a.v,2)])));
end

[min_dist,min_ind_b]=min(min_test_ind);
[~,min_ind_a]=min(vecnorm(loop_a.v-loop_b.v(:,min_ind_b)));

near_points_a.v=loop_a.v(:,min_ind_a);
near_points_a.uv=loop_a.uv(:,min_ind_a);

near_points_b.v=loop_b.v(:,min_ind_b);
near_points_b.uv=loop_b.uv(:,min_ind_b);



end


% figure;
% hold on;
% plot3(loop_a.v(1,:),loop_a.v(2,:),loop_a.v(3,:),'g');
% plot3(loop_b.v(1,:),loop_b.v(2,:),loop_b.v(3,:),'b');
% scatter3(near_points_b.v(1),near_points_b.v(2),near_points_b.v(3),'filled');
% scatter3(near_points_a.v(1),near_points_a.v(2),near_points_a.v(3),'filled');
% axis equal;
% hold off;
% 
% 
% figure;
% hold on;
% plot(loop_a.uv(1,:),loop_a.uv(2,:),'g');
% plot(loop_b.uv(1,:),loop_b.uv(2,:),'b');
% scatter(near_points_b.uv(1),near_points_b.uv(2),'filled');
% scatter(near_points_a.uv(1),near_points_a.uv(2),'filled');
% axis equal;
% hold off;


% % % function [min_dist_out,min_ind_loop_a.v,min_ind_loop_b.v]=find_shortest_path_between_loops(loop_a.v,loop_b.v)
% % % %find the  points with shortest straight (3D coords) between two
% % % %loops
% % % 
% % % 
% % % %to avoid memory problems splitt the curve into several parts
% % % max_track_part_length=2000;
% % % 
% % % num_parts_a=ceil(size(loop_a.v,2)/max_track_part_length);
% % % num_parts_b=ceil(size(loop_b.v,2)/max_track_part_length);
% % % 
% % % inds_a=1:size(loop_a.v,2);
% % % inds_b=1:size(loop_b.v,2);
% % % 
% % % %Split the first track a
% % % if num_parts_a>1
% % % loop_a.v_parts(num_parts_a).coords=[];
% % % loop_a.v_parts(num_parts_a).inds=[]; 
% % % x_inds=[1:max_track_part_length:size(loop_a.v,2) size(loop_a.v,2)];
% % % if x_inds(end-1)==x_inds(end)
% % % x_inds(end)=[];
% % % end
% % % for part_ind=1:num_parts_a
% % % loop_a.v_parts(part_ind).coords=loop_a.v(:,x_inds(part_ind):x_inds(part_ind+1)-1);
% % % loop_a.v_parts(part_ind).inds=inds_a(x_inds(part_ind):x_inds(part_ind+1)-1);
% % % end
% % % loop_a.v_parts(part_ind).coords=[loop_a.v_parts(part_ind).coords loop_a.v(:,end)];
% % % loop_a.v_parts(part_ind).inds=[loop_a.v_parts(part_ind).inds size(loop_a.v,2)];
% % % else
% % % loop_a.v_parts.coords=loop_a.v;
% % % loop_a.v_parts.inds=inds_a;
% % % end
% % % 
% % % %Split the second track b
% % % if num_parts_b>1
% % % loop_b.v_parts(num_parts_b).coords=[];
% % % loop_b.v_parts(num_parts_b).inds=[]; 
% % % x_inds=[1:max_track_part_length:size(loop_b.v,2) size(loop_b.v,2)];
% % % if x_inds(end-1)==x_inds(end)
% % % x_inds(end)=[];
% % % end
% % % for part_ind=1:num_parts_b
% % % loop_b.v_parts(part_ind).coords=loop_b.v(:,x_inds(part_ind):x_inds(part_ind+1)-1);
% % % loop_b.v_parts(part_ind).inds=inds_b(x_inds(part_ind):x_inds(part_ind+1)-1);
% % % end
% % % loop_b.v_parts(part_ind).coords=[loop_b.v_parts(part_ind).coords loop_b.v(:,end)];
% % % loop_b.v_parts(part_ind).inds=[loop_b.v_parts(part_ind).inds size(loop_b.v,2)];
% % % else
% % % loop_b.v_parts.coords=loop_b.v;
% % % loop_b.v_parts.inds=inds_b;
% % % end
% % % 
% % % part_distances=zeros(num_parts_a*num_parts_b,3);
% % % run_ind=1;
% % % %Calculate the distances
% % % for ind_a=1:num_parts_a
% % % for ind_b=1:num_parts_b
% % % loop_a.v_x_grid=repmat(loop_a.v_parts(ind_a).coords(1,:),[size(loop_b.v_parts(ind_b).coords,2) 1]);
% % % loop_a.v_y_grid=repmat(loop_a.v_parts(ind_a).coords(2,:),[size(loop_b.v_parts(ind_b).coords,2) 1]);
% % % loop_a.v_z_grid=repmat(loop_a.v_parts(ind_a).coords(3,:),[size(loop_b.v_parts(ind_b).coords,2) 1]);
% % % loop_b.v_x_grid=repmat(loop_b.v_parts(ind_b).coords(1,:)',[1 size(loop_a.v_parts(ind_a).coords,2)]);
% % % loop_b.v_y_grid=repmat(loop_b.v_parts(ind_b).coords(2,:)',[1 size(loop_a.v_parts(ind_a).coords,2)]);
% % % loop_b.v_z_grid=repmat(loop_b.v_parts(ind_b).coords(3,:)',[1 size(loop_a.v_parts(ind_a).coords,2)]);
% % % dist_mat=sqrt((loop_a.v_x_grid-loop_b.v_x_grid).^2+(loop_a.v_y_grid-loop_b.v_y_grid).^2+(loop_a.v_z_grid-loop_b.v_z_grid).^2);
% % % [min_dist_out,min_sub_ind]=min(dist_mat(:));
% % % [min_ind_loop_b.v,min_ind_loop_a.v] = ind2sub(size(dist_mat),min_sub_ind);
% % % part_distances(run_ind,:)=[loop_a.v_parts(ind_a).inds(min_ind_loop_a.v) loop_b.v_parts(ind_b).inds(min_ind_loop_b.v) min_dist_out];
% % % run_ind=run_ind+1;
% % % end
% % % end
% % % 
% % % [min_dist_out,min_ind]=min(part_distances(:,3));
% % % min_ind_loop_a.v=part_distances(min_ind,1);
% % % min_ind_loop_b.v=part_distances(min_ind,2);
% % % 
% % % % loop_a.v_x_grid=repmat(loop_a.v(1,:),[size(loop_b.v,2) 1]);
% % % % loop_a.v_y_grid=repmat(loop_a.v(2,:),[size(loop_b.v,2) 1]);
% % % % loop_a.v_z_grid=repmat(loop_a.v(3,:),[size(loop_b.v,2) 1]);
% % % % loop_b.v_x_grid=repmat(loop_b.v(1,:)',[1 size(loop_a.v,2)]);
% % % % loop_b.v_y_grid=repmat(loop_b.v(2,:)',[1 size(loop_a.v,2)]);
% % % % loop_b.v_z_grid=repmat(loop_b.v(3,:)',[1 size(loop_a.v,2)]);
% % % % dist_mat=sqrt((loop_a.v_x_grid-loop_b.v_x_grid).^2+(loop_a.v_y_grid-loop_b.v_y_grid).^2+(loop_a.v_z_grid-loop_b.v_z_grid).^2);
% % % % [min_dist_out,min_sub_ind]=min(dist_mat(:));
% % % % [min_ind_loop_b.v,min_ind_loop_a.v] = ind2sub(size(dist_mat),min_sub_ind);
% % % end

end