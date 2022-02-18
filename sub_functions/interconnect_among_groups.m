function coil_parts=interconnect_among_groups(coil_parts,opening_gap)
% interconect among groups genereting a single wire track

%local parameters
up_sample_factor=5;
%ratio_criteria_for_return_path_usage=1.5;

coil_parts(numel(coil_parts)).opening_cuts_among_groups.uv=[]; %save the opening cuts for plotting
coil_parts(numel(coil_parts)).wire_path=[];

for part_ind=1:numel(coil_parts)


winding_distance = calc_avg_winding_distance(coil_parts(part_ind).groups);
winding_distance(isnan(winding_distance))=0;


coil_parts(part_ind).opening_cuts_among_groups.uv=[]; %save the oping cuts for plotting

%upsample the loop with points to increase precision in distance evaluation
clear  upsampled_groups
upsampled_groups(numel(coil_parts(part_ind).connected_group)).v=[];
upsampled_groups(numel(coil_parts(part_ind).connected_group)).uv=[];
for group_ind=1:numel(coil_parts(part_ind).connected_group)
buff_track=upsample_loop(coil_parts(part_ind).connected_group(group_ind),up_sample_factor);
upsampled_groups(group_ind).uv = getfield( buff_track,'uv');
upsampled_groups(group_ind).v = getfield( buff_track,'v');
end



%%%% Group fusions%%%%

level_hierachy=cellfun(@numel,coil_parts(part_ind).level_positions);

for level_ind=max(level_hierachy):-1:0

levels_to_process=find(level_hierachy==level_ind);

%interconnect one level into a single group
for single_level_ind=1:numel(levels_to_process)
current_level=levels_to_process(single_level_ind);
groups_to_connect=coil_parts(part_ind).group_levels{current_level};
    

%select the current host group of the level
is_enclosing=zeros(1,numel(groups_to_connect));
if ~isempty(coil_parts(part_ind).level_positions{current_level})
current_top_group=coil_parts(part_ind).level_positions{current_level}(end);
groups_to_connect=[groups_to_connect current_top_group];
is_enclosing=[is_enclosing 1];
end

%make the n-1 interconnections in an optimal way resulating in one track
%for that level:
num_connections_to_do=numel(groups_to_connect)-1;
for connect_ind=1:num_connections_to_do
    
    

%get the tracks to connect
grouptracks_to_connect=upsampled_groups(groups_to_connect);

%select the return paths of those interconnected groups for later
%groups_return_paths=return_paths(groups_to_connect);

%find the mininmal distance positions between the groups
%and the points with mininmal distance
min_group_dists=zeros(numel(groups_to_connect));
min_pos_group1=cell(numel(groups_to_connect));
min_pos_group2=cell(numel(groups_to_connect));
for ind1=1:numel(groups_to_connect)
for ind2=1:numel(groups_to_connect)
if ind2<ind1
[min_group_dists(ind1,ind2),min_ind_first,min_ind_second]=find_shortest_path_between_loops(grouptracks_to_connect(ind1).v,grouptracks_to_connect(ind2).v);
min_pos_group1{ind1,ind2}=grouptracks_to_connect(ind1).uv(:,min_ind_first);
min_pos_group2{ind1,ind2}=grouptracks_to_connect(ind2).uv(:,min_ind_second);
end
end
end
min_group_dists=min_group_dists+min_group_dists';
min_group_dists(min_group_dists==0)=Inf;
min_pos_group1=cellfun(@(x,y) [x y], min_pos_group1, min_pos_group1', 'UniformOutput', false);
min_pos_group2=cellfun(@(x,y) [x y], min_pos_group2, min_pos_group2', 'UniformOutput', false);



%select the pair of groups with shortest respective distance
min_dist_couple=min_group_dists==min(min_group_dists(:));
min_dist_couple=find(min_dist_couple);
[couple_group1,couple_group2]=ind2sub(size(min_group_dists),min_dist_couple(1));


% %decide wether return paths will be employed for the group fusion
% % (here check either return of group1 or return of group2)
% use_first_return=use_return_path(couple_group2,couple_group1);
% use_second_return=use_return_path(couple_group1,couple_group2);
% 
% if ~use_first_return && ~use_second_return
% do_not_use_any_return=1;
% else
% do_not_use_any_return=0;
% end

%create and prepare the opening points from the shortest distance points
raw_cut_point_1=min_pos_group1{couple_group1,couple_group2};
raw_cut_point_2=min_pos_group2{couple_group1,couple_group2};
raw_cut_point_dist=min(min_group_dists(:));

middel_point=(raw_cut_point_1+raw_cut_point_2)./2;
track_1=grouptracks_to_connect(couple_group1).uv;
track_2=grouptracks_to_connect(couple_group2).uv;
[~,min_dists1]=sort(vecnorm(middel_point-track_1));
[~,min_dists2]=sort(vecnorm(middel_point-track_2));
min_dists1_shift=circshift(1:size(track_1,2),-1);
min_dists2_shift=circshift(1:size(track_2,2),-1);
longitudinal_vec1=track_1(:,min_dists1_shift(min_dists1(1)))-track_1(:,min_dists1(1));
longitudinal_vec2=track_2(:,min_dists2_shift(min_dists2(1)))-track_2(:,min_dists2(1));
if dot(longitudinal_vec1,longitudinal_vec2)<0
longitudinal_vec2=longitudinal_vec2.*(-1);
end
%make sure the min_dists are in the same winding!
longitudinal_avg=(longitudinal_vec1+longitudinal_vec2)./2;
% if  abs(mod(min_dists1(2),size(track_1,2))-mod(min_dists1(1),size(track_1,2)))>2 || abs(mod(min_dists2(2),size(track_2,2))-mod(min_dists2(1),size(track_2,2)))>2 
% transverse_vector=raw_cut_point_2-raw_cut_point_1;
% else
 transverse_vector=[longitudinal_avg(2);longitudinal_avg(1).*(-1)];
%end

transverse_vector=transverse_vector./vecnorm(transverse_vector).*opening_gap*10;
[cut1_x,cut1_y] = polyxpoly(track_1(1,:),track_1(2,:),[middel_point(1)+transverse_vector(1) middel_point(1)-transverse_vector(1)],...
[middel_point(2)+transverse_vector(2) middel_point(2)-transverse_vector(2)]);
[cut2_x,cut2_y] = polyxpoly(track_2(1,:),track_2(2,:),[middel_point(1)+transverse_vector(1) middel_point(1)-transverse_vector(1)],...
[middel_point(2)+transverse_vector(2) middel_point(2)-transverse_vector(2)]);
[~,min_inds_cut1]=sort(vecnorm(middel_point- [cut1_x cut1_y]'));
[~,min_inds_cut2]=sort(vecnorm(middel_point- [cut2_x cut2_y]'));

if ~isempty(cut1_x) && ~isempty(cut2_x) 
adj_cutpoint_1=[cut1_x(min_inds_cut1(1));cut1_y(min_inds_cut1(1))];
adj_cutpoint_2=[cut2_x(min_inds_cut2(1));cut2_y(min_inds_cut2(1))];
is_worse=vecnorm(adj_cutpoint_1-adj_cutpoint_2)>vecnorm(raw_cut_point_1-raw_cut_point_2)*1.3;
if  ~is_worse
cutpoint_1=adj_cutpoint_1;
cutpoint_2=adj_cutpoint_2;
else
cutpoint_1=raw_cut_point_1;
cutpoint_2=raw_cut_point_2;
end
else
cutpoint_1=raw_cut_point_1;
cutpoint_2=raw_cut_point_2;
end
  



% %mark the shortest group to group distances
% shortest_group_dist=min_group_dists==(repmat(min(min_group_dists,[],1),[size(min_group_dists,1) 1]))';
% %find the mininmal distance positions between the groups and the return
% %paths  (these distances will be used in the decision wether the cuts of the return paths will be employed)
% min_to_return_dists=zeros(numel(groups_to_connect));
% min_pos_group=cell(numel(groups_to_connect));
% min_pos_return=cell(numel(groups_to_connect));
% for ind1=1:numel(groups_to_connect)
% for ind2=1:numel(groups_to_connect)
% if ind1~=ind2
% [min_to_return_dists(ind1,ind2),min_to_retrun_ind_first,min_to_retrun_ind_second]=find_shortest_path_between_loops(grouptracks_to_connect(ind1).v,groups_return_paths(ind2).v);
% min_pos_group{ind1,ind2}=grouptracks_to_connect(ind1).uv(:,min_to_retrun_ind_first);
% min_pos_return{ind1,ind2}=groups_return_paths(ind2).uv(:,min_to_retrun_ind_second);
% end
% end
% end
% min_to_return_dists(min_to_return_dists==0)=Inf;
% %mark the return paths with more than one points
% return_path_has_direction=arrayfun(@(x) size(groups_return_paths(x).uv,2)>1 , 1:numel(groups_return_paths));
% return_path_has_direction=repmat(return_path_has_direction,[numel(return_path_has_direction) 1]);
% %calc the ratio of the distances between groups and groups and return paths
% connecetion_dist_ratio=min_to_return_dists./min_group_dists;
% %decide wether the return paths are used for the group connection
% use_return_path=connecetion_dist_ratio<ratio_criteria_for_return_path_usage & shortest_group_dist & return_path_has_direction ;



%do simple connection
[fused_group,coil_parts(part_ind).opening_cuts_among_groups(end+1).uv]=fuse_neighbouring_groups(grouptracks_to_connect(couple_group1),grouptracks_to_connect(couple_group2),...
    cutpoint_1,cutpoint_2,opening_gap,coil_parts(part_ind).coil_mesh,(abs(winding_distance(couple_group1))+abs(winding_distance(couple_group1)))/2,raw_cut_point_dist);


% if do_not_use_any_return
%     %do simple connection
%     [fused_group,opening_cuts(end+1).uv]=fuse_neighbouring_groups(grouptracks_to_connect(couple_group1),grouptracks_to_connect(couple_group2),...
%         min_pos_group1{couple_group1,couple_group2},min_pos_group2{couple_group1,couple_group2},opening_gap,parameterized_mesh);
%     else
%     %use a return path
%     if use_second_return
%     [fused_group,opening_cuts(end+1).uv]=fuse_group_by_neighbouring_return_path(grouptracks_to_connect(couple_group1),grouptracks_to_connect(couple_group2),opening_gap,groups_return_paths(couple_group2).uv,parameterized_mesh);
%     else    %use a return path of group 1
%     [fused_group,opening_cuts(end+1).uv]=fuse_group_by_neighbouring_return_path(grouptracks_to_connect(couple_group2),grouptracks_to_connect(couple_group1),opening_gap,groups_return_paths(couple_group1).uv,parameterized_mesh);
%     end
%     
% end

% %equlize the track
% fused_group = upsample_loop_uv(fused_group,1);

%overwrite fused track into both group point arrays
%delete in all the data containers one of the connected groups to avoid
%redundancy
%do not select here the host level group!
if is_enclosing(couple_group1)
upsampled_groups(groups_to_connect(couple_group1)).uv=fused_group;
[upsampled_groups(groups_to_connect(couple_group1)).v,upsampled_groups(groups_to_connect(couple_group1)).uv]=uv_to_xyz(fused_group,triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv'),triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'));
%upsampled_groups(groups_to_connect(couple_group2))=[];
is_enclosing(couple_group2)=[];
groups_to_connect(couple_group2)=[];
else
upsampled_groups(groups_to_connect(couple_group2)).uv=fused_group;
[upsampled_groups(groups_to_connect(couple_group2)).v,upsampled_groups(groups_to_connect(couple_group2)).uv]=uv_to_xyz(fused_group,triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv'),triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'));
%upsampled_groups(groups_to_connect(couple_group1))=[];
is_enclosing(couple_group1)=[];
groups_to_connect(couple_group1)=[];
end




end


end

end

%Select the full track as the final return
[~,is_final_ind]=max(arrayfun(@(x) size(upsampled_groups(x).v,2), 1:numel(upsampled_groups)));
full_track=upsampled_groups(is_final_ind(1));
%Shift the open ends to the boundarie of the coil
[~,min_ind]=min(vecnorm(full_track.uv-mean(full_track.uv,2)));
full_track.v=circshift(full_track.v,(-1)*min_ind,2);
full_track.uv=circshift(full_track.uv,(-1)*min_ind,2);

coil_parts(part_ind).opening_cuts_among_groups(1)=[];

%Assign the outputs
coil_parts(part_ind).wire_path=full_track;

end



end

function loop_to_loop_distance = calc_avg_winding_distance(group_container)
%calc the average distance between the contour lines for all groups
loop_to_loop_distance=zeros(1,numel(group_container));
for group_ind=1:numel(group_container)
loop_radi=numel(group_container(group_ind).loops);
for loop_ind=1:numel(group_container(group_ind).loops)
loop_radi(loop_ind)=mean(vecnorm(group_container(group_ind).loops(loop_ind).uv-mean(group_container(group_ind).loops(loop_ind).uv)));
end
loop_to_loop_distance(group_ind)=mean(abs(diff(loop_radi)));
end
end

function points_out = upsample_loop(points_in,upsample_factor)
%limit the point to point distane in the 2d representation
points_out=points_in;
seg_lengths=vecnorm(points_out.v(:,2:end)-points_out.v(:,1:end-1));
target_segment_length=mean(seg_lengths)/upsample_factor;
%for each segment calculate the of points to insert
point_num_to_insert=floor(seg_lengths./target_segment_length)-1;
segment_cells_v=arrayfun(@(x) [points_in.v(:,x-1) points_in.v(:,x)],2:size(points_in.v,2),'UniformOutput',false);
segment_cells_uv=arrayfun(@(x) [points_in.uv(:,x-1) points_in.uv(:,x)],2:size(points_in.uv,2),'UniformOutput',false);
for seg_ind=1:numel(seg_lengths)
seg_vec_v=[segment_cells_v{seg_ind}(:,2)-segment_cells_v{seg_ind}(:,1)]./(point_num_to_insert(seg_ind)+1);
points_to_insert_v=segment_cells_v{seg_ind}(:,1)+seg_vec_v.*(1:point_num_to_insert(seg_ind));
segment_cells_v{seg_ind}=[segment_cells_v{seg_ind}(:,1) points_to_insert_v];
seg_vec_uv=[segment_cells_uv{seg_ind}(:,2)-segment_cells_uv{seg_ind}(:,1)]./(point_num_to_insert(seg_ind)+1);
points_to_insert_uv=segment_cells_uv{seg_ind}(:,1)+seg_vec_uv.*(1:point_num_to_insert(seg_ind));
segment_cells_uv{seg_ind}=[segment_cells_uv{seg_ind}(:,1) points_to_insert_uv];
end
points_out.v=[segment_cells_v{:} points_in.v(:,end)];
points_out.uv=[segment_cells_uv{:} points_in.uv(:,end)];
end



function [min_dist_out,min_ind_loop_a,min_ind_loop_b]=find_shortest_path_between_loops_2d(loop_a,loop_b)
%find the  points with shortest straight (3D coords) between two
%loops
loop_a_u_grid=repmat(loop_a(1,:),[size(loop_b,2) 1]);
loop_a_v_grid=repmat(loop_a(2,:),[size(loop_b,2) 1]);
loop_b_u_grid=repmat(loop_b(1,:)',[1 size(loop_a,2)]);
loop_b_v_grid=repmat(loop_b(2,:)',[1 size(loop_a,2)]);
dist_mat=sqrt((loop_a_u_grid-loop_b_u_grid).^2+(loop_a_v_grid-loop_b_v_grid).^2);
[min_dist_out,min_sub_ind]=min(dist_mat(:));
[min_ind_loop_b,min_ind_loop_a] = ind2sub(size(dist_mat),min_sub_ind);
end

function [min_dist_out,min_ind_loop_a,min_ind_loop_b]=find_shortest_path_between_loops(loop_a,loop_b)
%find the  points with shortest straight (3D coords) between two
%loops


%to avoid memory problems splitt the curve into several parts
max_track_part_length=2000;

num_parts_a=ceil(size(loop_a,2)/max_track_part_length);
num_parts_b=ceil(size(loop_b,2)/max_track_part_length);

inds_a=1:size(loop_a,2);
inds_b=1:size(loop_b,2);

%Split the first track a
if num_parts_a>1
loop_a_parts(num_parts_a).coords=[];
loop_a_parts(num_parts_a).inds=[]; 
x_inds=[1:max_track_part_length:size(loop_a,2) size(loop_a,2)];
if x_inds(end-1)==x_inds(end)
x_inds(end)=[];
end
for part_ind=1:num_parts_a
loop_a_parts(part_ind).coords=loop_a(:,x_inds(part_ind):x_inds(part_ind+1)-1);
loop_a_parts(part_ind).inds=inds_a(x_inds(part_ind):x_inds(part_ind+1)-1);
end
loop_a_parts(part_ind).coords=[loop_a_parts(part_ind).coords loop_a(:,end)];
loop_a_parts(part_ind).inds=[loop_a_parts(part_ind).inds size(loop_a,2)];
else
loop_a_parts.coords=loop_a;
loop_a_parts.inds=inds_a;
end

%Split the second track b
if num_parts_b>1
loop_b_parts(num_parts_b).coords=[];
loop_b_parts(num_parts_b).inds=[]; 
x_inds=[1:max_track_part_length:size(loop_b,2) size(loop_b,2)];
if x_inds(end-1)==x_inds(end)
x_inds(end)=[];
end
for part_ind=1:num_parts_b
loop_b_parts(part_ind).coords=loop_b(:,x_inds(part_ind):x_inds(part_ind+1)-1);
loop_b_parts(part_ind).inds=inds_b(x_inds(part_ind):x_inds(part_ind+1)-1);
end
loop_b_parts(part_ind).coords=[loop_b_parts(part_ind).coords loop_b(:,end)];
loop_b_parts(part_ind).inds=[loop_b_parts(part_ind).inds size(loop_b,2)];
else
loop_b_parts.coords=loop_b;
loop_b_parts.inds=inds_b;
end

part_distances=zeros(num_parts_a*num_parts_b,3);
run_ind=1;
%Calculate the distances
for ind_a=1:num_parts_a
for ind_b=1:num_parts_b
loop_a_x_grid=repmat(loop_a_parts(ind_a).coords(1,:),[size(loop_b_parts(ind_b).coords,2) 1]);
loop_a_y_grid=repmat(loop_a_parts(ind_a).coords(2,:),[size(loop_b_parts(ind_b).coords,2) 1]);
loop_a_z_grid=repmat(loop_a_parts(ind_a).coords(3,:),[size(loop_b_parts(ind_b).coords,2) 1]);
loop_b_x_grid=repmat(loop_b_parts(ind_b).coords(1,:)',[1 size(loop_a_parts(ind_a).coords,2)]);
loop_b_y_grid=repmat(loop_b_parts(ind_b).coords(2,:)',[1 size(loop_a_parts(ind_a).coords,2)]);
loop_b_z_grid=repmat(loop_b_parts(ind_b).coords(3,:)',[1 size(loop_a_parts(ind_a).coords,2)]);
dist_mat=sqrt((loop_a_x_grid-loop_b_x_grid).^2+(loop_a_y_grid-loop_b_y_grid).^2+(loop_a_z_grid-loop_b_z_grid).^2);
[min_dist_out,min_sub_ind]=min(dist_mat(:));
[min_ind_loop_b,min_ind_loop_a] = ind2sub(size(dist_mat),min_sub_ind);
part_distances(run_ind,:)=[loop_a_parts(ind_a).inds(min_ind_loop_a) loop_b_parts(ind_b).inds(min_ind_loop_b) min_dist_out];
run_ind=run_ind+1;
end
end

[min_dist_out,min_ind]=min(part_distances(:,3));
min_ind_loop_a=part_distances(min_ind,1);
min_ind_loop_b=part_distances(min_ind,2);

% loop_a_x_grid=repmat(loop_a(1,:),[size(loop_b,2) 1]);
% loop_a_y_grid=repmat(loop_a(2,:),[size(loop_b,2) 1]);
% loop_a_z_grid=repmat(loop_a(3,:),[size(loop_b,2) 1]);
% loop_b_x_grid=repmat(loop_b(1,:)',[1 size(loop_a,2)]);
% loop_b_y_grid=repmat(loop_b(2,:)',[1 size(loop_a,2)]);
% loop_b_z_grid=repmat(loop_b(3,:)',[1 size(loop_a,2)]);
% dist_mat=sqrt((loop_a_x_grid-loop_b_x_grid).^2+(loop_a_y_grid-loop_b_y_grid).^2+(loop_a_z_grid-loop_b_z_grid).^2);
% [min_dist_out,min_sub_ind]=min(dist_mat(:));
% [min_ind_loop_b,min_ind_loop_a] = ind2sub(size(dist_mat),min_sub_ind);
end


function [fused_group,opening_cut]=fuse_neighbouring_groups(group_1,group_2,cutpoint_1,cutpoint_2,opening_gap,parameterized_mesh,winding_distance,raw_cut_point_dist)
%do simple connection



longitudinal_vector=cutpoint_1-cutpoint_2;
othorgonal_vector=[longitudinal_vector(2);longitudinal_vector(1).*(-1)];
othorgonal_vector=othorgonal_vector./norm(othorgonal_vector);
middle_point=(cutpoint_1+cutpoint_2)./2;
local_opening_gab=calc_local_opening_gab(parameterized_mesh,middle_point,othorgonal_vector,opening_gap,10);

% %create a circular opening cut
% opening_cut=[sin(0:(2*pi)/(20-1):2*pi); cos(0:(2*pi)/(20-1):2*pi)].*(local_opening_gab/2);
% opening_cut=opening_cut+middle_point;


longitudinal_vector=longitudinal_vector.*1.5;
othorgonal_vector=othorgonal_vector./(norm(othorgonal_vector)).*local_opening_gab;


%Check if a geodesic cut shape is necessary (in  case of bigger distances between the cut points)

opening_cut =[middle_point+longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2+othorgonal_vector/2];

% %create a circular opening cut
% opening_cut=[sin(0:(2*pi)/(20-1):2*pi); cos(0:(2*pi)/(20-1):2*pi)].*(local_opening_gab/2);
% opening_cut=opening_cut+middle_point;      
                            
                            
%circshift the inner group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_1.uv-cutpoint_1));
group_1.uv=circshift(group_1.uv,min_ind*(-1),2);

%circshift the outer group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_2.uv-cutpoint_2));
group_2.uv=circshift(group_2.uv,min_ind*(-1),2);



%     %add the "cut"points for a clean cut
%     [cut_x,cut_y] = polyxpoly(cut_points(1,:),cut_points(2,:),rectangle_cut(1,:),rectangle_cut(2,:));
%     opened_group{loop_num}=raw_loop;
%     %find out the order in which the cut points must inlcuded in the opened loop:
%     test_dists_1=vecnorm(raw_loop(:,1)-[cut_x(1);cut_y(1)]);
%     test_dists_2=vecnorm(raw_loop(:,1)-[cut_x(2);cut_y(2)]);
%     if    test_dists_1<test_dists_2
%         opened_group{loop_num}=[[cut_x(1) cut_y(1)]' opened_group{loop_num} [cut_x(2) cut_y(2)]' ] ;
%     else 
%         opened_group{loop_num}=[[cut_x(2) cut_y(2)]' opened_group{loop_num} [cut_x(1) cut_y(1)]' ] ;
%     end


%fuse both groups
%Check which fusing order is better:
track_combilength1=[group_1.uv group_2.uv];
track_combilength2=[group_2.uv group_1.uv];
track_combilength1=sum(vecnorm(track_combilength1(:,2:end)-track_combilength1(:,1:end-1)));
track_combilength2=sum(vecnorm(track_combilength2(:,2:end)-track_combilength2(:,1:end-1)));
if track_combilength1<track_combilength2
fused_group=[group_1.uv group_2.uv];
else
fused_group=[group_2.uv group_1.uv];
end
%Clear the points witihn the opening_cut to improve to connection
[in_cut_rect_ind,~] = inpolygon(fused_group(1,:),fused_group(2,:),opening_cut(1,:),opening_cut(2,:)); 
fused_group(:,in_cut_rect_ind)=[];


end
                            




function [fused_group,opening_cut]=fuse_group_by_neighbouring_return_path(group_1,group_2,opening_gap,return_points_group2,parameterized_mesh)
%find the cut rectangle using a return path of one of the groups


% find a new cut point by using the return path of the outer loop
return_points_group2=fliplr(return_points_group2);


return_points_check_points=return_points_group2;
[~,direction_point_ind1,~]=find_shortest_path_between_loops_2d(return_points_check_points,group_1.uv);
direction_point_1=return_points_check_points(:,direction_point_ind1);
return_points_check_points=return_points_check_points(:,setdiff(1:size(return_points_group2,2),direction_point_ind1));
[~,direction_point_ind2,~]=find_shortest_path_between_loops_2d(return_points_check_points,group_1.uv);
direction_point_2=return_points_check_points(:,direction_point_ind2);
cut_direction=direction_point_2-direction_point_1;
cut_direction=cut_direction./vecnorm(cut_direction)*10^7;
[cutx,cuty] = polyxpoly([ direction_point_1(1)+cut_direction(1) direction_point_1(1)-cut_direction(1)],...
[direction_point_1(2)+cut_direction(2) direction_point_1(2)-cut_direction(2)],...
group_1.uv(1,[1:end 1]),group_1.uv(2,[1:end 1]));
if ~isempty(cutx)
[~,min_ind,~]=find_shortest_path_between_loops_2d([cutx'; cuty'],return_points_group2);
cut_point_group1=[cutx(min_ind);cuty(min_ind)];
[~,cut_point_group2]=min(vecnorm(return_points_group2-cut_point_group1));
cut_point_group2=return_points_group2(:,cut_point_group2);
else
[~,~,direction_point_ind2]=find_shortest_path_between_loops_2d(return_points_check_points,group_1.uv);
direction_point_2=group_1.uv(:,direction_point_ind2);
cut_direction=direction_point_2-direction_point_1;
cut_direction=cut_direction./vecnorm(cut_direction)*10^7;
[cutx,cuty] = polyxpoly([ direction_point_1(1)+cut_direction(1) direction_point_1(1)-cut_direction(1)],...
[direction_point_1(2)+cut_direction(2) direction_point_1(2)-cut_direction(2)],...
group_1.uv(1,[1:end 1]),group_1.uv(2,[1:end 1]));
[~,min_ind,~]=find_shortest_path_between_loops_2d([cutx'; cuty'],return_points_group2);
cut_point_group1=[cutx(min_ind);cuty(min_ind)];
[~,cut_point_group2]=min(vecnorm(return_points_group2-cut_point_group1));
cut_point_group2=return_points_group2(:,cut_point_group2);
end

longitudinal_vector=(cut_point_group1-cut_point_group2).*1.2;
othorgonal_vector=[longitudinal_vector(2);longitudinal_vector(1).*(-1)];
othorgonal_vector=othorgonal_vector./norm(othorgonal_vector);
middle_point=(cut_point_group1+cut_point_group2)./2;
local_opening_gab=calc_local_opening_gab(parameterized_mesh,middle_point,othorgonal_vector,opening_gap,10);

othorgonal_vector=othorgonal_vector./(norm(othorgonal_vector)).*local_opening_gab;
opening_cut =[middle_point+longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2+othorgonal_vector/2];

% %create a circular opening cut
% opening_cut=[sin(0:(2*pi)/(20-1):2*pi); cos(0:(2*pi)/(20-1):2*pi)].*(local_opening_gab_a/2);
% opening_cut=opening_cut+(cut_point_group1+cut_point_group2)./2;

%circshift the inner group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_1.uv-cut_point_group1));
group_1.uv=circshift(group_1.uv,min_ind*(-1),2);

%circshift the outer group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_2.uv-cut_point_group2));
group_2.uv=circshift(group_2.uv,min_ind*(-1),2);

%fuse both groups
%Check which fusing order is better:
track_combilength1=[group_1.uv group_2.uv];
track_combilength2=[group_2.uv group_1.uv];
track_combilength1=sum(vecnorm(track_combilength1(:,2:end)-track_combilength1(:,1:end-1)));
track_combilength2=sum(vecnorm(track_combilength2(:,2:end)-track_combilength2(:,1:end-1)));
if track_combilength1<track_combilength2
fused_group=[group_1.uv group_2.uv];
else
fused_group=[group_2.uv group_1.uv];
end

%Clear the points witihn the opening_cut to improve to connection
[in_cut_rect_ind,~] = inpolygon(fused_group(1,:),fused_group(2,:),opening_cut(1,:),opening_cut(2,:)); 
fused_group(:,in_cut_rect_ind)=[];

end





% %for the top layer check all possible sequences of group interconnections
% %and chose the shortest one
% if ~isempty(level_positions{current_level}) || enforce_iterative_routing_flag
% routing_plan_exists=false;
% else
% routing_plan_exists=true;
% if numel(groups_to_connect)>1
% %use bruteforce: calculate all possible routes and choose the shortest one
% %which does not cut the boundary
% all_routes=perms(1:numel(groups_to_connect));
% routes_distances=zeros(1,size(all_routes,1));
% for route_ind=1:size(all_routes,1)
% first_groups=all_routes(route_ind,1:end-1);
% second_groups=all_routes(route_ind,2:end);
% route_sub_inds= sub2ind(size(min_group_dists),first_groups,second_groups);
% routes_distances(route_ind)=sum(min_group_dists(route_sub_inds));
% end
% [~,best_route_ind]=min(routes_distances);
% routing_plan=all_routes(best_route_ind,:);
% end
% end

    
% routing_plan_exists=0;
% %for the top level use the routing plan
% if routing_plan_exists
% couple_group1=routing_plan(connect_ind);
% couple_group2=routing_plan(connect_ind+1);
% else %use the next closest group
% %select the pair of groups with shortest respective distance
% min_dist_couple=min_group_dists==min(min_group_dists(:));
% min_dist_couple=find(min_dist_couple);
% [couple_group1,couple_group2]=ind2sub(size(min_group_dists),min_dist_couple(1));
% end
