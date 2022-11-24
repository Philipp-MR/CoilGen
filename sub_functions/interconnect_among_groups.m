function coil_parts=interconnect_among_groups(coil_parts,input)
% interconect among groups genereting a single wire track


opening_gap=input.interconnection_cut_width;

%local parameters
cut_height_ratio=1/10;
additional_points_to_remove=2;
%ratio_criteria_for_return_path_usage=1.5;

coil_parts(numel(coil_parts)).opening_cuts_among_groups= []; %save the opening cuts for plotting
coil_parts(numel(coil_parts)).wire_path=[];

for part_ind=1:numel(coil_parts)

connected_group_buff=coil_parts(part_ind).connected_group;

%gather all return points to dismiss them for search of the group cut points
all_cut_points=[];
for group_ind=1:numel(connected_group_buff)
all_cut_points=[all_cut_points connected_group_buff(group_ind).return_path.uv];
end


planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.v);

% winding_distance = calc_avg_winding_distance(coil_parts(part_ind).groups);
% winding_distance(isnan(winding_distance))=0;


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

if num_connections_to_do>0

coil_parts(part_ind).opening_cuts_among_groups(num_connections_to_do).cut1=[];
coil_parts(part_ind).opening_cuts_among_groups(num_connections_to_do).cut2=[];

for connect_ind=1:num_connections_to_do


%get the tracks to connect
grouptracks_to_connect=connected_group_buff(groups_to_connect);

%remove the return_path for the search of mutual group cuts
grouptracks_to_connect_without_returns=connected_group_buff(groups_to_connect);
for ind1=1:numel(groups_to_connect)
for ind2=1:numel(groups_to_connect)
if ind2<ind1
[grouptracks_to_connect_without_returns(ind1).uv,grouptracks_to_connect_without_returns(ind1).v]=remove_points_from_loop(grouptracks_to_connect(ind1),all_cut_points,additional_points_to_remove);
[grouptracks_to_connect_without_returns(ind2).uv,grouptracks_to_connect_without_returns(ind2).v]=remove_points_from_loop(grouptracks_to_connect(ind2),all_cut_points,additional_points_to_remove);
end
end
end

%select the return paths of those interconnected groups for later
%groups_return_paths=return_paths(groups_to_connect);
min_group_dists=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_group_inds1=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_group_inds2=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group1=cell(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group2=cell(numel(groups_to_connect),numel(groups_to_connect));
grouptracks_to_connect_without_returns(numel(groups_to_connect)).unrolled_coords=[];

if input.surface_is_cylinder_flag %in case of cylinder, find the cut points within the "unrolled domain" for getting rid of overlapps with the boardes of printed ciruit boards


%Calculate the "unrolled" 2D coordinates for the cylinder
for group_ind=1:numel(groups_to_connect)
grouptracks_to_connect_without_returns(group_ind).unrolled_coords=[atan2(grouptracks_to_connect_without_returns(group_ind).v(2,:),grouptracks_to_connect_without_returns(group_ind).v(1,:)); grouptracks_to_connect_without_returns(group_ind).v(3,:)];
end
%Find the mutual nearest positions for the "unrolled" groups
min_group_dists=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_group_inds=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group=cell(numel(groups_to_connect),numel(groups_to_connect));
for group_ind_1=1:numel(groups_to_connect)
for group_ind_2=1:numel(groups_to_connect)
if group_ind_1~=group_ind_2
group_a=grouptracks_to_connect_without_returns(group_ind_1).unrolled_coords;
group_b=grouptracks_to_connect_without_returns(group_ind_2).unrolled_coords;
near_ind=zeros(1,size(group_a,2));
near_dist=zeros(1,size(group_a,2));
for point_ind=1:size(group_a,2)
[near_dist(point_ind),near_ind(point_ind)]=min(((group_b(1,:)-group_a(1,point_ind)).^2+(group_b(2,:)-group_a(2,point_ind)).^2).^(1/2));
end
[total_min_dist,total_min_ind]=min(near_dist);
min_group_inds(group_ind_1,group_ind_2)=total_min_ind;
min_pos_group{group_ind_1,group_ind_2}=grouptracks_to_connect_without_returns(group_ind_1).uv(:,total_min_ind);
min_group_dists(group_ind_1,group_ind_2)=total_min_dist;
else
min_pos_group{group_ind_1,group_ind_2}=[nan; nan];
min_group_dists(group_ind_1,group_ind_2)=inf;
end
end
end
min_pos_group1=cellfun(@(x,y) [x y], min_pos_group, min_pos_group', 'UniformOutput', false);
min_pos_group2=cellfun(@(x,y) [x y], min_pos_group', min_pos_group, 'UniformOutput', false);
%select the pair of groups with shortest respective distance
min_dist_couple=min_group_dists==min(min_group_dists(:));
min_dist_couple=find(min_dist_couple);
[couple_group1,couple_group2]=ind2sub(size(min_group_dists),min_dist_couple(1));
%build opening cutshapes for both groups
cut_shape_1=build_cut_circle(min_pos_group1{couple_group1,couple_group2}(:,1),opening_gap);
cut_shape_2=build_cut_circle(min_pos_group2{couple_group1,couple_group2}(:,1),opening_gap);
coil_parts(part_ind).opening_cuts_among_groups(connect_ind).cut1=cut_shape_1;
coil_parts(part_ind).opening_cuts_among_groups(connect_ind).cut2=cut_shape_2;
%Open both groups
opend_group_1=open_group(grouptracks_to_connect(couple_group1),cut_shape_1,cut_shape_1);
opend_group_2=open_group(grouptracks_to_connect(couple_group2),cut_shape_2,cut_shape_2);

else
%find the mininmal distance positions between the groups
%and the points with mininmal distance
for ind1=1:numel(groups_to_connect)
for ind2=1:numel(groups_to_connect)
if ind2<ind1
%remove the return_path for the search of mutual group cuts
[min_group_dists(ind1,ind2),near_points_first,min_group_inds1(ind1,ind2),near_points_second,min_group_inds2(ind1,ind2)]=find_min_mutual_loop_distance(grouptracks_to_connect_without_returns(ind1),grouptracks_to_connect_without_returns(ind2),false);
%[min_group_dists(ind1,ind2),near_points_first,min_group_inds1(ind1,ind2),near_points_second,min_group_inds2(ind1,ind2)]=find_min_mutual_loop_distance(grouptracks_to_connect(ind1),grouptracks_to_connect(ind2));
min_pos_group1{ind1,ind2}=near_points_first.uv;
min_pos_group2{ind1,ind2}=near_points_second.uv;
end
end
end
%min_group_dists=min_group_dists+min_group_dists';
min_group_dists(min_group_dists==0)=Inf;
min_pos_group1=cellfun(@(x,y) [x y], min_pos_group1, min_pos_group1', 'UniformOutput', false);
min_pos_group2=cellfun(@(x,y) [x y], min_pos_group2, min_pos_group2', 'UniformOutput', false);

%select the pair of groups with shortest respective distance
min_dist_couple=min_group_dists==min(min_group_dists(:));
min_dist_couple=find(min_dist_couple);
[couple_group1,couple_group2]=ind2sub(size(min_group_dists),min_dist_couple(1));

% %create and prepare the opening points from the shortest distance points
% raw_cut_point_1=min_pos_group1{couple_group1,couple_group2};
% raw_cut_point_2=min_pos_group2{couple_group1,couple_group2};
% raw_cut_point_dist=min(min_group_dists(:));

%build opening cutshapes for both groups
cut_width_1=calc_local_opening_gab(grouptracks_to_connect_without_returns(couple_group1),min_pos_group1{couple_group1,couple_group2},[],opening_gap);
cut_width_2=calc_local_opening_gab(grouptracks_to_connect_without_returns(couple_group2),min_pos_group2{couple_group1,couple_group2},[],opening_gap);


% cut_width_1=calc_local_opening_gab(grouptracks_to_connect(couple_group1),min_group_inds1(couple_group1,couple_group2),opening_gap);
% cut_width_2=calc_local_opening_gab(grouptracks_to_connect(couple_group2),min_group_inds2(couple_group1,couple_group2),opening_gap);
% cut_shape_1=build_cut_rectangle(grouptracks_to_connect(couple_group1),min_pos_group1{couple_group1,couple_group2},min_group_inds1(couple_group1,couple_group2),cut_width_1,cut_height_ratio);
% cut_shape_2=build_cut_rectangle(grouptracks_to_connect(couple_group2),min_pos_group2{couple_group1,couple_group2},min_group_inds2(couple_group1,couple_group2),cut_width_2,cut_height_ratio);
cut_shape_1=build_cut_circle(min_pos_group1{couple_group1,couple_group2},cut_width_1);
cut_shape_2=build_cut_circle(min_pos_group2{couple_group1,couple_group2},cut_width_2);


%save the cutshapes for later plotting
coil_parts(part_ind).opening_cuts_among_groups(connect_ind).cut1=cut_shape_1;
coil_parts(part_ind).opening_cuts_among_groups(connect_ind).cut2=cut_shape_2;


%Open both groups
opend_group_1=open_group(grouptracks_to_connect(couple_group1),cut_shape_1,cut_shape_2);
opend_group_2=open_group(grouptracks_to_connect(couple_group2),cut_shape_2,cut_shape_1);

end


% %check wether group1 is inside group2 to find out wether the outer turn or
% %the inner turn should be opened
% level_positions=coil_parts(part_ind).level_positions;
% level_positions{cellfun(@isempty, level_positions)}=0;
% level_positions=cell2mat(level_positions);
% level_pos_group1=level_positions(arrayfun(@(x) ismember(groups_to_connect(couple_group1),coil_parts(part_ind).group_levels{x}),1:numel(coil_parts(part_ind).group_levels)));
% level_pos_group2=level_positions(arrayfun(@(x) ismember(groups_to_connect(couple_group2),coil_parts(part_ind).group_levels{x}),1:numel(coil_parts(part_ind).group_levels)));
% 
% %Open both groups
% if level_pos_group1<level_pos_group2
% opend_group_1=open_group(grouptracks_to_connect(couple_group1),cut_shape_1,coil_parts(part_ind).groups(groups_to_connect(couple_group1)),'inner');
% opend_group_2=open_group(grouptracks_to_connect(couple_group2),cut_shape_2,coil_parts(part_ind).groups(groups_to_connect(couple_group2)),'outer');
% elseif level_pos_group1==level_pos_group2
% opend_group_1=open_group(grouptracks_to_connect(couple_group1),cut_shape_1,coil_parts(part_ind).groups(groups_to_connect(couple_group1)),'outer');
% opend_group_2=open_group(grouptracks_to_connect(couple_group2),cut_shape_2,coil_parts(part_ind).groups(groups_to_connect(couple_group2)),'outer');
% elseif level_pos_group1>level_pos_group2   
% opend_group_1=open_group(grouptracks_to_connect(couple_group1),cut_shape_1,coil_parts(part_ind).groups(groups_to_connect(couple_group1)),'outer');
% opend_group_2=open_group(grouptracks_to_connect(couple_group2),cut_shape_2,coil_parts(part_ind).groups(groups_to_connect(couple_group2)),'inner');
% end



%fuse both groups
%Check which fusing order is better:

% % if strcmp(input.group_interconnection_method,'crossed')

track_combilength1=[opend_group_1 opend_group_2];
track_combilength2=[opend_group_2 opend_group_1];
track_combilength1=sum(vecnorm(track_combilength1(:,2:end)-track_combilength1(:,1:end-1)));
track_combilength2=sum(vecnorm(track_combilength2(:,2:end)-track_combilength2(:,1:end-1)));
if track_combilength1<track_combilength2
fused_group=[opend_group_1 opend_group_2];
else
fused_group=[opend_group_2 opend_group_1];
end

% % else % do a strainght interconnection from the middle of the respective opening cuts
% % 
% % end


% %equlize the track
% fused_group = upsample_loop_uv(fused_group,1);


%overwrite fused track into both group point arrays
%delete in all the data containers one of the connected groups to avoid
%redundancy
%do not select here the host level group!
if is_enclosing(couple_group1)
connected_group_buff(groups_to_connect(couple_group1));
connected_group_buff(groups_to_connect(couple_group1)).uv=fused_group;
[connected_group_buff(groups_to_connect(couple_group1)).v,connected_group_buff(groups_to_connect(couple_group1)).uv]=uv_to_xyz(fused_group,planary_mesh,curved_mesh);
%upsampled_groups(groups_to_connect(couple_group2))=[];
is_enclosing(couple_group2)=[];
groups_to_connect(couple_group2)=[];
else
connected_group_buff(groups_to_connect(couple_group2)).uv=fused_group;
[connected_group_buff(groups_to_connect(couple_group2)).v,connected_group_buff(groups_to_connect(couple_group2)).uv]=uv_to_xyz(fused_group,planary_mesh,curved_mesh);
%upsampled_groups(groups_to_connect(couple_group1))=[];
is_enclosing(couple_group1)=[];
groups_to_connect(couple_group1)=[];
end

end

end


end

end

%Select the full track as the final return
[~,is_final_ind]=max(arrayfun(@(x) size(connected_group_buff(x).v,2), 1:numel(connected_group_buff)));
full_track=connected_group_buff(is_final_ind(1));
%Shift the open ends to the boundarie of the coil
[~,min_ind]=max(vecnorm(full_track.uv-mean(full_track.uv,2)));
full_track.v=circshift(full_track.v,(-1)*min_ind,2);
full_track.uv=circshift(full_track.uv,(-1)*min_ind,2);

%Assign the outputs
coil_parts(part_ind).wire_path=full_track;

end



end



