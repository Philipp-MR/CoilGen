function coil_parts=interconnect_among_groups(coil_parts,input)
%Interconect among groups genereting a single wire track



%local parameters
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
for group_ind=1:numel(groups_to_connect)
[grouptracks_to_connect_without_returns(group_ind).uv,grouptracks_to_connect_without_returns(group_ind).v]=remove_points_from_loop(grouptracks_to_connect(group_ind),all_cut_points,additional_points_to_remove);
grouptracks_to_connect_without_returns(group_ind).unrolled_coords=[atan2(grouptracks_to_connect_without_returns(group_ind).v(2,:),grouptracks_to_connect_without_returns(group_ind).v(1,:)); grouptracks_to_connect_without_returns(group_ind).v(3,:)];
end

%select the return paths of those interconnected groups for later
%groups_return_paths=return_paths(groups_to_connect);
min_group_dists=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_group_inds=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group.v=cell(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group.uv=cell(numel(groups_to_connect),numel(groups_to_connect));


%find the mininmal distance positions between the groups
%and the points with mininmal distance
for ind1=1:numel(groups_to_connect)
for ind2=1:numel(groups_to_connect)
if ind2~=ind1
% if input.surface_is_cylinder_flag %in case of cylinder, find the cut points within the "unrolled domain" for getting rid of overlapps with the boardes of printed ciruit boards
    % group_a=grouptracks_to_connect_without_returns(ind1).unrolled_coords;
    % group_b=grouptracks_to_connect_without_returns(ind2).unrolled_coords;
    group_a=grouptracks_to_connect_without_returns(ind1);
    group_b=grouptracks_to_connect_without_returns(ind2);
    near_ind=zeros(1,size(group_a.v,2));
    near_dist=zeros(1,size(group_a.v,2));
    for point_ind=1:size(group_a.v,2) % use "for" loop to avoid memory overflow
    %[near_dist(point_ind),near_ind(point_ind)]=min(((group_b.v(1,:)-group_a.v(1,point_ind)).^2+(group_b.v(2,:)-group_a.v(2,point_ind)).^2).^(1/2));
    [near_dist(point_ind),near_ind(point_ind)]=min(((group_b.v(1,:)-group_a.v(1,point_ind)).^2+(group_b.v(2,:)-group_a.v(2,point_ind)).^2+(group_b.v(3,:)-group_a.v(3,point_ind)).^2).^(1/2));
    end
    [total_min_dist,total_min_ind]=min(near_dist);
    min_group_inds(ind1,ind2)=total_min_ind;
   
    min_pos_group.uv{ind1,ind2}=grouptracks_to_connect_without_returns(ind1).uv(:,total_min_ind);
    min_pos_group.v{ind1,ind2}=grouptracks_to_connect_without_returns(ind1).v(:,total_min_ind);
    min_group_dists(ind1,ind2)=total_min_dist;
    min_group_dists(min_group_dists==0)=Inf;
% else
%     %remove the return_path for the search of mutual group cuts
%     [min_group_dists(ind1,ind2),near_points_first,~,~,~]=find_min_mutual_loop_distance(group_a,group_b,false);
%     min_pos_group.uv{ind1,ind2}=near_points_first.uv;
% end
end
end
end



%Select the pair of groups with shortest respective distance
min_dist_couple=min_group_dists==min(min_group_dists(:));
min_dist_couple=find(min_dist_couple);
[couple_group1,couple_group2]=ind2sub(size(min_group_dists),min_dist_couple(1));

%Open the loop
[opened_group_1,cut_shape_1,~]=open_loop_with_3d_sphere(grouptracks_to_connect(couple_group1),min_pos_group.v{couple_group1,couple_group2}(:,1),input.interconnection_cut_width);
[opened_group_2,cut_shape_2,~]=open_loop_with_3d_sphere(grouptracks_to_connect(couple_group2),min_pos_group.v{couple_group2,couple_group1}(:,1),input.interconnection_cut_width);

%Save the cutshapes for later plotting
coil_parts(part_ind).opening_cuts_among_groups(connect_ind).cut1=cut_shape_1;
coil_parts(part_ind).opening_cuts_among_groups(connect_ind).cut2=cut_shape_2;




%fuse both groups
%Check which fusing order is better:

% % if strcmp(input.group_interconnection_method,'crossed')

track_combilength1=[opened_group_1.v opened_group_2.v];
track_combilength2=[opened_group_2.v opened_group_1.v];
track_combilength1=sum(vecnorm(track_combilength1(:,2:end)-track_combilength1(:,1:end-1)));
track_combilength2=sum(vecnorm(track_combilength2(:,2:end)-track_combilength2(:,1:end-1)));
if track_combilength1<track_combilength2
fused_group.v=[opened_group_1.v opened_group_2.v];
fused_group.uv=[opened_group_1.uv opened_group_2.uv];
else
fused_group.v=[opened_group_2.v opened_group_1.v];
fused_group.uv=[opened_group_2.uv opened_group_1.uv];
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
connected_group_buff(groups_to_connect(couple_group1)).uv=fused_group.uv;
connected_group_buff(groups_to_connect(couple_group1)).v=fused_group.v;
is_enclosing(couple_group2)=[];
groups_to_connect(couple_group2)=[];
else
connected_group_buff(groups_to_connect(couple_group2)).uv=fused_group.uv;
connected_group_buff(groups_to_connect(couple_group2)).v=fused_group.v;
is_enclosing(couple_group1)=[];
groups_to_connect(couple_group1)=[];
end

end
end
end
end

%Select the full track as the final return
[~,is_final_ind]=max(arrayfun(@(x) size(connected_group_buff(x).v,2), 1:numel(connected_group_buff))); % select the biggest part..
full_track=connected_group_buff(is_final_ind(1));
%Shift the open ends to the boundarie of the coil
[~,min_ind]=max(full_track.v(3,:));
full_track.v=circshift(full_track.v,(-1)*min_ind,2);
full_track.uv=circshift(full_track.uv,(-1)*min_ind,2);

%Assign the outputs
coil_parts(part_ind).wire_path=full_track;

end

end