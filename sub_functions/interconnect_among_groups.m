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

%gather all return points to dismiss them for search of the group cut points
all_cut_points=[];
for group_ind=1:numel(coil_parts(part_ind).connected_group)
all_cut_points=[all_cut_points coil_parts(part_ind).connected_group(group_ind).return_path.uv];
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
grouptracks_to_connect=coil_parts(part_ind).connected_group(groups_to_connect);

%remove the return_path for the search of mutual group cuts
grouptracks_to_connect_without_returns=coil_parts(part_ind).connected_group(groups_to_connect);
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

%find the mininmal distance positions between the groups
%and the points with mininmal distance
min_group_dists=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_group_inds1=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_group_inds2=zeros(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group1=cell(numel(groups_to_connect),numel(groups_to_connect));
min_pos_group2=cell(numel(groups_to_connect),numel(groups_to_connect));
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
opend_group_1=open_group(grouptracks_to_connect(couple_group1),cut_shape_1);
opend_group_2=open_group(grouptracks_to_connect(couple_group2),cut_shape_2);

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
coil_parts(part_ind).connected_group(groups_to_connect(couple_group1));
coil_parts(part_ind).connected_group(groups_to_connect(couple_group1)).uv=fused_group;
[coil_parts(part_ind).connected_group(groups_to_connect(couple_group1)).v,coil_parts(part_ind).connected_group(groups_to_connect(couple_group1)).uv]=uv_to_xyz(fused_group,planary_mesh,curved_mesh);
%upsampled_groups(groups_to_connect(couple_group2))=[];
is_enclosing(couple_group2)=[];
groups_to_connect(couple_group2)=[];
else
coil_parts(part_ind).connected_group(groups_to_connect(couple_group2)).uv=fused_group;
[coil_parts(part_ind).connected_group(groups_to_connect(couple_group2)).v,coil_parts(part_ind).connected_group(groups_to_connect(couple_group2)).uv]=uv_to_xyz(fused_group,planary_mesh,curved_mesh);
%upsampled_groups(groups_to_connect(couple_group1))=[];
is_enclosing(couple_group1)=[];
groups_to_connect(couple_group1)=[];
end

end

end


end

end

%Select the full track as the final return
[~,is_final_ind]=max(arrayfun(@(x) size(coil_parts(part_ind).connected_group(x).v,2), 1:numel(coil_parts(part_ind).connected_group)));
full_track=coil_parts(part_ind).connected_group(is_final_ind(1));
%Shift the open ends to the boundarie of the coil
[~,min_ind]=max(vecnorm(full_track.uv-mean(full_track.uv,2)));
full_track.v=circshift(full_track.v,(-1)*min_ind,2);
full_track.uv=circshift(full_track.uv,(-1)*min_ind,2);

%Assign the outputs
coil_parts(part_ind).wire_path=full_track;

end



end



