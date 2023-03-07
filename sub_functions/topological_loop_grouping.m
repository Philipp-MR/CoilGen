function coil_parts = topological_loop_grouping(coil_parts,input)
%Group the contour loops in topological order
%check for all loop enclosements of other loops



coil_parts(numel(coil_parts)).groups=[];
coil_parts(numel(coil_parts)).group_levels=[];
coil_parts(numel(coil_parts)).level_positions=[];


for part_ind=1:numel(coil_parts)

%create cell array for each loop that give the embodied loops
num_total_loops=numel(coil_parts(part_ind).contour_lines);
loop_in_loop_mat=zeros(num_total_loops);
for loop_to_test=1:num_total_loops
for loop_num=1:num_total_loops
if loop_to_test~=loop_num
loop_in_loop_mat(loop_to_test,loop_num)=check_mutual_loop_inclusion(coil_parts(part_ind).contour_lines(loop_num).uv,coil_parts(part_ind).contour_lines(loop_to_test).uv);
end
end
end
loop_in_loop_mat=loop_in_loop_mat.*~diag(1:num_total_loops);
lower_loops=cell(1,num_total_loops);
for loop_to_test=1:num_total_loops
lower_loops{loop_to_test}=find(loop_in_loop_mat(:,loop_to_test));
end
higher_loops=cell(1,num_total_loops);
for loop_to_test=1:num_total_loops
higher_loops{loop_to_test}=find(loop_in_loop_mat(loop_to_test,:));
end
%assign '0' to loops that have no lower loops
empty_cells=find(cellfun(@isempty,lower_loops));
for ssss=1:numel(empty_cells)
lower_loops(empty_cells(ssss))={0};
end

%convert the list of lower loops into parallel levels
coil_parts(part_ind).group_levels=cell(1,num_total_loops);
for aaaa=1:num_total_loops
coil_parts(part_ind).group_levels{aaaa}=find(cellfun(@(x) isequal(lower_loops{aaaa},x),lower_loops));
end
%check the posibillity of top level is composed out of a single group
%that means that there is a loop that contains all other groups
is_global_top_loop=arrayfun(@(x) numel(higher_loops{x})==num_total_loops-1,1:numel(higher_loops));
%delete the repetition in the parallel levels and the singular levels
coil_parts(part_ind).group_levels=cellfun(@str2num,unique(cellfun(@num2str,coil_parts(part_ind).group_levels,'UniformOutput',0)),'UniformOutput',0); %convert level indice to char; check for uniqunes; convert back to num
coil_parts(part_ind).group_levels=coil_parts(part_ind).group_levels(cellfun(@numel,coil_parts(part_ind).group_levels)~=1 | arrayfun(@(x) is_global_top_loop(coil_parts(part_ind).group_levels{x}(1)),1:numel(coil_parts(part_ind).group_levels))==1); % remove levels with only one member except it is the singular top level
%creating the loop groups (containing still the loops of the inner groups)
overlapping_loop_groups_num=horzcat(coil_parts(part_ind).group_levels{:});
overlapping_loop_groups=num2cell(horzcat(coil_parts(part_ind).group_levels{:}));
for pppp=1:numel(overlapping_loop_groups)
        overlapping_loop_groups{pppp}=[overlapping_loop_groups{pppp} higher_loops{overlapping_loop_groups_num(pppp)}];
end
%build the group topology by checking the loop content of a certain group
%to see if it is sup-set of the loop content of another group
group_in_group_mat=zeros(numel(overlapping_loop_groups),numel(overlapping_loop_groups));
for aaaa=1:numel(overlapping_loop_groups)
for bbbb=1:numel(overlapping_loop_groups)
group_in_group_mat(aaaa,bbbb)=all(ismember(overlapping_loop_groups{aaaa}, overlapping_loop_groups{bbbb}));
end
group_in_group_mat(aaaa,aaaa)=0;
end
group_is_subgroup_of=cell(1,numel(overlapping_loop_groups));
for aaaa=1:numel(overlapping_loop_groups)
group_is_subgroup_of{aaaa}=find(group_in_group_mat(aaaa,:)==1);
end
group_contains_following_group=cell(1,numel(overlapping_loop_groups));
for aaaa=1:numel(overlapping_loop_groups)
group_contains_following_group{aaaa}=find(group_in_group_mat(:,aaaa)==1);
end
%remove loops from group if they also belong to a respective subgroup
coil_parts(part_ind).loop_groups=overlapping_loop_groups;
for aaaa=1:numel(overlapping_loop_groups)
loops_to_remove=[];
for bbbb=1:numel(group_contains_following_group{aaaa})
loops_to_remove=[loops_to_remove overlapping_loop_groups{group_contains_following_group{aaaa}(bbbb)}];
end
coil_parts(part_ind).loop_groups{aaaa}=setdiff(overlapping_loop_groups{aaaa},loops_to_remove);
end
%order the groups on the number of loops
renamed_group_levels=coil_parts(part_ind).group_levels;
[~,sort_ind]=sort(cellfun(@numel,coil_parts(part_ind).loop_groups),'descend');
coil_parts(part_ind).loop_groups=coil_parts(part_ind).loop_groups(sort_ind);
%renumerate (=rename) the groups (also in the levels)
for aaaa=1:numel(coil_parts(part_ind).group_levels)
for bbbb=1:numel(coil_parts(part_ind).group_levels{aaaa})
for cccc=1:numel(coil_parts(part_ind).loop_groups)
if ismember(coil_parts(part_ind).group_levels{aaaa}(bbbb),coil_parts(part_ind).loop_groups{cccc})
renamed_group_levels{aaaa}(bbbb)=cccc;
end
end
end
end
%resort parallel_levels to new group names
[~,sort_ind_level]=sort(cellfun(@min,renamed_group_levels));
renamed_group_levels=renamed_group_levels(sort_ind_level);
coil_parts(part_ind).group_levels = cellfun(@sort,renamed_group_levels,'UniformOutput',0);

%find for each parallel level the groups that contain that level (to build the "level_positions")
%first find the loops of unique loops of a certain level
loops_per_level=cell(1,numel(coil_parts(part_ind).group_levels));
for aaaa=1:numel(coil_parts(part_ind).group_levels)
for bbbb=1:numel(coil_parts(part_ind).group_levels{aaaa})
loops_per_level{aaaa}=[loops_per_level{aaaa} coil_parts(part_ind).loop_groups{coil_parts(part_ind).group_levels{aaaa}(bbbb)}];
end
end

level_enclosed_by_loop=cell(1,numel(coil_parts(part_ind).group_levels));
for aaaa=1:numel(coil_parts(part_ind).group_levels)
for bbbb=1:numel(loops_per_level{aaaa})
for cccc=1:numel(coil_parts(part_ind).contour_lines)
if loop_in_loop_mat(cccc,loops_per_level{aaaa}(bbbb))==1
level_enclosed_by_loop{aaaa}=[level_enclosed_by_loop{aaaa} cccc];
end
end
end
end

coil_parts(part_ind).level_positions=cell(1,numel(coil_parts(part_ind).group_levels));
for aaaa=1:numel(level_enclosed_by_loop)
for bbbb=1:numel(coil_parts(part_ind).loop_groups)
if any(ismember(level_enclosed_by_loop{aaaa},coil_parts(part_ind).loop_groups{bbbb}))
coil_parts(part_ind).level_positions{aaaa}=[coil_parts(part_ind).level_positions{aaaa} bbbb];
end
end
end
for aaaa=1:numel(coil_parts(part_ind).level_positions)
coil_parts(part_ind).level_positions{aaaa}=setdiff(coil_parts(part_ind).level_positions{aaaa},coil_parts(part_ind).group_levels{aaaa});
end

%sort the level_positions according to their rank
rank_of_group=zeros(1,numel(coil_parts(part_ind).loop_groups));
for aaaa=1:numel(coil_parts(part_ind).loop_groups)
rank_of_group(aaaa)=numel(coil_parts(part_ind).level_positions{cellfun(@(x) ismember(aaaa,x),coil_parts(part_ind).group_levels)});
end
for aaaa=1:numel(coil_parts(part_ind).group_levels)
    [~,sort_ind_rank]=sort(rank_of_group(coil_parts(part_ind).level_positions{aaaa}));
	coil_parts(part_ind).level_positions{aaaa}=coil_parts(part_ind).level_positions{aaaa}(sort_ind_rank);
end

%build the group container
for iiii=1:numel(coil_parts(part_ind).loop_groups)
%sort the loops in each group according to the rank
[~,sort_ind_loops]=sort(cellfun(@(x) numel(higher_loops{x}),num2cell(coil_parts(part_ind).loop_groups{iiii})),'descend');
kkkk=0;
for jjjj=sort_ind_loops
kkkk=kkkk+1;
coil_parts(part_ind).groups(iiii).loops(kkkk).number_points=size(coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).uv,2);
coil_parts(part_ind).groups(iiii).loops(kkkk).v=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).v;
coil_parts(part_ind).groups(iiii).loops(kkkk).uv=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).uv;
coil_parts(part_ind).groups(iiii).loops(kkkk).potential=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).potential;
coil_parts(part_ind).groups(iiii).loops(kkkk).current_orientation=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).current_orientation;
end
% %sort the loops in each group according to their rank
% [~,sort_ind_loops]=sort(cellfun(@numel,lower_loops(loop_groups{iiii})));
% kkkk=0;
% for jjjj=sort_ind_loops
% kkkk=kkkk+1;
% coil_parts(part_ind).groups(iiii).loops(kkkk).number_points=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).number_points;
% coil_parts(part_ind).groups(iiii).loops(kkkk).v=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).v;
% coil_parts(part_ind).groups(iiii).loops(kkkk).uv=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).uv;
% coil_parts(part_ind).groups(iiii).loops(kkkk).potential=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).potential;
% coil_parts(part_ind).groups(iiii).loops(kkkk).current_orientation=coil_parts(part_ind).contour_lines(coil_parts(part_ind).loop_groups{iiii}(jjjj)).current_orientation;
% end
end


% %Sort the groups wihtin their level accoring to theri average z-position
% if input.sort_groups_along_z
% avg_z_value=zeros(1,numel(coil_parts(part_ind).loop_groups));
% old_group_inds=1:numel(coil_parts(part_ind).groups);
% for group_ind=1:numel(coil_parts(part_ind).groups)
% all_points=[];
% for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
% all_points=[all_points coil_parts(part_ind).groups(group_ind).loops(loop_ind).v];
% end
% avg_z_value(group_ind)=sum(sum(all_points.*[0.05 0 1]',1))./size(all_points,2);
% end
% [~,new_group_inds]=sort(avg_z_value);
% avg_z_value=avg_z_value(new_group_inds);
% avg_z_value=avg_z_value-min(avg_z_value);
% avg_z_value=avg_z_value./max(avg_z_value);
% %Update the level container
% for level_ind=1:numel(coil_parts(part_ind).group_levels)
% coil_parts(part_ind).group_levels{level_ind}=new_group_inds(coil_parts(part_ind).group_levels{level_ind});
% coil_parts(part_ind).level_positions{level_ind}=new_group_inds(coil_parts(part_ind).level_positions{level_ind});
% end
% %Upate the group structs
% coil_parts(part_ind).loop_groups=coil_parts(part_ind).loop_groups(new_group_inds);
% coil_parts(part_ind).groups=coil_parts(part_ind).groups(new_group_inds);

% loop_groups_old=coil_parts(part_ind).loop_groups;
% groups_old=coil_parts(part_ind).groups;
% for group_ind=1:numel(coil_parts(part_ind).groups)
% coil_parts(part_ind).loop_groups(group_ind)=loop_groups_old(new_group_inds(group_ind));
% coil_parts(part_ind).groups(group_ind)=groups_old(new_group_inds(group_ind));
% end
% end

end


% figure;
% color_map= cool(numel(coil_parts(part_ind).groups));
% hold on;
% for group_ind=1:numel(coil_parts(part_ind).groups)
% for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
% plot3(coil_parts.groups(group_ind).loops(loop_ind).v(1,:),coil_parts.groups(group_ind).loops(loop_ind).v(2,:),coil_parts.groups(group_ind).loops(loop_ind).v(3,:),'color',color_map(group_ind,:));
% end
% end
% axis equal;
% hold off;


end
