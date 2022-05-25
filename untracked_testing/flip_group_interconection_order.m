function group_container=flip_group_interconection_order(group_container,group_orientation)
for group_ind=1:numel(group_container)
%flip the order of loops within the groups if indicated by group_ordering
if group_orientation(group_ind)==-1
v_cell=fliplr({group_container(group_ind).loops.v});
uv_cell=fliplr({group_container(group_ind).loops.uv});
point_num_cell=fliplr({group_container(group_ind).loops.number_points});
potential_cell=fliplr({group_container(group_ind).loops.potential});
for vec_ind=1:numel(group_container(group_ind).loops)
group_container(group_ind).loops(vec_ind).v=v_cell{vec_ind};
group_container(group_ind).loops(vec_ind).uv=uv_cell{vec_ind};
group_container(group_ind).loops(vec_ind).number_points=point_num_cell{vec_ind};
group_container(group_ind).loops(vec_ind).potential=potential_cell{vec_ind};
end
end
end
end
