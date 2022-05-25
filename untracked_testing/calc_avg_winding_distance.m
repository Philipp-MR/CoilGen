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