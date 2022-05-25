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