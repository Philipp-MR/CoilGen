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