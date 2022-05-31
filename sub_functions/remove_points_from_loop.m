function [loop_out_uv,loop_out_v]= remove_points_from_loop(loop,points_to_remove,boundary_threshold)
%remove points with idenditcal uv coordinates from a loop, even with some
%additional more points around 

rep_u=repmat(loop.uv(1,:),[size(points_to_remove,2) 1]);
rep_v=repmat(loop.uv(2,:),[size(points_to_remove,2) 1]);

rep_u2=repmat(points_to_remove(1,:),[size(loop.uv,2) 1]);
rep_v2=repmat(points_to_remove(2,:),[size(loop.uv,2) 1]);


identical_point_inds=find((rep_u==rep_u2') & (rep_v==rep_v2')); 

[~,identical_point_inds] = ind2sub(size(rep_u),identical_point_inds);

below_inds=(min(identical_point_inds)-boundary_threshold):min(identical_point_inds);
below_inds(below_inds<0)=below_inds(below_inds<0)+size(loop.uv,2);

abow_inds=max(identical_point_inds):(max(identical_point_inds)+boundary_threshold);
abow_inds(abow_inds>size(loop.uv,2))=abow_inds(abow_inds>size(loop.uv,2))-size(loop.uv,2);

%add more points as a "boundary thresshold"
full_point_inds_to_remove=[below_inds identical_point_inds' abow_inds];

[inds_out,~] = setdiff(1:size(loop.uv,2),full_point_inds_to_remove,'stable');

loop_out_uv=loop.uv(:,inds_out);
loop_out_v=loop.v(:,inds_out);


end

