function loop_normal=calc_mean_loop_normal(group,coil_mesh)


all_loop_normals=zeros(3,numel(group.loops));
for loop_ind=1:numel(group.loops)
group_center=mean(group.loops(loop_ind).v,2);
loop_vecs=group.loops(loop_ind).v(:,2:end)-group.loops(loop_ind).v(:,1:end-1);
center_vecs=group.loops(loop_ind).v(:,1:end-1)-group_center;
all_loop_normals(:,loop_ind)=mean(cross(loop_vecs,center_vecs),2);
end

loop_normal=mean(all_loop_normals,2);
loop_normal=loop_normal./norm(loop_normal);

%make sure the loop normal points outwards seen from the coordiante center
if sum(loop_normal.*group_center)<0
loop_normal=loop_normal.*(-1);
end
%loop_normal=loop_normal.*sign(group.loops(1).potential);



% all_loops=[];
% for contour_ind=1:numel(group.loops)
% all_loops=[all_loops group.loops(contour_ind).uv];
% end
% curved_mesh=triangulation(coil_mesh.faces',coil_mesh.vertices');
% %Calculate vertex normals for later
% face_normal=faceNormal(curved_mesh);
% %make sure that they are pointing to the outisde of the surface
% if mean(dot([vertexNormal(curved_mesh)]',[curved_mesh.Points-mean(curved_mesh.Points)]'))<0
% face_normal=face_normal.*(-1);
% end
% [target_triangle_normal,~] = pointLocation(triangulation(coil_mesh.faces',coil_mesh.uv'),all_loops(1,:)',all_loops(2,:)');
% loop_normal=mean(face_normal(target_triangle_normal(~isnan(target_triangle_normal)),:),1)';
% loop_normal=loop_normal./norm(loop_normal);

end
