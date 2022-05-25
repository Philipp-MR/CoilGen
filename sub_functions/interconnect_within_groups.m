function coil_parts=interconnect_within_groups(coil_parts,input)

cut_plane_definition='nearest'; % nearest or B0
cut_height_ratio=1/2; %the ratio of height to widh of the indivuial cut_shapes

coil_parts(numel(coil_parts)).connected_group=[];
coil_parts(numel(coil_parts)).return_paths=[];

for part_ind=1:numel(coil_parts)

planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.v);

coil_parts(numel(coil_parts)).connected_group(numel(coil_parts(part_ind).groups)).spiral_in_uv=[];
coil_parts(numel(coil_parts)).connected_group(numel(coil_parts(part_ind).groups)).spiral_out_uv=[];
coil_parts(numel(coil_parts)).connected_group(numel(coil_parts(part_ind).groups)).uv=[];
coil_parts(numel(coil_parts)).connected_group(numel(coil_parts(part_ind).groups)).return_path_uv=[];

%take the cut selection if it is given in the input
switch numel(input.force_cut_selection)
    case 0
force_cut_selection=repmat({'none'},1,numel(coil_parts(part_ind).groups));
    case 1
if strcmp(input.force_cut_selection{1},'high')
force_cut_selection=repmat({'high'},1,numel(coil_parts(part_ind).groups));
else
force_cut_selection=repmat({'none'},1,numel(coil_parts(part_ind).groups));
end
    otherwise
force_cut_selection=input.force_cut_selection;
end


%Generate cutshapes, open and interconnect the loops within each group
cut_shape(numel(coil_parts(part_ind).groups)).group=[];

for group_ind=1:numel(coil_parts(part_ind).groups)

if numel(coil_parts(part_ind).groups(group_ind).loops)==1

%if the groups consists of only one loop; it is not neccessary to open it
coil_parts(part_ind).groups(group_ind).opened_loop(group_ind).uv=coil_parts(part_ind).groups(group_ind).loops.uv;
coil_parts(part_ind).groups(group_ind).opened_loop(group_ind).v=coil_parts(part_ind).groups(group_ind).loops.v;
coil_parts(part_ind).groups(group_ind).cutshape.uv=[nan; nan];
coil_parts(part_ind).connected_group(group_ind).return_path.uv=[nan; nan];
coil_parts(part_ind).connected_group(group_ind).uv=coil_parts(part_ind).groups(group_ind).loops.uv;
coil_parts(part_ind).connected_group(group_ind).v=coil_parts(part_ind).groups(group_ind).loops.v;
coil_parts(part_ind).connected_group(group_ind).spiral_in.uv=coil_parts(part_ind).groups(group_ind).loops.uv;
coil_parts(part_ind).connected_group(group_ind).spiral_in.v=coil_parts(part_ind).groups(group_ind).loops.v;
coil_parts(part_ind).connected_group(group_ind).spiral_out.uv=coil_parts(part_ind).groups(group_ind).loops.uv;
coil_parts(part_ind).connected_group(group_ind).spiral_out.v=coil_parts(part_ind).groups(group_ind).loops.v;

else

coil_parts(part_ind).groups(group_ind).opened_loop(numel(coil_parts(part_ind).groups(group_ind))).uv=[];
coil_parts(part_ind).groups(group_ind).cutshape(numel(coil_parts(part_ind).groups(group_ind))).uv=[];

%generate the cutshapes for all the loops within the group
cut_shape(group_ind).group=generate_group_cutshapes(coil_parts(part_ind).groups(group_ind),coil_parts(part_ind).group_centers.v(:,group_ind),coil_parts(part_ind).coil_mesh,input.interconnection_cut_width,input.b_0_direction,cut_plane_definition,cut_height_ratio);

%choose either low or high cutshape
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
switch force_cut_selection{group_ind}
    case 'high'
cut_shape_points=cut_shape(group_ind).group(loop_ind).high_cutshape;
    case 'low'
cut_shape_points=cut_shape(group_ind).group(loop_ind).low_cutshape;
    otherwise
cut_shape_points=cut_shape(group_ind).group(loop_ind).high_cutshape;
end
coil_parts(part_ind).groups(group_ind).cutshape(loop_ind).uv=cut_shape_points;
coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv=open_loop(coil_parts(part_ind).groups(group_ind).loops(loop_ind),cut_shape_points);
end

%build the interconnected group by adding the opened loops

coil_parts(part_ind).connected_group(group_ind).spiral_in.uv=[];
coil_parts(part_ind).connected_group(group_ind).spiral_out.uv=[];
coil_parts(part_ind).connected_group(group_ind).uv=[];
coil_parts(part_ind).connected_group(group_ind).return_path.uv=[];
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
coil_parts(part_ind).connected_group(group_ind).uv=[coil_parts(part_ind).connected_group(group_ind).uv coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv];
coil_parts(part_ind).connected_group(group_ind).spiral_out.uv=[coil_parts(part_ind).connected_group(group_ind).spiral_out.uv coil_parts(part_ind).groups(group_ind).opened_loop(numel(coil_parts(part_ind).groups(group_ind).loops)+1-loop_ind).uv];
end
coil_parts(part_ind).connected_group(group_ind).spiral_in.uv=coil_parts(part_ind).connected_group(group_ind).spiral_in.uv;
%add the return path
for loop_ind=numel(coil_parts(part_ind).groups(group_ind).loops):-1:1
coil_parts(part_ind).connected_group(group_ind).return_path.uv=[coil_parts(part_ind).connected_group(group_ind).return_path.uv mean(coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv(:,[1 end]),2)];
end
coil_parts(part_ind).connected_group(group_ind).uv=[coil_parts(part_ind).connected_group(group_ind).uv coil_parts(part_ind).connected_group(group_ind).return_path.uv];

%close the connected group
coil_parts(part_ind).connected_group(group_ind).uv=[coil_parts(part_ind).connected_group(group_ind).uv coil_parts(part_ind).connected_group(group_ind).uv(:,1)];

%transform to the curved domain
[coil_parts(part_ind).connected_group(group_ind).v,~]=uv_to_xyz(coil_parts(part_ind).connected_group(group_ind).uv,planary_mesh,curved_mesh);
[coil_parts(part_ind).connected_group(group_ind).spiral_in.v,~]=uv_to_xyz(coil_parts(part_ind).connected_group(group_ind).spiral_in.uv,planary_mesh,curved_mesh);
[coil_parts(part_ind).connected_group(group_ind).spiral_out.v,~]=uv_to_xyz(coil_parts(part_ind).connected_group(group_ind).spiral_out.uv,planary_mesh,curved_mesh);

end

end

end

% %Plot results 2D
% for part_ind=1:numel(coil_parts)
% figure;
% for group_ind=1:numel(coil_parts(part_ind).groups)
% hold on;
% plot(coil_parts(part_ind).connected_group(group_ind).uv(1,:),coil_parts(part_ind).connected_group(group_ind).uv(2,:));
% for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
% plot(coil_parts(part_ind).groups(group_ind).cutshape(loop_ind).uv(1,:),coil_parts(part_ind).groups(group_ind).cutshape(loop_ind).uv(2,:),'r');
% end
% end
% axis equal;
% xlabel('u');ylabel('v'); 
% hold off;
% end
% %Plot results 3D
% figure;
% hold on;
% for part_ind=1:numel(coil_parts)
% for group_ind=1:numel(coil_parts(part_ind).groups)
% plot3(coil_parts(part_ind).connected_group(group_ind).v(1,:),coil_parts(part_ind).connected_group(group_ind).v(2,:),coil_parts(part_ind).connected_group(group_ind).v(3,:));
% end
% end
% axis equal;
% xlabel('x');ylabel('y'); zlabel('z'); 
% hold off;



end