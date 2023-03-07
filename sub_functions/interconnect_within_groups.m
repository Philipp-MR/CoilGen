function coil_parts=interconnect_within_groups(coil_parts,input)

cut_plane_definition='nearest'; % nearest or B0
cut_height_ratio=1/2; %the ratio of height to widh of the indivuial cut_shapes

coil_parts(numel(coil_parts)).connected_group=[];
coil_parts(numel(coil_parts)).return_paths=[];

for part_ind=1:numel(coil_parts)

planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.v);
coil_parts(numel(coil_parts)).connected_group(numel(coil_parts(part_ind).groups)).uv=[];

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
cut_position(numel(coil_parts(part_ind).groups)).group=[];

%Sort the force cut selection wihtin their level accoring to theri average z-position
avg_z_value=zeros(1,numel(coil_parts(part_ind).loop_groups));
old_group_inds=1:numel(coil_parts(part_ind).groups);
for group_ind=1:numel(coil_parts(part_ind).groups)
all_points=[];
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
all_points=[all_points coil_parts(part_ind).groups(group_ind).loops(loop_ind).v];
end
avg_z_value(group_ind)=sum(sum(all_points.*[0.05 0 1]',1))./size(all_points,2);
end
[~,new_group_inds]=sort(avg_z_value);




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
cut_position(group_ind).group=find_group_cut_position(coil_parts(part_ind).groups(group_ind),coil_parts(part_ind).group_centers.v(:,group_ind),coil_parts(part_ind).coil_mesh,input.b_0_direction,cut_plane_definition);

%choose either low or high cutshape
force_cut_selection=force_cut_selection(new_group_inds);
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
switch force_cut_selection{group_ind}
    case 'high'
cut_position_used=cut_position(group_ind).group(loop_ind).high_cut.v;
    case 'low'
cut_position_used=cut_position(group_ind).group(loop_ind).low_cut.v;
    otherwise
cut_position_used=cut_position(group_ind).group(loop_ind).high_cut.v;
end

%Open the loop
[opened_loop,coil_parts(part_ind).groups(group_ind).cutshape(loop_ind).uv,~]=open_loop_with_3d_sphere(coil_parts(part_ind).groups(group_ind).loops(loop_ind),cut_position_used,input.interconnection_cut_width);
coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv=opened_loop.uv;
coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).v=opened_loop.v;

end

%Build the interconnected group by adding the opened loops
coil_parts(part_ind).connected_group(group_ind).spiral_in.uv=[];
coil_parts(part_ind).connected_group(group_ind).spiral_in.v=[];
coil_parts(part_ind).connected_group(group_ind).spiral_out.uv=[];
coil_parts(part_ind).connected_group(group_ind).spiral_out.v=[];
coil_parts(part_ind).connected_group(group_ind).uv=[];
coil_parts(part_ind).connected_group(group_ind).v=[];
coil_parts(part_ind).connected_group(group_ind).return_path.uv=[];
coil_parts(part_ind).connected_group(group_ind).return_path.v=[];
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops)
coil_parts(part_ind).connected_group(group_ind).uv=[coil_parts(part_ind).connected_group(group_ind).uv coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv];
coil_parts(part_ind).connected_group(group_ind).v=[coil_parts(part_ind).connected_group(group_ind).v coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).v];
coil_parts(part_ind).connected_group(group_ind).spiral_in.uv=[coil_parts(part_ind).connected_group(group_ind).spiral_in.uv coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv];
coil_parts(part_ind).connected_group(group_ind).spiral_in.v=[coil_parts(part_ind).connected_group(group_ind).spiral_in.v coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).v];
coil_parts(part_ind).connected_group(group_ind).spiral_out.uv=[coil_parts(part_ind).connected_group(group_ind).spiral_out.uv coil_parts(part_ind).groups(group_ind).opened_loop(numel(coil_parts(part_ind).groups(group_ind).loops)+1-loop_ind).uv];
coil_parts(part_ind).connected_group(group_ind).spiral_out.v=[coil_parts(part_ind).connected_group(group_ind).spiral_out.v coil_parts(part_ind).groups(group_ind).opened_loop(numel(coil_parts(part_ind).groups(group_ind).loops)+1-loop_ind).v];
end
%Add the return path
for loop_ind=numel(coil_parts(part_ind).groups(group_ind).loops):-1:1
coil_parts(part_ind).connected_group(group_ind).return_path.uv=[coil_parts(part_ind).connected_group(group_ind).return_path.uv mean(coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).uv(:,[1 end]),2)];
coil_parts(part_ind).connected_group(group_ind).return_path.v=[coil_parts(part_ind).connected_group(group_ind).return_path.v mean(coil_parts(part_ind).groups(group_ind).opened_loop(loop_ind).v(:,[1 end]),2)];
end
coil_parts(part_ind).connected_group(group_ind).uv=[coil_parts(part_ind).connected_group(group_ind).uv coil_parts(part_ind).connected_group(group_ind).return_path.uv];
coil_parts(part_ind).connected_group(group_ind).v=[coil_parts(part_ind).connected_group(group_ind).v coil_parts(part_ind).connected_group(group_ind).return_path.v];
%Close the connected group
coil_parts(part_ind).connected_group(group_ind).uv=[coil_parts(part_ind).connected_group(group_ind).uv coil_parts(part_ind).connected_group(group_ind).uv(:,1)];
coil_parts(part_ind).connected_group(group_ind).v=[coil_parts(part_ind).connected_group(group_ind).v coil_parts(part_ind).connected_group(group_ind).v(:,1)];

% %Transform to the curved domain
% [coil_parts(part_ind).connected_group(group_ind).v,~]=uv_to_xyz(coil_parts(part_ind).connected_group(group_ind).uv,planary_mesh,curved_mesh);
% [coil_parts(part_ind).connected_group(group_ind).spiral_in.v,~]=uv_to_xyz(coil_parts(part_ind).connected_group(group_ind).spiral_in.uv,planary_mesh,curved_mesh);
% [coil_parts(part_ind).connected_group(group_ind).spiral_out.v,~]=uv_to_xyz(coil_parts(part_ind).connected_group(group_ind).spiral_out.uv,planary_mesh,curved_mesh);

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