function coil_parts= generate_cylindrical_pcb_print_(coil_parts,input)
%Generate a 2D pattern that can be rolled around a cylinder, @Philipp Amrein, Uniklinik
%Freiburg 2022

angular_shrink_factor=1;
pcb_track_width=0.005;


if input.surface_is_cylinder_flag & input.make_cylndrical_pcb

for part_ind=1:numel(coil_parts)

if strcmp(input.pcb_interconnection_method,'spiral_in_out') % first layer spiral in; second layer spiral out

%Rotate the wire on the cylinder 
rot_mat=calc_3d_rotation_matrix_by_vector([input.cylinder_mesh_parameter_list(5) input.cylinder_mesh_parameter_list(6) input.cylinder_mesh_parameter_list(7)]',input.cylinder_mesh_parameter_list(8));
aligned_wire_path=(rot_mat*coil_parts(part_ind).wire_path.v);

%for each point calculate the phi angle
phi_coord=atan2(aligned_wire_path(2,:),aligned_wire_path(1,:));
layout_2d=[phi_coord; aligned_wire_path(3,:)];

%open the wire on the end of angular wrap-arround;
%mark the segments of the wire which have angular wraparounds
positive_wrap=find(diff(layout_2d(1,:))>1.75*pi);
negative_wrap=find(diff(layout_2d(1,:))<(1.75*pi)*(-1));
positive_wrap(positive_wrap==size(layout_2d,2))=[]; %take care of indice overflow problems
negative_wrap(positive_wrap==size(layout_2d,2))=[];
full_wrap_spart_inds=sort([1 size(layout_2d,2) positive_wrap+1 negative_wrap+1] );

%Shift the following points after a phase wrap to the other side in order
%to obtain the clean cuts
layout_2d(:,positive_wrap+1)=layout_2d(:,positive_wrap+1)-[ones(1,numel(positive_wrap)).*2.*pi; zeros(1,numel(positive_wrap))];
layout_2d(:,negative_wrap+1)=layout_2d(:,negative_wrap+1)+[ones(1,numel(negative_wrap)).*2.*pi;  zeros(1,numel(positive_wrap))];

%Generate the wire segmential parts between the wraps
wire_part(numel(full_wrap_spart_inds)-1).uv=[];
wire_part(numel(full_wrap_spart_inds)-1).ind1=[];
wire_part(numel(full_wrap_spart_inds)-1).ind2=[];


for point_ind=1:numel(full_wrap_spart_inds)-1
wire_part(point_ind).uv=layout_2d(:,full_wrap_spart_inds(point_ind)+1:full_wrap_spart_inds(point_ind+1));
wire_part(point_ind).ind1=full_wrap_spart_inds(point_ind)+1;
wire_part(point_ind).ind2=full_wrap_spart_inds(point_ind+1);
end

%Calculate clean cuts on pi,-pi wraps
cut_rectangle=[-pi*angular_shrink_factor pi*angular_shrink_factor pi*angular_shrink_factor -pi*angular_shrink_factor; ...
                                        max(aligned_wire_path(3,:))+abs(max(aligned_wire_path(3,:))).*0.1 ...
                                        max(aligned_wire_path(3,:))+abs(max(aligned_wire_path(3,:))).*0.1 ...
                                        min(aligned_wire_path(3,:))-abs(max(aligned_wire_path(3,:))).*0.1 ...
                                        min(aligned_wire_path(3,:))-abs(max(aligned_wire_path(3,:))).*0.1];
cut_rectangle=[cut_rectangle cut_rectangle(:,1)];


%Add the cut points for a clean cut
for wrap_ind=1:numel(wire_part)
intersection_cut=find_segment_intersections(wire_part(wrap_ind).uv,cut_rectangle);
is_real_cut_ind=find(~isnan([intersection_cut(:).segment_inds]));
if ~isempty(is_real_cut_ind)
wire_part_points=wire_part(wrap_ind).uv;
uv_point=intersection_cut(is_real_cut_ind(1)).uv;
cut_segment_ind=intersection_cut(is_real_cut_ind(1)).segment_inds;
if cut_segment_ind~=1
wire_part(wrap_ind).uv=[wire_part_points(:,1:cut_segment_ind) uv_point wire_part_points(:,cut_segment_ind+1:end-1)];
end
end
end
%add a clean cut also for the second open end of the wire_parts
for wrap_ind=2:numel(wire_part)
if wire_part(wrap_ind-1).uv(1,end)>0
wire_part(wrap_ind).uv=[wire_part(wrap_ind-1).uv(:,end)-[2*pi 0]' wire_part(wrap_ind).uv];
else
wire_part(wrap_ind).uv=[wire_part(wrap_ind-1).uv(:,end)+[2*pi 0]' wire_part(wrap_ind).uv];
end
end


%Generate the track shapes for the indivial wire parts




else % first windings; second layer return paths


end




end





end

%rotate the alinged wire path around the z-axis in order to minimize the
%number of cuts which are neccessary around the angular wrap-arround  



end




% figure;
% hold on;
% %plot(layout_2d(1,:),layout_2d(2,:),'b');
% plot(cut_rectangle(1,:),cut_rectangle(2,:),'r');
% %plot(layout_2d(1,:),layout_2d(2,:),'b');
% scatter(layout_2d(1,positive_wrap),layout_2d(2,positive_wrap),[],[1 0 0],'filled');
% scatter(layout_2d(1,negative_wrap),layout_2d(2,negative_wrap),[],[0 1 0],'filled');
% for wire_part_ind=1:numel(wire_part)
% plot(wire_part(wire_part_ind).uv(1,:),wire_part(wire_part_ind).uv(2,:),'linewidth',2)
% end
% hold off