function coil_parts= generate_cylindrical_pcb_print(coil_parts,input)
%Generate a 2D pattern that can be rolled around a cylinder, @Philipp Amrein, Uniklinik
%Freiburg 2022


pcb_track_width=input.conductor_cross_section_width;
cylinder_radius=input.cylinder_mesh_parameter_list(2);
rot_mat=calc_3d_rotation_matrix_by_vector([input.cylinder_mesh_parameter_list(5) input.cylinder_mesh_parameter_list(6) input.cylinder_mesh_parameter_list(7)]',input.cylinder_mesh_parameter_list(8));



if input.surface_is_cylinder_flag & input.make_cylndrical_pcb

if ~strcmp(input.pcb_interconnection_method,'spiral_in_out') % first layer spiral in; second layer spiral out

for part_ind=1:numel(coil_parts)

%Calucate the boundaries of the un-rolled cylinder
% rot_cylinder_vertices=rot_mat*coil_parts(part_ind).coil_mesh.vertices;
% phi_coords_mesh=atan2(rot_cylinder_vertices(2,:),rot_cylinder_vertices(1,:));
%unrolled_cylinder=[phi_coords_mesh.*cylinder_radius; rot_cylinder_vertices(3,:)];


%Rotate the wire on the cylinder 
aligned_wire_path=(rot_mat*coil_parts(part_ind).wire_path.v);

%for each point calculate the phi angle
phi_coord=atan2(aligned_wire_path(2,:),aligned_wire_path(1,:));
layout_2d=[phi_coord.*cylinder_radius; aligned_wire_path(3,:)];
if layout_2d(1,1)~=layout_2d(1,end)&layout_2d(2,1)~=layout_2d(2,end)
layout_2d=[layout_2d layout_2d(:,1)];
end


segment_starts=sort([1 find(diff(coil_parts(part_ind).points_to_shift)==1) find(diff(coil_parts(part_ind).points_to_shift)==-1)+1]);


%Generate the track shapes for the indivial wire parts
warning('off','all');
pcb_parts(numel(segment_starts)).track_shape=[];
pcb_parts(numel(segment_starts)).polygon_track=[];
for wire_part_ind=1:numel(pcb_parts)
if wire_part_ind<numel(pcb_parts)
wire_part_inds=segment_starts(wire_part_ind):segment_starts(wire_part_ind+1);
else
wire_part_inds=segment_starts(wire_part_ind):size(layout_2d,2);
end
segment_points=layout_2d(:,wire_part_inds);
long_vecs=segment_points(:,2:end)-segment_points(:,1:end-1);
long_vecs=[long_vecs long_vecs(:,end)];
long_vecs=long_vecs./repmat(vecnorm(long_vecs),[2 1]);
ortho_vecs=[long_vecs(2,:); long_vecs(1,:)*(-1)];
ortho_vecs=ortho_vecs./repmat(vecnorm(ortho_vecs),[2 1]);
pcb_parts(wire_part_ind).track_shape=[segment_points+ortho_vecs.*(pcb_track_width/2) fliplr(segment_points)-fliplr(ortho_vecs).*(pcb_track_width/2)];
pcb_parts(wire_part_ind).track_shape=[pcb_parts(wire_part_ind).track_shape pcb_parts(wire_part_ind).track_shape(:,1)];
pcb_parts(wire_part_ind).polygon_track=polyshape(pcb_parts(wire_part_ind).track_shape(1,:),pcb_parts(wire_part_ind).track_shape(2,:));
end
warning('on','all');
coil_parts(part_ind).pcb_tracks.upper_layer.group_layouts.wire_parts=pcb_parts;


end

else       %Generate the pcb form the spiral in/out tracks


upper_layer(numel(coil_parts)).group_layouts=[];
lower_layer(numel(coil_parts)).group_layouts=[];

for part_ind=1:numel(coil_parts)

%Calucate the boundaries of the un-rolled cylinder
rot_cylinder_vertices=rot_mat*coil_parts(part_ind).coil_mesh.vertices;
phi_coords_mesh=atan2(rot_cylinder_vertices(2,:),rot_cylinder_vertices(1,:));
unrolled_cylinder=[phi_coords_mesh.*cylinder_radius; rot_cylinder_vertices(3,:)];

upper_layer(part_ind).group_layouts(numel(coil_parts(part_ind).connected_group)).wire_parts=[];
lower_layer(part_ind).group_layouts(numel(coil_parts(part_ind).connected_group)).wire_parts=[];

for group_ind=1:numel(coil_parts(part_ind).connected_group)

track_spiral_in=coil_parts(part_ind).connected_group(group_ind).spiral_in.v;
track_spiral_out=coil_parts(part_ind).connected_group(group_ind).spiral_out.v;
aligned_wire_path_spiral_in=rot_mat*track_spiral_in;
aligned_wire_path_spiral_out=rot_mat*track_spiral_out;

%for each point calculate the phi angle
phi_coord_spiral_in=atan2(aligned_wire_path_spiral_in(2,:),aligned_wire_path_spiral_in(1,:));
phi_coord_spiral_out=atan2(aligned_wire_path_spiral_out(2,:),aligned_wire_path_spiral_out(1,:));
layout_2d_spiral_in=[phi_coord_spiral_in; aligned_wire_path_spiral_in(3,:)];
layout_2d_spiral_out=[phi_coord_spiral_out; aligned_wire_path_spiral_out(3,:)];

%Add common points at the start and end of both in/out layers
point_1=(layout_2d_spiral_in(:,1)+layout_2d_spiral_out(:,end))./2;
point_2=(layout_2d_spiral_in(:,end)+layout_2d_spiral_out(:,1))./2;
center_position=(layout_2d_spiral_in(:,1)+layout_2d_spiral_in(:,end)+layout_2d_spiral_out(:,1)+layout_2d_spiral_out(:,end))./4;
%Shit those points a little to avoid overlapps with the interconnections
point_1=point_1+(point_1-center_position).*(input.pcb_spiral_end_shift_factor/100);
point_2=point_2+(point_2-center_position).*(input.pcb_spiral_end_shift_factor/100);
layout_2d_spiral_in=[point_1 point_1 layout_2d_spiral_in(:,2:end-1) point_2 point_2];
layout_2d_spiral_out=[point_2 point_2 layout_2d_spiral_out(:,2:end-1) point_1 point_1];

for group_layer={'upper' 'lower'}


if strcmp(group_layer,'upper')
%Rotate the wire on the cylinder
layout_2d=layout_2d_spiral_in;
else
layout_2d=layout_2d_spiral_out;
end

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
layout_2d(:,negative_wrap+1)=layout_2d(:,negative_wrap+1)+[ones(1,numel(negative_wrap)).*2.*pi;  zeros(1,numel(negative_wrap))];

%Generate a cut border on the pi,-pi phase wrap
cut_rectangle=[-pi pi pi -pi; ...
                                        max(aligned_wire_path_spiral_in(3,:))+abs(max(aligned_wire_path_spiral_in(3,:))).*0.1 ...
                                        max(aligned_wire_path_spiral_in(3,:))+abs(max(aligned_wire_path_spiral_in(3,:))).*0.1 ...
                                        min(aligned_wire_path_spiral_in(3,:))-abs(max(aligned_wire_path_spiral_in(3,:))).*0.1 ...
                                        min(aligned_wire_path_spiral_in(3,:))-abs(max(aligned_wire_path_spiral_in(3,:))).*0.1];
cut_rectangle=[cut_rectangle cut_rectangle(:,1)];

%Scale the geometrie to the cylinder radius
layout_2d(1,:)=layout_2d(1,:).*cylinder_radius;
cut_rectangle(1,:)=cut_rectangle(1,:).*cylinder_radius;


%Generate the wire segmential parts between the wraps
clear pcb_parts
pcb_parts(numel(full_wrap_spart_inds)-1).uv=[];
pcb_parts(numel(full_wrap_spart_inds)-1).ind1=[];
pcb_parts(numel(full_wrap_spart_inds)-1).ind2=[];
pcb_parts(numel(full_wrap_spart_inds)-1).track_shape=[];
pcb_parts(numel(full_wrap_spart_inds)-1).polygon_track=[];

for point_ind=1:numel(full_wrap_spart_inds)-1
pcb_parts(point_ind).uv=layout_2d(:,full_wrap_spart_inds(point_ind)+1:full_wrap_spart_inds(point_ind+1));
pcb_parts(point_ind).ind1=full_wrap_spart_inds(point_ind)+1;
pcb_parts(point_ind).ind2=full_wrap_spart_inds(point_ind+1);
end

%Add the cut points for a clean cut
for wrap_ind=1:numel(pcb_parts)
intersection_cut=find_segment_intersections(pcb_parts(wrap_ind).uv,cut_rectangle);
is_real_cut_ind=find(~isnan([intersection_cut(:).segment_inds]));
if ~isempty(is_real_cut_ind)
wire_part_points=pcb_parts(wrap_ind).uv;
uv_point=intersection_cut(is_real_cut_ind(1)).uv;
cut_segment_ind=intersection_cut(is_real_cut_ind(1)).segment_inds;
if cut_segment_ind~=1
pcb_parts(wrap_ind).uv=[wire_part_points(:,1:cut_segment_ind) uv_point wire_part_points(:,cut_segment_ind+1:end-1)];
end
end
end
%add a clean cut also for the second open end of the wire_parts
for wrap_ind=2:numel(pcb_parts)
if pcb_parts(wrap_ind-1).uv(1,end)>0
pcb_parts(wrap_ind).uv=[pcb_parts(wrap_ind-1).uv(:,end)-[2*pi*cylinder_radius 0]' pcb_parts(wrap_ind).uv];
else
pcb_parts(wrap_ind).uv=[pcb_parts(wrap_ind-1).uv(:,end)+[2*pi*cylinder_radius 0]' pcb_parts(wrap_ind).uv];
end
end

pcb_parts(arrayfun(@(x) size(pcb_parts(x).uv,2),1:numel(pcb_parts))<2)=[]; %delete fragments


%Generate the track shapes for the indivial wire parts
warning('off','all');
for wire_part_ind=1:numel(pcb_parts)
if size(pcb_parts(wire_part_ind).uv,2)>5
smoothed_track=[pcb_parts(wire_part_ind).uv(:,1) smooth_track_by_folding(pcb_parts(wire_part_ind).uv(:,2:end-1),3) pcb_parts(wire_part_ind).uv(:,end)];
else
smoothed_track=pcb_parts(wire_part_ind).uv;
end
long_vecs=smoothed_track(:,2:end)-smoothed_track(:,1:end-1);
long_vecs=[long_vecs long_vecs(:,end)];
long_vecs=long_vecs./repmat(vecnorm(long_vecs),[2 1]);
ortho_vecs=[long_vecs(2,:); long_vecs(1,:)*(-1)];
ortho_vecs=ortho_vecs./repmat(vecnorm(ortho_vecs),[2 1]);
pcb_parts(wire_part_ind).track_shape=[smoothed_track+ortho_vecs.*(pcb_track_width/2) fliplr(smoothed_track)-fliplr(ortho_vecs).*(pcb_track_width/2)];
pcb_parts(wire_part_ind).track_shape=[pcb_parts(wire_part_ind).track_shape pcb_parts(wire_part_ind).track_shape(:,1)];
pcb_parts(wire_part_ind).polygon_track=polyshape(pcb_parts(wire_part_ind).track_shape(1,:),pcb_parts(wire_part_ind).track_shape(2,:));
end
warning('on','all');

%write the outputs
if strcmp(group_layer,'upper')
upper_layer(part_ind).group_layouts(group_ind).wire_parts=pcb_parts;
%writetable(array2table(wire_part.track_shape','VariableNames',{'x','y'}),"upper_layer_part"+num2str(part_ind)+"_"+"group"+num2str(group_ind)+"_"+"wire_part"+num2str(wire_part_ind));
else
lower_layer(part_ind).group_layouts(group_ind).wire_parts=pcb_parts;
%writetable(array2table(wire_part.track_shape','VariableNames',{'x','y'}),"lower_layer_part"+num2str(part_ind)+"_"+"group"+num2str(group_ind)+"_"+"wire_part"+num2str(wire_part_ind));
end

end

end

coil_parts(part_ind).pcb_tracks.upper_layer=upper_layer;
coil_parts(part_ind).pcb_tracks.lower_layer=lower_layer;

%Save the tracks as a vector file
%save_pcb_tracks_as_svg(coil_parts(part_ind).pcb_tracks,input.field_shape_function,'pcb_layout',part_ind,unrolled_cylinder,input.output_directory);

end

end





end

end

function save_pcb_tracks_as_svg(input_poly,coil_name,filename,part_num,unrolled_mesh_points,output_directory)
%Save the pcb tracked as a svg file

mesh_boundary_points=[max(unrolled_mesh_points(1,:)) max(unrolled_mesh_points(1,:)) min(unrolled_mesh_points(1,:)) min(unrolled_mesh_points(1,:)); max(unrolled_mesh_points(2,:)) min(unrolled_mesh_points(2,:)) min(unrolled_mesh_points(2,:)) max(unrolled_mesh_points(2,:))];
mesh_boundary_points=[mesh_boundary_points mesh_boundary_points(:,1)].*1000;
mesh_size=[max(unrolled_mesh_points(1,:)) min(unrolled_mesh_points(1,:)) max(unrolled_mesh_points(2,:)) min(unrolled_mesh_points(2,:))].*1000;

if ispc
filename_full=strcat(output_directory,'\',coil_name,'_',filename,'_part',erase(num2str(part_num),{'.' '^' ','}),'.svg');
else
filename_full=strcat(output_directory,'/',coil_name,'_',filename,'_part',erase(num2str(part_num),{'.' '^' ','}),'.svg');
end

% Create output SVG header
fid = fopen(filename_full,'w');
fprintf(fid,'<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n');
fprintf(fid,['<svg xmlns="http://www.w3.org/2000/svg" version="1.200000" width="' num2str(mesh_size(2)-mesh_size(1)) 'mm" height="' num2str(mesh_size(4)-mesh_size(3)) 'mm" xmlns:xlink="http://www.w3.org/1999/xlink">\n']);
%fprintf(fid,['<svg xmlns="http://www.w3.org/2000/svg" version="1.200000" width="100%%" height="100%%" viewBox="0 0 ' int2str(im_size(2) + 1) ' ' int2str(im_size(1) + 1) '" xmlns:xlink="http://www.w3.org/1999/xlink">']);

%for the upper layer (spiral in)
for group_ind=1:numel(input_poly.upper_layer.group_layouts)
for wire_part_ind=1:numel(input_poly.upper_layer.group_layouts(group_ind).wire_parts)
x1 = input_poly.upper_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track.Vertices(:,1).*1000;
y1 = input_poly.upper_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track.Vertices(:,2).*1000;
fprintf(fid,'<polygon points="');
for i = 1:length( input_poly.upper_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track.Vertices)
% Print to SVG
fprintf(fid,strcat(num2str(x1(i)),',',num2str(y1(i))," "));
end
fprintf(fid,'" stroke="none" fill="blue" stroke-width="0" shape-rendering="crispEdges" />\n');
end
end

%for the lower layer (spiral out)
for group_ind=1:numel(input_poly.lower_layer.group_layouts)
for wire_part_ind=1:numel(input_poly.lower_layer.group_layouts(group_ind).wire_parts)
x1 = input_poly.lower_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track.Vertices(:,1).*1000;
y1 = input_poly.lower_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track.Vertices(:,2).*1000;
fprintf(fid,'<polygon points="');
for i = 1:length( input_poly.lower_layer.group_layouts(group_ind).wire_parts(wire_part_ind).polygon_track.Vertices)
% Print to SVG
fprintf(fid,strcat(num2str(x1(i)),',',num2str(y1(i))," "));
end
fprintf(fid,'" stroke="none" fill="red" stroke-width="0" shape-rendering="crispEdges" />\n');
end
end

%Print the outer bounary of the mesh
for i = 1:size(mesh_boundary_points,2)-1
fprintf(fid,strcat('<line x1="',num2str(mesh_boundary_points(1,i)),'" y1="',num2str(mesh_boundary_points(2,i)),'" x2="',num2str(mesh_boundary_points(1,i+1)),'" y2="',num2str(mesh_boundary_points(2,i+1)),'" stroke="red" stroke-width="0.1"/>\n'));
end

fprintf(fid,'"/>\n');


% End SVG
%fprintf(fid,'</g>\n');
fprintf(fid,'</svg>');


end



function [P,intersect_edge_inds] = InterX(L1,varargin)
    if nargin == 1
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i)
        P = zeros(2,0);
        intersect_edge_inds=[];
        return; 
    end
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    intersect_edge_inds=unique(sort([i' j],2),'rows');
    %intersect_edge_inds=[i' j];
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
            
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
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