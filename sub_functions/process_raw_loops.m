function coil_parts= process_raw_loops(coil_parts,input)
%take care of loops crossing boundaries or off boundary loops


min_point_loop_number=input.min_point_loop_number;
area_perimeter_deletion_ratio=input.area_perimeter_deletion_ratio;
max_allowed_angle_within_coil_track=input.max_allowed_angle_within_coil_track;
min_allowed_angle_within_coil_track=input.min_allowed_angle_within_coil_track;
tiny_segment_length_percentage=input.tiny_segment_length_percentage;

for part_ind=1:numel(coil_parts)

planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.v);


%primitive deleting small loops
loops_to_delete=arrayfun(@(x) size(coil_parts(part_ind).contour_lines(x).uv,2)<min_point_loop_number,1:numel(coil_parts(part_ind).contour_lines));
coil_parts(part_ind).contour_lines(loops_to_delete)=[];


%Smooth the loops
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
%coil_parts(part_ind).contour_lines(loop_ind).uv = equilize_point_distances_spline(coil_parts(part_ind).contour_lines(loop_ind).uv,1);
if coil_parts(part_ind).contour_lines(loop_ind).uv(1,end)~=coil_parts(part_ind).contour_lines(loop_ind).uv(1,1) & coil_parts(part_ind).contour_lines(loop_ind).uv(2,end)~=coil_parts(part_ind).contour_lines(loop_ind).uv(2,1) 
coil_parts(part_ind).contour_lines(loop_ind).uv=[coil_parts(part_ind).contour_lines(loop_ind).uv coil_parts(part_ind).contour_lines(loop_ind).uv(:,1)]; %close the loops
end
end

%Generate the curved coordinates
for loop_num=1:numel(coil_parts(part_ind).contour_lines)
[coil_parts(part_ind).contour_lines(loop_num).v,coil_parts(part_ind).contour_lines(loop_num).uv]=uv_to_xyz(coil_parts(part_ind).contour_lines(loop_num).uv,planary_mesh,curved_mesh);
end


end




end


% % % %first delete loops which are either to small or 
% % % %have a very high boundary-length-to-surface ratio
% % % %loops which such ratios will not contribute significantly to the target
% % % %field since of narrow opposing currents on the boundary
% % % loop_surface_areas=zeros(1,numel(coil_parts(part_ind).contour_lines));
% % % loop_perimeter=zeros(1,numel(coil_parts(part_ind).contour_lines));
% % % 
% % % warning off all;
% % % for loop_ind=1:numel(coil_parts(part_ind).contour_lines) 
% % %     loop_points=coil_parts(part_ind).contour_lines(loop_ind).uv;
% % %     poly_object=polyshape(loop_points');
% % %     loop_surface_areas(loop_ind)=area(poly_object);
% % %     loop_perimeter(loop_ind)=perimeter(poly_object);
% % % end
% % % warning on all;
% % % perimeter_to_area_ratio=loop_perimeter.*loop_perimeter./loop_surface_areas;
% % % loops_to_delete=find(perimeter_to_area_ratio>area_perimeter_deletion_ratio*mean(perimeter_to_area_ratio));
% % % coil_parts(part_ind).contour_lines(loops_to_delete)=[];

% % % % %check for sharp angles which might be removed
% % % % %delete points which build a degenerate angle of almost zero i.e. the point
% % % % %track does not change when these points will be removed
% % % % 
% % % % track.point_inds=[];
% % % % track.vecs=[];
% % % % track.angles=[];
% % % %  
% % % % 
% % % % deleted_loops=0;
% % % % 
% % % % for loop_ind=1:numel(coil_parts(part_ind).contour_lines) 
% % % %     point_to_delete=ones(1,size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2));
% % % %     while any(point_to_delete)  
% % % % 	track.angle_point_inds=[[circshift(1:size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2),1)]' [1:size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2)]' [circshift(1:size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2),-1)]'];
% % % %     first_vecs=[coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,2))-coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,1))];
% % % %     second_vecs=[coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,3))-coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,2))];
% % % %     first_vecs=first_vecs./vecnorm(first_vecs,2,1);
% % % %     second_vecs=second_vecs./vecnorm(second_vecs,2,1);
% % % %     dot_p=dot(first_vecs, second_vecs);
% % % %     dot_p(dot_p>1)=1;
% % % %     dot_p(dot_p<-1)=-1;
% % % %     track.angles= acos(dot_p)./pi.*180;
% % % %     %Remove the points of the track that build to sharp angles, dull angles
% % % %     %or tiny legs
% % % %     point_to_delete=abs(track.angles)>max_allowed_angle_within_coil_track | abs(track.angles)<min_allowed_angle_within_coil_track ;
% % % %     if ~any(point_to_delete)
% % % %         break;
% % % %     end
% % % %     coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,point_to_delete)=[];
% % % %     if isempty(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv) || size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2)<2 % delete the loop entry if there are no points left
% % % %     coil_parts(part_ind).contour_lines(loop_ind-deleted_loops)=[];
% % % %     deleted_loops=deleted_loops+1;
% % % %     break;
% % % %     end
% % % %     %Recalculate the angles
% % % % 	track.angle_point_inds=[[circshift(1:size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2),1)]' [1:size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2)]' [circshift(1:size(coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv,2),-1)]'];
% % % %     first_vecs=[coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,2))-coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,1))];
% % % %     second_vecs=[coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,3))-coil_parts(part_ind).contour_lines(loop_ind-deleted_loops).uv(:,track.angle_point_inds(:,2))];
% % % %     first_vecs=first_vecs./vecnorm(first_vecs,2,1);
% % % %     second_vecs=second_vecs./vecnorm(second_vecs,2,1);
% % % %     dot_p=dot(first_vecs, second_vecs);
% % % %     dot_p(dot_p>1)=1;
% % % %     dot_p(dot_p<-1)=-1;
% % % %     track.angles= acos(dot_p)./pi.*180;
% % % %     point_to_delete=abs(track.angles)>max_allowed_angle_within_coil_track | abs(track.angles)<min_allowed_angle_within_coil_track;
% % % %     end
% % % % end
% % % % 
% % % % 
% % % % %delete points which build a segment which is extremly short
% % % % %limit the point to point distane in the 2d representation
% % % % for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
% % % % coil_parts(part_ind).contour_lines(loop_ind).uv=[coil_parts(part_ind).contour_lines(loop_ind).uv coil_parts(part_ind).contour_lines(loop_ind).uv(:,1)];
% % % % seq_dist=vecnorm(coil_parts(part_ind).contour_lines(loop_ind).uv(:,2:end)-coil_parts(part_ind).contour_lines(loop_ind).uv(:,1:end-1));
% % % % seq_is_to_short=find([false seq_dist<mean(seq_dist)*tiny_segment_length_percentage/100]);
% % % % %Replace the two nodes of the short segments with one node in the middle
% % % % if ~isempty(seq_is_to_short)
% % % % num_replaced_segments=0;
% % % % for repl_ind=seq_is_to_short
% % % % replacement_point=[coil_parts(part_ind).contour_lines(loop_ind).uv(:,repl_ind-num_replaced_segments)+coil_parts(part_ind).contour_lines(loop_ind).uv(:,repl_ind-1-num_replaced_segments)]./2;
% % % % coil_parts(part_ind).contour_lines(loop_ind).uv(:,repl_ind-num_replaced_segments-1)=replacement_point;
% % % % coil_parts(part_ind).contour_lines(loop_ind).uv(:,repl_ind-num_replaced_segments)=[];
% % % % num_replaced_segments=num_replaced_segments+1;
% % % % end
% % % % end
% % % % end
    












