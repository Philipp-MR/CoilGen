function coil_parts= process_raw_loops(coil_parts,min_point_loop_number,area_perimeter_deletion_ratio,max_allowed_angle_within_coil_track,min_allowed_angle_within_coil_track,tiny_segment_length_percentage)
%take care of loops crossing boundaries or off boundary loops




for part_ind=1:numel(coil_parts)

%first delete loops which are either to small or 
%have a very high boundary-length-to-surface ratio
%loops which such ratios will not contribute significantly to the target
%field since of narrow opposing currents on the boundary
loop_surface_areas=zeros(1,numel(coil_parts(part_ind).contour_lines));
loop_perimeter=zeros(1,numel(coil_parts(part_ind).contour_lines));

warning off all;
for loop_ind=1:numel(coil_parts(part_ind).contour_lines) 
    loop_points=coil_parts(part_ind).contour_lines(loop_ind).uv;
    poly_object=polyshape(loop_points');
    loop_surface_areas(loop_ind)=area(poly_object);
    loop_perimeter(loop_ind)=perimeter(poly_object);
end
warning on all;
perimeter_to_area_ratio=loop_perimeter.*loop_perimeter./loop_surface_areas;
loops_to_delete=find(perimeter_to_area_ratio>area_perimeter_deletion_ratio*mean(perimeter_to_area_ratio));
coil_parts(part_ind).contour_lines(loops_to_delete)=[];


%primitive deleting small loops
loops_to_delete=arrayfun(@(x) size(coil_parts(part_ind).contour_lines(x).uv,2)<min_point_loop_number,1:numel(coil_parts(part_ind).contour_lines));
coil_parts(part_ind).contour_lines(loops_to_delete)=[];





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
    

%Smooth the loops
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
%coil_parts(part_ind).contour_lines(loop_ind).uv = equilize_point_distances_spline(coil_parts(part_ind).contour_lines(loop_ind).uv,1);
coil_parts(part_ind).contour_lines(loop_ind).uv=[coil_parts(part_ind).contour_lines(loop_ind).uv coil_parts(part_ind).contour_lines(loop_ind).uv(:,1)]; %close the loops
end

%Generate the curved coordinates
coil_parts(part_ind).contour_lines=create_curved_contour_lines(coil_parts(part_ind).contour_lines,triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv'),triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices'));

end


function contour_lines_curved=create_curved_contour_lines(contour_lines,planary_mesh,curved_mesh)
%Create the contour lines on the curved surface
contour_lines_curved=contour_lines;
%contour_lines_curved = rmfield(contour_lines_curved,'uv');
% for loop_num=1:numel(contour_lines)
% contour_lines_curved(loop_num).point_coordinates=zeros(3,size(contour_lines(loop_num).uv,2));
% end
for loop_num=1:numel(contour_lines)
[contour_lines_curved(loop_num).point_coordinates,contour_lines(loop_num).uv]=uv_to_xyz(contour_lines(loop_num).uv,planary_mesh,curved_mesh);
%disp(strcat('Loop Number',' ',num2str(loop_num),' Point Number ',num2str(point_num)));      
end

end


function sampled_points = equilize_point_distances_spline(unsampled_points,upsampling_factor)
%Equilize the point the distances, keeping the same number of points
number_points=size(unsampled_points,2);
t=0:(1/(number_points*upsampling_factor)):1;
t = t(:);
nt = numel(t);
ndim = 2;
%delete non-unique points
ind_point_x=diff(unsampled_points(1,:))==0;
ind_point_y=diff(unsampled_points(2,:))==0;
unique_points=~(ind_point_x|ind_point_y);
unsampled_points=unsampled_points(:,unique_points')';
pt = NaN(nt,ndim);
chordlen = sqrt(sum(diff(unsampled_points,[],1).^2,2));
chordlen = chordlen/sum(chordlen);
cumarc = [0;cumsum(chordlen)];
spline_interpol = cell(1,ndim);
spld = spline_interpol;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for spl_dim = 1:ndim
spline_interpol{spl_dim} = spline(cumarc,unsampled_points(:,spl_dim));
nc = numel(spline_interpol{spl_dim}.coefs);
if nc < 4
spline_interpol{spl_dim}.coefs = [zeros(1,4-nc),spline_interpol{spl_dim}.coefs];
spline_interpol{spl_dim}.order = 4;
end
xp = spline_interpol{spl_dim};
xp.coefs = xp.coefs*diffarray;
xp.order = 3;
spld{spl_dim} = xp;
end
cumarc = spline_interpol{1}.breaks;
n = numel(cumarc);
polyarray = zeros(ndim,3);
segment_length = zeros(n-1,1);
opts = odeset('reltol',1.e-9);
for iii = 1:spline_interpol{1}.pieces
for jjj = 1:ndim
polyarray(jjj,:) = spld{jjj}.coefs(iii,:);
end
[tout,yout] = ode45(@(t,y) segkernel(t,y,ndim,polyarray),[0,chordlen(iii)],0,opts); %#ok
segment_length(iii) = yout(end);
end
totalsplinelength = sum(segment_length);
cumulative_segment_length = [0;cumsum(segment_length)];
[junk,tbins] = histc(t*totalsplinelength,cumulative_segment_length); %#ok
tbins((tbins <= 0) | (t <= 0)) = 1;
tbins((tbins >= n) | (t >= 1)) = n - 1;
s = totalsplinelength*t;
opts = odeset('reltol',1.e-9,'events',@ode_events);
ti = t;
for ii = 1:nt
si = s(ii) - cumulative_segment_length(tbins(ii));
for jj = 1:ndim
polyarray(jj,:) = spld{jj}.coefs(tbins(ii),:);
end
[tout,yout,te,ye] = ode45(@(t,y) segkernel(t,y,ndim,polyarray),[0,chordlen(tbins(ii))],-si,opts); %#ok
if ~isempty(te)
ti(ii) = te(1) + cumarc(tbins(ii));
else
if abs(yout(1)) < abs(yout(end))
ti(ii) = tout(1) + cumarc(tbins(ii));
else
ti(ii) = tout(end) + cumarc(tbins(ii));
end
end
end
for L = 1:ndim
pt(:,L) = ppval(spline_interpol{L},ti);
end
dudt = zeros(nt,ndim);
for L = 1:ndim
dudt(:,L) = ppval(spld{L},ti);
end
function val = segkernel(t,y,dim_num,poly_array) %#ok
val = zeros(size(t));
for k = 1:dim_num
val = val + polyval(poly_array(k,:),t).^2;
end
val = sqrt(val);
end
function [value,isterminal,direction] = ode_events(t,y) %#ok
value = y;
isterminal = ones(size(y));
direction = ones(size(y));
end
sampled_points=transpose(pt);
end


end







