function coil_parts = calculate_opening_directions(coil_parts,b_0_direction,opening_gap)
%define the realtive alignment of the opening cuts


expand_factor=2;  % factor for which the cut areas will be scaled in length


coil_parts(numel(coil_parts)).rectangle_cuts=[];

for part_ind=1:numel(coil_parts)

% the cuts for the differnt groups shall be on a consistent line (if possible)
clear cut_points high_cuts low_cuts
cut_points(numel(coil_parts(part_ind).groups))=struct();
high_cuts(numel(coil_parts(part_ind).groups))=struct();
low_cuts(numel(coil_parts(part_ind).groups))=struct();

%define the realtive alignment of the opening cuts
% the cuts for the differnt groups shall be on a consistent line (if possible)

%calc the boundaries of the volume which contains all points
all_points_v=[];
all_points_uv=[];
for group_ind=1:numel(coil_parts(part_ind).groups)
for loop_ind=1:numel(coil_parts(part_ind).groups(group_ind).loops  )
all_points_v=[all_points_v coil_parts(part_ind).groups(group_ind).loops(loop_ind).point_coordinates];
all_points_uv=[all_points_uv coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv];
end
end
total_center_v=mean(all_points_v,2);
total_center_uv=mean(all_points_uv,2);
diag_points_v=[max(all_points_v(1,:)) max(all_points_v(2,:)) max(all_points_v(3,:)); min(all_points_v(1,:)) min(all_points_v(2,:)) min(all_points_v(3,:))];
vol_diagonal_v=norm(diag_points_v(2,:)-diag_points_v(1,:));
diag_points_uv=[max(all_points_uv(1,:)) max(all_points_uv(2,:)); min(all_points_uv(1,:)) min(all_points_uv(2,:))];
vol_diagonal_uv=norm(diag_points_uv(2,:)-diag_points_uv(1,:));


for group_ind=1:numel(coil_parts(part_ind).groups)

loop_num=numel(coil_parts(part_ind).groups(group_ind).loops);
    

cut_points(group_ind).v(loop_num).coords=[];
cut_points(group_ind).uv(loop_num).coords=[];
high_cuts(group_ind).v(loop_num).coords=[];
high_cuts(group_ind).uv(loop_num).coords=[];
low_cuts(group_ind).v(loop_num).coords=[];
low_cuts(group_ind).uv(loop_num).coords=[];


%loop_normal=calc_mean_loop_normal(groups,group_ind,group_center,total_center_v);
loop_normal=calc_mean_loop_normal2(coil_parts(part_ind).groups,group_ind,coil_parts(part_ind).coil_mesh);

%test if loop normal and b0_direction are not independent enough
if norm(cross(b_0_direction./norm(b_0_direction),loop_normal))<0.05
%select a vector from the first point to the group center as the second
%vector for the cut plane


alternative_center_normal=b_0_direction;
[~,coord_vec_min_ind]=min([dot(b_0_direction,[1;0;0]) dot(b_0_direction,[0;1;0]) dot(b_0_direction,[0;0;1])]);

switch coord_vec_min_ind
case 1
alternative_center_normal=cross(b_0_direction,[1, 0, 0])';
alternative_center_normal=alternative_center_normal./norm(alternative_center_normal);
case 2 
alternative_center_normal=cross(b_0_direction,[0, 1, 0])';
alternative_center_normal=alternative_center_normal./norm(alternative_center_normal);
case 3
alternative_center_normal=cross(b_0_direction,[0, 0, 1])';
alternative_center_normal=alternative_center_normal./norm(alternative_center_normal);
end
%rotate the second vector around the b0 direction
%v1=rotationVectorToMatrix(b_0_direction.*opening_cut_rotation_angle)*v1;
cut_plane_direction=cross(b_0_direction,alternative_center_normal);
cut_plane_direction=cut_plane_direction./norm(cut_plane_direction);
else
cut_plane_direction=cross(b_0_direction,loop_normal);
cut_plane_direction=cut_plane_direction./norm(cut_plane_direction);
end

%find the cut points of the cut plane with the loops
for loop_ind=1:loop_num
for point_ind=2:size(coil_parts(part_ind).groups(group_ind).loops(loop_ind).point_coordinates,2)
point_a=coil_parts(part_ind).groups(group_ind).loops(loop_ind).point_coordinates(:,point_ind-1)';
point_b=coil_parts(part_ind).groups(group_ind).loops(loop_ind).point_coordinates(:,point_ind)';
[cut_p,cut_flag]=plane_line_intersect(cut_plane_direction,coil_parts(part_ind).group_centers.v(:,group_ind)',point_a,point_b);
if cut_flag==1 %test if there is a cut between the line points
cut_points(group_ind).v(loop_ind).coords=[cut_points(group_ind).v(loop_ind).coords cut_p'];
%build the corresponding uv point
cut_point_ratio=norm(point_a-cut_p)/norm(point_a-point_b);
point_a_uv=coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv(:,point_ind-1);
point_b_uv=coil_parts(part_ind).groups(group_ind).loops(loop_ind).uv(:,point_ind);
cut_poin_uv=point_a_uv+(point_b_uv-point_a_uv)*cut_point_ratio;
cut_points(group_ind).uv(loop_ind).coords=[cut_points(group_ind).uv(loop_ind).coords cut_poin_uv];
end
end
end

%sort the cutpoints into high and low among the cut plane
high_cuts(group_ind).v(loop_num).coords=[];
low_cuts(group_ind).uv(loop_num).coords=[];

%define a "seperation" plane which is the cut plane rotated for 90°
seperaton_plane_normal=(cut_points(group_ind).v(1).coords(:,1)-cut_points(group_ind).v(1).coords(:,2))./norm(cut_points(group_ind).v(1).coords(:,1)-cut_points(group_ind).v(1).coords(:,2));
%Make sure that the orientation is well defined
if dot(seperaton_plane_normal,b_0_direction) <0
seperaton_plane_normal=seperaton_plane_normal.*(-1);
end
%Seperate the cut points by checking wether they are above or below the
%separtion plane
for loop_ind=1:loop_num
    
cut_orientations=sign(dot(cut_points(group_ind).v(loop_ind).coords-coil_parts(part_ind).group_centers.v(:,group_ind),repmat(seperaton_plane_normal,[1 size(cut_points(group_ind).v(loop_ind).coords,2)])));
[~,cut_dist_order]=sort(vecnorm(cut_points(group_ind).v(loop_ind).coords-coil_parts(part_ind).group_centers.v(:,group_ind)));
high_cut_ind=intersect(cut_dist_order,find(cut_orientations==1),'stable');
low_cut_ind=intersect(cut_dist_order,find(cut_orientations==-1),'stable');

high_cuts(group_ind).v(loop_ind).coords=cut_points(group_ind).v(loop_ind).coords(:,high_cut_ind(1));
low_cuts(group_ind).v(loop_ind).coords=cut_points(group_ind).v(loop_ind).coords(:,low_cut_ind(1));
high_cuts(group_ind).uv(loop_ind).coords=cut_points(group_ind).uv(loop_ind).coords(:,high_cut_ind(1));
low_cuts(group_ind).uv(loop_ind).coords=cut_points(group_ind).uv(loop_ind).coords(:,low_cut_ind(1));


% if dot(cut_points(group_ind).v(loop_ind).coords(:,1)-coil_parts(part_ind).group_centers.v(:,group_ind),seperaton_plane_normal)>0
% high_cuts(group_ind).v(loop_ind).coords=cut_points(group_ind).v(loop_ind).coords(:,1);
% low_cuts(group_ind).v(loop_ind).coords=cut_points(group_ind).v(loop_ind).coords(:,2);
% high_cuts(group_ind).uv(loop_ind).coords=cut_points(group_ind).uv(loop_ind).coords(:,1);
% low_cuts(group_ind).uv(loop_ind).coords=cut_points(group_ind).uv(loop_ind).coords(:,2);
% else
% high_cuts(group_ind).v(loop_ind).coords=cut_points(group_ind).v(loop_ind).coords(:,2);
% low_cuts(group_ind).v(loop_ind).coords=cut_points(group_ind).v(loop_ind).coords(:,1);
% high_cuts(group_ind).uv(loop_ind).coords=cut_points(group_ind).uv(loop_ind).coords(:,2);
% low_cuts(group_ind).uv(loop_ind).coords=cut_points(group_ind).uv(loop_ind).coords(:,1);
% end
 end

%sort again for distance to the the group 3d center
[~,sort_inds]=sort(vecnorm([high_cuts(group_ind).uv(:).coords]-coil_parts(part_ind).group_centers.uv(:,group_ind)));
high_cuts(group_ind).v=high_cuts(group_ind).v(sort_inds);
high_cuts(group_ind).uv=high_cuts(group_ind).uv(sort_inds);
[~,sort_inds]=sort(vecnorm([low_cuts(group_ind).uv(:).coords]-coil_parts(part_ind).group_centers.uv(:,group_ind)));
low_cuts(group_ind).v=low_cuts(group_ind).v(sort_inds);
low_cuts(group_ind).uv=low_cuts(group_ind).uv(sort_inds);
end

%assign the final output
coil_parts(part_ind).rectangle_cuts.high=create_cut_rectangle(coil_parts(part_ind).groups,high_cuts,coil_parts(part_ind).coil_mesh,opening_gap,coil_parts(part_ind).group_centers,vol_diagonal_uv,expand_factor);
coil_parts(part_ind).rectangle_cuts.low=create_cut_rectangle(coil_parts(part_ind).groups,low_cuts,coil_parts(part_ind).coil_mesh,opening_gap,coil_parts(part_ind).group_centers,vol_diagonal_uv,expand_factor);


end

function [I,check]=plane_line_intersect(n,V0,P0,P1)
I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
        if N == 0           % The segment lies in plane
            check=2;
            return
        else
            check=0;       %no intersection
            return
        end
end
%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;
if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end
end

function rectangle_points=create_cut_rectangle(group_container,cut_points,coil_mesh,opening_gap,group_center,vol_diagonal,expand_factor)
%createting the "cut rectangles" for the groups of loops, which whom the
%loops will be cut open later
    
rectangle_points(numel(cut_points))=struct();
for group_inds=1:numel(cut_points)
middle_path_points=[cut_points(group_inds).uv(:).coords];
cut_points_left=zeros(2,size(middle_path_points,2));
cut_points_right=zeros(2,size(middle_path_points,2));
%calculate the side points for the cut rectangle
if size(middle_path_points,2)>1
longitudinal_vector=middle_path_points(:,2:end)-middle_path_points(:,1:end-1);
%longitudinal_vector=[longitudinal_vector(:,1) longitudinal_vector];
longitudinal_vector=[ longitudinal_vector longitudinal_vector(:,end)];
othorgonal_vector=[longitudinal_vector(2,:); longitudinal_vector(1,:)*(-1)];
othorgonal_vector=othorgonal_vector./(vecnorm(othorgonal_vector));
for hh=1:size(middle_path_points,2)
local_opening_width=calc_local_opening_gab(coil_mesh,middle_path_points(:,hh),othorgonal_vector(:,hh),opening_gap,vol_diagonal);
cut_points_right(:,hh)=middle_path_points(:,hh)-othorgonal_vector(:,hh)*local_opening_width/2;
cut_points_left(:,hh)=middle_path_points(:,hh)+othorgonal_vector(:,hh)*local_opening_width/2;
end

%exclude point s from left and right cut_points that go in the opposit
%direction than the previous ones:
%for the right points
while true
first_vecs=[cut_points_right(:,2:end-1)-cut_points_right(:,1:end-2)];
second_vecs=[cut_points_right(:,3:end)-cut_points_right(:,2:end-1)];
first_vecs=first_vecs./vecnorm(first_vecs,2,1);
second_vecs=second_vecs./vecnorm(second_vecs,2,1);
dot_p=dot(first_vecs, second_vecs);
dot_p(dot_p>1)=1;
dot_p(dot_p<-1)=-1;
track.angles= acos(dot_p)./pi.*180;
point_to_delete=abs(track.angles)>70;
point_to_delete=[false point_to_delete false];
if ~any(point_to_delete)
    break;
end
cut_points_right(:,point_to_delete)=[];
end
%for the left points
while true
first_vecs=[cut_points_left(:,2:end-1)-cut_points_left(:,1:end-2)];
second_vecs=[cut_points_left(:,3:end)-cut_points_left(:,2:end-1)];
first_vecs=first_vecs./vecnorm(first_vecs,2,1);
second_vecs=second_vecs./vecnorm(second_vecs,2,1);
dot_p=dot(first_vecs, second_vecs);
dot_p(dot_p>1)=1;
dot_p(dot_p<-1)=-1;
track.angles= acos(dot_p)./pi.*180;
point_to_delete=abs(track.angles)>70;
point_to_delete=[false point_to_delete false];
if ~any(point_to_delete)
    break;
end
cut_points_left(:,point_to_delete)=[];
end

%expand the ends of rect points a little to assure overlapp on the outer loops
mean_loop_dist=mean(vecnorm(middle_path_points(:,2:end)-middle_path_points(:,1:end-1)));
cut_points_right(:,end)=cut_points_right(:,end)+[cut_points_right(:,end)-cut_points_right(:,end-1)]./vecnorm([cut_points_right(:,end)-cut_points_right(:,end-1)]).*mean_loop_dist.*expand_factor;
cut_points_left(:,end)=cut_points_left(:,end)+[cut_points_left(:,end)-cut_points_left(:,end-1)]./vecnorm([cut_points_left(:,end)-cut_points_left(:,end-1)]).*mean_loop_dist.*expand_factor;
%compose the cut rectingle with the group center as starting and end point
rect_points= [group_center.uv(:,group_inds)'; cut_points_right'; fliplr(cut_points_left)'; group_center.uv(:,group_inds)']';
else
%if there is only a single loop - create a quadratic cut rectangle
opend_loop_points=group_container(group_inds).loops(1).uv(:,1:end-1);
[~,min_inds]=sort(vecnorm(opend_loop_points-middle_path_points,2,1));    
longitudinal_vector=opend_loop_points(:,min_inds(2))-opend_loop_points(:,min_inds(1));
othorgonal_vector=[longitudinal_vector(2) longitudinal_vector(1)*(-1)]';
othorgonal_vector=othorgonal_vector./(norm(othorgonal_vector));
local_opening_width=calc_local_opening_gab(coil_mesh,middle_path_points,othorgonal_vector,opening_gap,vol_diagonal);
longitudinal_vector=longitudinal_vector./norm(longitudinal_vector).*local_opening_width/2;
othorgonal_vector=othorgonal_vector./(norm(othorgonal_vector)).*local_opening_width/2;
rect_points=[middle_path_points+longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_path_points+longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_path_points-longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_path_points-longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_path_points+longitudinal_vector/2+othorgonal_vector/2];
end
rectangle_points(group_inds).uv=rect_points;
end



end



function loop_normal=calc_mean_loop_normal(group_container,group_ind,group_center,total_center)
single_group_center=group_center.v(:,group_ind);
cross_p=zeros(3,numel(group_container(group_ind).loops));
for iiii=1:numel(group_container(group_ind).loops)
loop_vecs=group_container(group_ind).loops(iiii).point_coordinates(:,2:end)-group_container(group_ind).loops(iiii).point_coordinates(:,1:end-1);
if group_container(group_ind).loops(iiii).current_orientation==-1
loop_vecs=loop_vecs.*(-1);
end
vecs_to_group_center=group_container(group_ind).loops(iiii).point_coordinates(:,2:end)-single_group_center;
cross_p(:,iiii)=mean(cross(loop_vecs,vecs_to_group_center),2);
end
loop_normal=mean(cross_p,2);
loop_normal=loop_normal./norm(loop_normal);
%mirror the normal if does not point toward the center
loop_orientation=dot(total_center-single_group_center,loop_normal);
if loop_orientation<0
loop_normal=loop_normal.*(-1);
end
end


function loop_normal=calc_mean_loop_normal2(group_container,group_ind,coil_mesh)
all_loops=[];
for contour_ind=1:numel(group_container(group_ind).loops)
all_loops=[all_loops group_container(group_ind).loops(contour_ind).uv];
end
curved_mesh=triangulation(coil_mesh.faces',coil_mesh.vertices');
%Calculate vertex normals for later
face_normal=faceNormal(curved_mesh);
%make sure that they are pointing to the outisde of the surface
if mean(dot([vertexNormal(curved_mesh)]',[curved_mesh.Points-mean(curved_mesh.Points)]'))<0
face_normal=face_normal.*(-1);
end
[target_triangle_normal,~] = pointLocation(triangulation(coil_mesh.faces',coil_mesh.uv'),all_loops(1,:)',all_loops(2,:)');
loop_normal=mean(face_normal(target_triangle_normal(~isnan(target_triangle_normal)),:),1)';
loop_normal=loop_normal./norm(loop_normal);
end


% function contour_lines= turn_loops(contour_lines,planary_mesh_matlab_format)
% %Turn the loops in a unifying way by point shifting
% mesh_boundaries=[max(planary_mesh_matlab_format.Points(:,1)); max(planary_mesh_matlab_format.Points(:,2))];
% for loop_inds=1:numel(contour_lines)
% point_inds=vecnorm(contour_lines(loop_inds).uv-mesh_boundaries);
% [~,min_ind]=min(point_inds);
% contour_lines(loop_inds).uv=circshift(contour_lines(loop_inds).uv,min_ind*(-1)+1,2);
% contour_lines(loop_inds).point_coordinates=circshift(contour_lines(loop_inds).point_coordinates,min_ind*(-1)+1,2);
% end
% end

    
% figure;
% hold on;
% for group_ind=1:numel(group_container)
% high_cut_p=[high_cuts(group_ind).v(:).coords];
% low_cut_p=[low_cuts(group_ind).v(:).coords];
% for loop_ind=1:numel(group_container(group_ind).loops)
% plot3(group_container(group_ind).loops(loop_ind).point_coordinates(1,:),...
% group_container(group_ind).loops(loop_ind).point_coordinates(2,:),...
% group_container(group_ind).loops(loop_ind).point_coordinates(3,:));
% end
% plot3(high_cut_p(1,:),high_cut_p(2,:),high_cut_p(3,:),'r');
% plot3(low_cut_p(1,:),low_cut_p(2,:),low_cut_p(3,:),'g');
% end
% hold off;


% hold on; 
% for group_ind=1:numel(group_container)
% for loop_ind=1:numel(group_container(group_ind).loops)
%  scatter3(      group_container(group_ind).loops(loop_ind).point_coordinates(1,group_container(group_ind).loops(loop_ind).highest_lowest_point_ind),...
%                     group_container(group_ind).loops(loop_ind).point_coordinates(2,group_container(group_ind).loops(loop_ind).highest_lowest_point_ind),...
%                     group_container(group_ind).loops(loop_ind).point_coordinates(3,group_container(group_ind).loops(loop_ind).highest_lowest_point_ind),'g','filled');
% end
% end
% hold off;

% group_container(group_ind).is_perpendicular_to_b0=zeros(1,numel(group_container(group_ind).loops));
% for loop_ind=1:numel(group_container(group_ind).loops)
% loop_points=group_container(group_ind).loops(loop_ind).point_coordinates;
% rotated_points=rotationMatrix*loop_points;
% b0_projection_height=rotated_points(1,:);
% %check if all loop points have simular projection heigth which means
% %that the loop is perpendicular to the b0_direction
% aprox_loop_radius=mean(vecnorm(group_container(group_ind).loops(loop_ind).point_coordinates-group_center.v(:,group_ind)));
% if std(b0_projection_height)/aprox_loop_radius<0.1
% group_container(group_ind).loops(loop_ind).highest_lowest_point_ind=[nan nan];
% group_container(group_ind).is_perpendicular_to_b0(loop_ind)=1;
% else
% [~,b0_projection_height_sort_inds]=sort(b0_projection_height);
% group_container(group_ind).loops(loop_ind).highest_lowest_point_ind=[b0_projection_height_sort_inds(1) b0_projection_height_sort_inds(end)];
% end
% end

%define the orientation of the cut plane
%find vector orthorgonal to B0 and the normal of vector of mesh at the
%group center

% if abs(dot(loop_normal,b_0_direction))>0.9
% [~,vi_max_Ind]=max(v1);
% loop_normal=v1.*point_vol_size(1,vi_max_Ind);
% end

% rot_angle=acos(dot(v1, v2) / (norm(v1) * norm(v2)));
% rotationVector = cross(v1,v2);
% rotationVector=rotationVector./norm(rotationVector).*rot_angle;
% rotationMatrix = rotationVectorToMatrix(rotationVector);

end

