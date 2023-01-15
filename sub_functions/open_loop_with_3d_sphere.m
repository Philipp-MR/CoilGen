function [opened_loop,uv_cut,cut_points]=open_loop_with_3d_sphere(curve_points_in,sphere_center,sphere_diameter)
%Opening a loop with by overlapping it with a 3D sphere with given radia
%and given center position
%@Philipp Amrein, Uniklinik Freiburg,2023


%remove doubled points from the curve
points_to_delete=vecnorm(curve_points_in.v(:,[2:end 1])-curve_points_in.v(:,[1:end]))<10^(-10);
curve_points_in.v(:,points_to_delete)=[];
curve_points_in.uv(:,points_to_delete)=[];
curve_points_in.number_points=size(curve_points_in.v,2);

%Add a point within the curve which has the shortest distance to the sphere
[curve_points,~]=add_nearest_ref_point_to_curve(curve_points_in,sphere_center);

inside_sphere_ind=vecnorm(curve_points.v-sphere_center)<sphere_diameter/2;
%circshift the path so that the starting location is outside the sphere
if any(inside_sphere_ind)&~all(inside_sphere_ind) %throw error if no overlapp with the sphere
min_ind=min(find(~inside_sphere_ind))-1;
curve_points.v=circshift(curve_points.v,-(min_ind),2);
curve_points.uv=circshift(curve_points.uv,-(min_ind),2);
inside_sphere_ind=vecnorm(curve_points.v-sphere_center)<sphere_diameter/2;
if inside_sphere_ind(end) %shift again to avid problems at the end of the curve
curve_points.v=circshift(curve_points.v,-1,2);
curve_points.uv=circshift(curve_points.uv,-1,2);
end
else
error('Opening of loop not possible, no overlapp between cut sphere and loop')
end

inside_sphere_ind=vecnorm(curve_points.v-sphere_center)<sphere_diameter/2;

%In case of multiple cuts with the sphere select the part of the curve
%which is closer to the sphere center
if sum(abs(diff(inside_sphere_ind)))>2 %when there are multiple parts..
parts_start=find(diff(inside_sphere_ind)==1)+1;
parts_end=find(diff(inside_sphere_ind)==-1);
parts_avg_dist=zeros(1,numel(parts_start));
for part_ind=1:numel(parts_start)
parts_avg_dist(part_ind)=mean(vecnorm(curve_points.v(:,parts_start(part_ind):parts_end(part_ind))-sphere_center));
end
[~,nearest_part]=min(parts_avg_dist);
inside_sphere_ind_unique=false(size(inside_sphere_ind));
inside_sphere_ind_unique(parts_start(nearest_part):parts_end(nearest_part))=true;
else
inside_sphere_ind_unique=inside_sphere_ind;
end


%Choose the bigger part(old)
% [~,bigger_over_part]=max(find(diff(inside_sphere_ind)==-1)-find(diff(inside_sphere_ind)==1));
% part_start=find(diff(inside_sphere_ind)==1);
% part_start=part_start(bigger_over_part)+1;
% part_end=find(diff(inside_sphere_ind)==-1);
% part_end=part_end(bigger_over_part);
% inside_sphere_ind=false(1,size(curve_points.v,2));
% inside_sphere_ind(part_start:part_end)=true;

%Find the positions where the curve enters the sphere
first_sphere_penetration_locations=find(abs(diff(inside_sphere_ind_unique)));
second_sphere_penetration_locations=find(abs(diff(inside_sphere_ind_unique)))+1;

first_distances=vecnorm(curve_points.v(:,first_sphere_penetration_locations)-sphere_center);
second_distances=vecnorm(curve_points.v(:,second_sphere_penetration_locations)-sphere_center);
% first_distances=vecnorm(curve_points.v(:,first_sphere_penetration_locations)-sphere_center)./(sphere_diameter/2);
% second_distances=vecnorm(curve_points.v(:,second_sphere_penetration_locations)-sphere_center)./(sphere_diameter/2);


sphrere_crossing_vecs.v=curve_points.v(:,second_sphere_penetration_locations)-curve_points.v(:,first_sphere_penetration_locations);
sphrere_crossing_vecs.uv=curve_points.uv(:,second_sphere_penetration_locations)-curve_points.uv(:,first_sphere_penetration_locations);


%Caluclate the penetration points by means of interpolation of weighted mean for the radial
%distance
repeated_radia=ones(1,numel(first_distances))*sphere_diameter/2;
cut_points.v=curve_points.v(:,first_sphere_penetration_locations)+sphrere_crossing_vecs.v.*repmat((repeated_radia-first_distances)./(second_distances-first_distances),[3 1]);
cut_points.uv=curve_points.uv(:,first_sphere_penetration_locations)+sphrere_crossing_vecs.uv.*repmat((repeated_radia-first_distances)./(second_distances-first_distances),[2 1]);


%Open the loop; Check which parts of the curve are inside or outside the sphere
%inside_sphere_ind=vecnorm(curve_points.v-sphere_center)<sphere_diameter/2;
shift_ind=(min(find(inside_sphere_ind_unique==1))-1)*(-1);
curve_points.v=circshift(curve_points.v,shift_ind,2);
curve_points.uv=circshift(curve_points.uv,shift_ind,2);
inside_sphere_ind_unique=circshift(inside_sphere_ind_unique,shift_ind);
curve_points.v=curve_points.v(:,~inside_sphere_ind_unique);
curve_points.uv=curve_points.uv(:,~inside_sphere_ind_unique);

%Build the "opened" loop with the cut_points as open ends
%Remove curve points which are still inside the sphere
finished_loop_case1.v=[cut_points.v(:,1) curve_points.v cut_points.v(:,end)];
finished_loop_case2.v=[cut_points.v(:,end) curve_points.v cut_points.v(:,1)];
finished_loop_case1.uv=[cut_points.uv(:,1) curve_points.uv cut_points.uv(:,end)];
finished_loop_case2.uv=[cut_points.uv(:,end) curve_points.uv cut_points.uv(:,1)];
mean_dist_1=sum(vecnorm(finished_loop_case1.v(:,2:end)-finished_loop_case1.v(:,1:end-1)));
mean_dist_2=sum(vecnorm(finished_loop_case2.v(:,2:end)-finished_loop_case2.v(:,1:end-1)));
if mean_dist_1<mean_dist_2
opened_loop.v=finished_loop_case1.v;
opened_loop.uv=finished_loop_case1.uv;
else
opened_loop.v=finished_loop_case2.v;
opened_loop.uv=finished_loop_case2.uv;
end

%Generate the 2d contour of the cut shape for later plotting
radius_2d=vecnorm(opened_loop.uv(:,1)-opened_loop.uv(:,end))/2;
uv_cut=[sin([0:50]./(50/(2*pi)));cos([0:50]./(50/(2*pi)))].*radius_2d+(opened_loop.uv(:,1)+opened_loop.uv(:,end))./2;


function [curve_track_out,near_points]=add_nearest_ref_point_to_curve(curve_track_in,target_point)
%Calculate the mutal nearest positions and segment indices between two
%loops 
% Copyright: Philipp Amrein, Uniklinik Freiburg Mai 2022
curve_track=curve_track_in;
if curve_track.v(1,1)~=curve_track.v(1,end)|curve_track.v(2,1)~=curve_track.v(2,end)|curve_track.v(3,1)~=curve_track.v(3,end)
curve_track.v=[curve_track.v curve_track.v(:,1)];
curve_track.uv=[curve_track.uv curve_track.uv(:,1)];
end
seg_starts.v=curve_track.v(:,1:end-1);
seg_starts.uv=curve_track.uv(:,1:end-1);
seg_ends.v=curve_track.v(:,2:end);
seg_ends.uv=curve_track.uv(:,2:end);
t=dot(target_point-seg_starts.v,seg_ends.v-seg_starts.v)./sum((seg_ends.v-seg_starts.v).^2,1);
t(t<0)=0;
t(t>1)=1;
all_near_points.v=seg_starts.v+(seg_ends.v-seg_starts.v).*repmat(t,[3 1]);
all_near_points.uv=seg_starts.uv+(seg_ends.uv-seg_starts.uv).*repmat(t,[2 1]);
all_dists=vecnorm(all_near_points.v-target_point);
[~,min_ind_seq]=min(all_dists);
near_points.v=all_near_points.v(:,min_ind_seq);
near_points.uv=all_near_points.uv(:,min_ind_seq);
%Add the near point within the curve
curve_track_out.v=[curve_track_in.v(:,1:min_ind_seq) near_points.v];
curve_track_out.uv=[curve_track_in.uv(:,1:min_ind_seq) near_points.uv];
if min_ind_seq~=size(curve_track_in.v,2)
curve_track_out.v=[curve_track_out.v curve_track_in.v(:,min_ind_seq+1:end)];
curve_track_out.uv=[curve_track_out.uv curve_track_in.uv(:,min_ind_seq+1:end)];
end
end


end