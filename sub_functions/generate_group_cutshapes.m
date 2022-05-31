function cut_shapes=generate_group_cutshapes(loop_group,group_center,mesh,cut_width,b_0_direction,cut_plane_definition,cut_height_ratio)

%define the cut plane orientation for the group
%find the opening shapes and cut points
%seperated in higher and lower cut points
%check wether a forced cut seletion is given
%generate cirular cut shape and delte overlapping points

loop_normal=calc_mean_loop_normal(loop_group,mesh);

%test if loop normal and b0_direction are not independent enough
if norm(cross(b_0_direction./norm(b_0_direction),loop_normal))<0.01

if strcmp('nearest',cut_plane_definition)
[~,min_ind]=min(vecnorm(loop_group.loops(1).v-group_center));
alternative_cut_direction=loop_group.loops(1).v(:,min_ind)-group_center;
alternative_cut_direction=alternative_cut_direction./(vecnorm(alternative_cut_direction));
cut_plane_direction=cross(alternative_cut_direction,loop_normal);
cut_plane_direction=cut_plane_direction./norm(cut_plane_direction);
else
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
cut_plane_direction=cross(b_0_direction,alternative_center_normal);
cut_plane_direction=cut_plane_direction./norm(cut_plane_direction);

end

else
cut_plane_direction=cross(b_0_direction,loop_normal);
cut_plane_direction=cut_plane_direction./norm(cut_plane_direction);
end


%Calcualte the cut shapes
%find the cut points of the cut plane with the loops
cut_shapes(numel(loop_group.loops)).cut_point=[];
cut_shapes(numel(loop_group.loops)).cut_point=[];
cut_shapes(numel(loop_group.loops)).high_cut=[];
cut_shapes(numel(loop_group.loops)).low_cut=[];
cut_shapes(numel(loop_group.loops)).high_cutshape=[];
cut_shapes(numel(loop_group.loops)).low_cutshape=[];
cut_shapes(numel(loop_group.loops)).cut_widh_high=[];
cut_shapes(numel(loop_group.loops)).cut_widh_low=[];


for loop_ind=1:numel(loop_group.loops)


cut_shapes(loop_ind).cut_point.v=[];
cut_shapes(loop_ind).cut_point.uv=[];
cut_shapes(loop_ind).cut_point.segment_ind=[];
for point_ind=2:size(loop_group.loops(loop_ind).v,2)
point_a=loop_group.loops(loop_ind).v(:,point_ind-1)';
point_b=loop_group.loops(loop_ind).v(:,point_ind)';
[cut_p,cut_flag]=plane_line_intersect(cut_plane_direction,group_center',point_a,point_b);
if cut_flag==1 %test if there is a cut between the line points
cut_shapes(loop_ind).cut_point.v=[cut_shapes(loop_ind).cut_point.v cut_p'];
cut_shapes(loop_ind).cut_point.segment_ind=[cut_shapes(loop_ind).cut_point.segment_ind point_ind-1];
%build the corresponding uv point
cut_point_ratio=norm(point_a-cut_p)/norm(point_a-point_b);
point_a_uv=loop_group.loops(loop_ind).uv(:,point_ind-1);
point_b_uv=loop_group.loops(loop_ind).uv(:,point_ind);
cut_poin_uv=point_a_uv+(point_b_uv- point_a_uv)*cut_point_ratio;
cut_shapes(loop_ind).cut_point.uv=[cut_shapes(loop_ind).cut_point.uv cut_poin_uv];
end
end

%delete repeating degenerate cut points 
is_repeating_cutpoint=all(abs([[1 1]' diff(cut_shapes(loop_ind).cut_point.uv,1,2)])<10^(-10),1);
cut_shapes(loop_ind).cut_point.v(:,is_repeating_cutpoint)=[];
cut_shapes(loop_ind).cut_point.uv(:,is_repeating_cutpoint)=[];
cut_shapes(loop_ind).cut_point.segment_ind(is_repeating_cutpoint)=[];

%seperated in higher and lower cut points:
%first: use the 2d representation of the loop 

%Since cutpoints are always alternating between cutpoints "in" and "out" ,
%(seen from the orientation of the cut plane); 


%choose for the first pair of cuts the pair for which the combinded distance to the geometric center of the group is the
%smallest
if loop_ind == 1

pair_inds=repelem(1:size(cut_shapes(loop_ind).cut_point.v,2)/2,2); %create indices of the cut pairs

first_pair_dists_to_center=vecnorm(cut_shapes(loop_ind).cut_point.v-group_center);
first_pair_dists_to_center=first_pair_dists_to_center(1:end-1)+first_pair_dists_to_center(2:end);
first_pair_dists_to_center=first_pair_dists_to_center(logical(mod(1:numel(first_pair_dists_to_center),2)));
[~,first_pair_dists_to_center]=min(first_pair_dists_to_center);
first_pair=find(pair_inds==first_pair_dists_to_center);

%find the direction for which high and low cuts are seperated
cut_direction=cut_shapes(loop_ind).cut_point.v(:,first_pair(2))-cut_shapes(loop_ind).cut_point.v(:,first_pair(1));
cut_direction=cut_direction./vecnorm(cut_direction);

%project the coordinates of the cut pairs
[~,min_ind]=min(sum(cut_shapes(loop_ind).cut_point.v(:,first_pair).*cut_direction,1));

cut_shapes(loop_ind).high_cut.segment_ind=cut_shapes(loop_ind).cut_point.segment_ind(first_pair(min_ind));
cut_shapes(loop_ind).low_cut.segment_ind=cut_shapes(loop_ind).cut_point.segment_ind(first_pair([1 2]~=min_ind));
cut_shapes(loop_ind).high_cut.v=cut_shapes(loop_ind).cut_point.v(:,first_pair(min_ind));
cut_shapes(loop_ind).low_cut.v=cut_shapes(loop_ind).cut_point.v(:,first_pair([1 2]~=min_ind));
cut_shapes(loop_ind).high_cut.uv=cut_shapes(loop_ind).cut_point.uv(:,first_pair(min_ind));
cut_shapes(loop_ind).low_cut.uv=cut_shapes(loop_ind).cut_point.uv(:,first_pair([1 2]~=min_ind));

center_first_cut=(cut_shapes(loop_ind).high_cut.v+cut_shapes(loop_ind).low_cut.v)./2;

else
%choose for every following loop the pair of interconnections for which has the smallest distance to the previous pair
pair_inds=repelem(1:size(cut_shapes(loop_ind).cut_point.v,2)/2,2); %create indices of the cut pairs

dists_to_center=vecnorm(cut_shapes(loop_ind).cut_point.v-center_first_cut);
dists_to_center=dists_to_center(1:end-1)+dists_to_center(2:end);
dists_to_center=dists_to_center(logical(mod(1:numel(dists_to_center),2)));
[~,min_ind]=min(dists_to_center);
next_pair=find(pair_inds==min_ind);

%project the coordinates of the cut pairs
[~,min_ind]=min(sum(cut_shapes(loop_ind).cut_point.v(:,first_pair).*cut_direction,1));
cut_shapes(loop_ind).high_cut.segment_ind=cut_shapes(loop_ind).cut_point.segment_ind(first_pair(min_ind));
cut_shapes(loop_ind).low_cut.segment_ind=cut_shapes(loop_ind).cut_point.segment_ind(first_pair([1 2]~=min_ind));
cut_shapes(loop_ind).high_cut.v=cut_shapes(loop_ind).cut_point.v(:,first_pair(min_ind));
cut_shapes(loop_ind).low_cut.v=cut_shapes(loop_ind).cut_point.v(:,first_pair([1 2]~=min_ind));
cut_shapes(loop_ind).high_cut.uv=cut_shapes(loop_ind).cut_point.uv(:,first_pair(min_ind));
cut_shapes(loop_ind).low_cut.uv=cut_shapes(loop_ind).cut_point.uv(:,first_pair([1 2]~=min_ind));


end

%finaly build cut shapes for the high and low position
cut_shapes(loop_ind).cut_widh_high=calc_local_opening_gab(loop_group.loops(loop_ind),cut_shapes(loop_ind).high_cut.segment_ind+1,cut_shapes(loop_ind).high_cut.segment_ind,cut_width);
cut_shapes(loop_ind).cut_widh_low=calc_local_opening_gab(loop_group.loops(loop_ind),cut_shapes(loop_ind).high_cut.segment_ind+1,cut_shapes(loop_ind).high_cut.segment_ind,cut_width);

cut_shapes(loop_ind).high_cutshape=build_cut_rectangle(loop_group.loops(loop_ind),cut_shapes(loop_ind).high_cut.uv,cut_shapes(loop_ind).high_cut.segment_ind,cut_shapes(loop_ind).cut_widh_high,cut_height_ratio);
cut_shapes(loop_ind).low_cutshape=build_cut_rectangle(loop_group.loops(loop_ind),cut_shapes(loop_ind).low_cut.uv,cut_shapes(loop_ind).low_cut.segment_ind,cut_shapes(loop_ind).cut_widh_low,cut_height_ratio);



end


end