function cut_position=find_group_cut_position(loop_group,group_center,mesh,b_0_direction,cut_plane_definition)

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
cut_position(numel(loop_group.loops)).cut_point=[];
cut_position(numel(loop_group.loops)).cut_point=[];
cut_position(numel(loop_group.loops)).high_cut=[];
cut_position(numel(loop_group.loops)).low_cut=[];
cut_position(numel(loop_group.loops)).high_cutshape=[];
cut_position(numel(loop_group.loops)).low_cutshape=[];
cut_position(numel(loop_group.loops)).cut_widh_high=[];
cut_position(numel(loop_group.loops)).cut_widh_low=[];


for loop_ind=1:numel(loop_group.loops)


cut_position(loop_ind).cut_point.v=[];
cut_position(loop_ind).cut_point.uv=[];
cut_position(loop_ind).cut_point.segment_ind=[];
for point_ind=2:size(loop_group.loops(loop_ind).v,2)
point_a=loop_group.loops(loop_ind).v(:,point_ind-1)';
point_b=loop_group.loops(loop_ind).v(:,point_ind)';
[cut_p,cut_flag]=plane_line_intersect(cut_plane_direction,group_center',point_a,point_b);
if cut_flag==1 %test if there is a cut between the line points
cut_position(loop_ind).cut_point.v=[cut_position(loop_ind).cut_point.v cut_p'];
cut_position(loop_ind).cut_point.segment_ind=[cut_position(loop_ind).cut_point.segment_ind point_ind-1];
%build the corresponding uv point
cut_point_ratio=norm(point_a-cut_p)/norm(point_a-point_b);
point_a_uv=loop_group.loops(loop_ind).uv(:,point_ind-1);
point_b_uv=loop_group.loops(loop_ind).uv(:,point_ind);
cut_poin_uv=point_a_uv+(point_b_uv- point_a_uv)*cut_point_ratio;
cut_position(loop_ind).cut_point.uv=[cut_position(loop_ind).cut_point.uv cut_poin_uv];
end
end

%delete repeating degenerate cut points 
is_repeating_cutpoint=all(abs([[1 1]' diff(cut_position(loop_ind).cut_point.uv,1,2)])<10^(-10),1);
cut_position(loop_ind).cut_point.v(:,is_repeating_cutpoint)=[];
cut_position(loop_ind).cut_point.uv(:,is_repeating_cutpoint)=[];
cut_position(loop_ind).cut_point.segment_ind(is_repeating_cutpoint)=[];


%Take care of the exception that the cut plane has not intersection with
%the loop; In this case, apply alterative calcuation of cut_position by projection:
if isempty(cut_position(loop_ind).cut_point.v)
[~,max_ind]=max(loop_group.loops(loop_ind).v(1,:)*cut_plane_direction(1)+loop_group.loops(loop_ind).v(2,:)*cut_plane_direction(2)+loop_group.loops(loop_ind).v(3,:)*cut_plane_direction(3));
[~,min_ind]=min(loop_group.loops(loop_ind).v(1,:)*cut_plane_direction(1)+loop_group.loops(loop_ind).v(2,:)*cut_plane_direction(2)+loop_group.loops(loop_ind).v(3,:)*cut_plane_direction(3));
cut_position(loop_ind).cut_point.v=[loop_group.loops(loop_ind).v(:,min_ind) loop_group.loops(loop_ind).v(:,max_ind)];
cut_position(loop_ind).cut_point.uv=[loop_group.loops(loop_ind).uv(:,min_ind) loop_group.loops(loop_ind).uv(:,max_ind)];
if max_ind~=size(loop_group.loops(loop_ind).v,2)
cut_position(loop_ind).cut_point.segment_ind=[min_ind max_ind-1];
else
cut_position(loop_ind).cut_point.segment_ind=[min_ind max_ind];
end
end


%seperated in higher and lower cut points:
%first: use the 2d representation of the loop 

%Since cutpoints are always alternating between cutpoints "in" and "out" ,
%(seen from the orientation of the cut plane); 

if loop_ind == 1


%Select the first pair of the cuts, the cuts with largest distance to group
%center
[~,cut_sort_ind] = sort(vecnorm(cut_position(loop_ind).cut_point.v-group_center));
first_pair=[cut_sort_ind(1) cut_sort_ind(2)]; 

%find the direction for which high and low cuts are seperated
cut_direction=cut_position(loop_ind).cut_point.v(:,first_pair(2))-cut_position(loop_ind).cut_point.v(:,first_pair(1));
cut_direction=cut_direction./vecnorm(cut_direction);
if sum(b_0_direction.*cut_direction)<0
    cut_direction=cut_direction.*(-1);
end


%project the coordinates of the cut pairs
[~,min_ind]=min(sum(cut_position(loop_ind).cut_point.v(:,first_pair).*cut_direction,1));
high_ind=first_pair(min_ind);
low_ind=first_pair(first_pair~=high_ind);

high_cut_primer=cut_position(loop_ind).cut_point.v(:,high_ind);
low_cut_primer=cut_position(loop_ind).cut_point.v(:,low_ind);


cut_position(loop_ind).high_cut.segment_ind=cut_position(loop_ind).cut_point.segment_ind(high_ind);
cut_position(loop_ind).low_cut.segment_ind=cut_position(loop_ind).cut_point.segment_ind(low_ind);
cut_position(loop_ind).high_cut.v=cut_position(loop_ind).cut_point.v(:,high_ind);
cut_position(loop_ind).low_cut.v=cut_position(loop_ind).cut_point.v(:,low_ind);
cut_position(loop_ind).high_cut.uv=cut_position(loop_ind).cut_point.uv(:,high_ind);
cut_position(loop_ind).low_cut.uv=cut_position(loop_ind).cut_point.uv(:,low_ind);
center_first_cut=(cut_position(loop_ind).high_cut.v+cut_position(loop_ind).low_cut.v)./2;

else %now for the follwing inner loops


%Choose the following cut pairs regarding their distance to the high and
%low cut of the first loop
high_dists=vecnorm(cut_position(loop_ind).cut_point.v-high_cut_primer);
low_dists=vecnorm(cut_position(loop_ind).cut_point.v-low_cut_primer);

[~,high_min_ind]=min(high_dists);
[~,low_min_ind]=min(low_dists);

cut_position(loop_ind).high_cut.segment_ind=cut_position(loop_ind).cut_point.segment_ind(high_min_ind);
cut_position(loop_ind).low_cut.segment_ind=cut_position(loop_ind).cut_point.segment_ind(low_min_ind);
cut_position(loop_ind).high_cut.v=cut_position(loop_ind).cut_point.v(:,high_min_ind);
cut_position(loop_ind).low_cut.v=cut_position(loop_ind).cut_point.v(:,low_min_ind);
cut_position(loop_ind).high_cut.uv=cut_position(loop_ind).cut_point.uv(:,high_min_ind);
cut_position(loop_ind).low_cut.uv=cut_position(loop_ind).cut_point.uv(:,low_min_ind);


end



end


end