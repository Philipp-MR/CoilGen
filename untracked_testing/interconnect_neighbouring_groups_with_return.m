function [fused_group,opening_cut]=interconnect_neighbouring_groups_with_return(group_1,group_2,opening_gap,return_points_group2,parameterized_mesh)
%find the cut rectangle using a return path of one of the groups


% find a new cut point by using the return path of the outer loop
return_points_group2=fliplr(return_points_group2);


return_points_check_points=return_points_group2;
[~,direction_point_ind1,~]=find_shortest_path_between_loops_2d(return_points_check_points,group_1.uv);
direction_point_1=return_points_check_points(:,direction_point_ind1);
return_points_check_points=return_points_check_points(:,setdiff(1:size(return_points_group2,2),direction_point_ind1));
[~,direction_point_ind2,~]=find_shortest_path_between_loops_2d(return_points_check_points,group_1.uv);
direction_point_2=return_points_check_points(:,direction_point_ind2);
cut_direction=direction_point_2-direction_point_1;
cut_direction=cut_direction./vecnorm(cut_direction)*10^7;
[cutx,cuty] = polyxpoly([ direction_point_1(1)+cut_direction(1) direction_point_1(1)-cut_direction(1)],...
[direction_point_1(2)+cut_direction(2) direction_point_1(2)-cut_direction(2)],...
group_1.uv(1,[1:end 1]),group_1.uv(2,[1:end 1]));
if ~isempty(cutx)
[~,min_ind,~]=find_shortest_path_between_loops_2d([cutx'; cuty'],return_points_group2);
cut_point_group1=[cutx(min_ind);cuty(min_ind)];
[~,cut_point_group2]=min(vecnorm(return_points_group2-cut_point_group1));
cut_point_group2=return_points_group2(:,cut_point_group2);
else
[~,~,direction_point_ind2]=find_shortest_path_between_loops_2d(return_points_check_points,group_1.uv);
direction_point_2=group_1.uv(:,direction_point_ind2);
cut_direction=direction_point_2-direction_point_1;
cut_direction=cut_direction./vecnorm(cut_direction)*10^7;
[cutx,cuty] = polyxpoly([ direction_point_1(1)+cut_direction(1) direction_point_1(1)-cut_direction(1)],...
[direction_point_1(2)+cut_direction(2) direction_point_1(2)-cut_direction(2)],...
group_1.uv(1,[1:end 1]),group_1.uv(2,[1:end 1]));
[~,min_ind,~]=find_shortest_path_between_loops_2d([cutx'; cuty'],return_points_group2);
cut_point_group1=[cutx(min_ind);cuty(min_ind)];
[~,cut_point_group2]=min(vecnorm(return_points_group2-cut_point_group1));
cut_point_group2=return_points_group2(:,cut_point_group2);
end

longitudinal_vector=(cut_point_group1-cut_point_group2).*1.2;
othorgonal_vector=[longitudinal_vector(2);longitudinal_vector(1).*(-1)];
othorgonal_vector=othorgonal_vector./norm(othorgonal_vector);
middle_point=(cut_point_group1+cut_point_group2)./2;
local_opening_gab=calc_local_opening_gab(parameterized_mesh,middle_point,othorgonal_vector,opening_gap,10);

othorgonal_vector=othorgonal_vector./(norm(othorgonal_vector)).*local_opening_gab;
opening_cut =[middle_point+longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2+othorgonal_vector/2];

% %create a circular opening cut
% opening_cut=[sin(0:(2*pi)/(20-1):2*pi); cos(0:(2*pi)/(20-1):2*pi)].*(local_opening_gab_a/2);
% opening_cut=opening_cut+(cut_point_group1+cut_point_group2)./2;

%circshift the inner group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_1.uv-cut_point_group1));
group_1.uv=circshift(group_1.uv,min_ind*(-1),2);

%circshift the outer group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_2.uv-cut_point_group2));
group_2.uv=circshift(group_2.uv,min_ind*(-1),2);

%fuse both groups
%Check which fusing order is better:
track_combilength1=[group_1.uv group_2.uv];
track_combilength2=[group_2.uv group_1.uv];
track_combilength1=sum(vecnorm(track_combilength1(:,2:end)-track_combilength1(:,1:end-1)));
track_combilength2=sum(vecnorm(track_combilength2(:,2:end)-track_combilength2(:,1:end-1)));
if track_combilength1<track_combilength2
fused_group=[group_1.uv group_2.uv];
else
fused_group=[group_2.uv group_1.uv];
end

%Clear the points witihn the opening_cut to improve to connection
[in_cut_rect_ind,~] = inpolygon(fused_group(1,:),fused_group(2,:),opening_cut(1,:),opening_cut(2,:)); 
fused_group(:,in_cut_rect_ind)=[];

end