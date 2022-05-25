function [fused_group,opening_cut]=interconect_neighbouring_groups(group_1,group_2,cutpoint_1,cutpoint_2,opening_gap,parameterized_mesh,winding_distance,raw_cut_point_dist)
%do simple connection



longitudinal_vector=cutpoint_1-cutpoint_2;
othorgonal_vector=[longitudinal_vector(2);longitudinal_vector(1).*(-1)];
othorgonal_vector=othorgonal_vector./norm(othorgonal_vector);
middle_point=(cutpoint_1+cutpoint_2)./2;
local_opening_gab=calc_local_opening_gab(parameterized_mesh,middle_point,othorgonal_vector,opening_gap,10);

% %create a circular opening cut
% opening_cut=[sin(0:(2*pi)/(20-1):2*pi); cos(0:(2*pi)/(20-1):2*pi)].*(local_opening_gab/2);
% opening_cut=opening_cut+middle_point;


longitudinal_vector=longitudinal_vector.*1.5;
othorgonal_vector=othorgonal_vector./(norm(othorgonal_vector)).*local_opening_gab;


%Check if a geodesic cut shape is necessary (in  case of bigger distances between the cut points)

opening_cut =[middle_point+longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2-othorgonal_vector/2 ...
                                middle_point-longitudinal_vector/2+othorgonal_vector/2 ...
                                middle_point+longitudinal_vector/2+othorgonal_vector/2];

% %create a circular opening cut
% opening_cut=[sin(0:(2*pi)/(20-1):2*pi); cos(0:(2*pi)/(20-1):2*pi)].*(local_opening_gab/2);
% opening_cut=opening_cut+middle_point;      
                            
                            
%circshift the inner group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_1.uv-cutpoint_1));
group_1.uv=circshift(group_1.uv,min_ind*(-1),2);

%circshift the outer group that the open ends are the return path
[~,min_ind]=min(vecnorm(group_2.uv-cutpoint_2));
group_2.uv=circshift(group_2.uv,min_ind*(-1),2);



%     %add the "cut"points for a clean cut
%     [cut_x,cut_y] = polyxpoly(cut_points(1,:),cut_points(2,:),rectangle_cut(1,:),rectangle_cut(2,:));
%     opened_group{loop_num}=raw_loop;
%     %find out the order in which the cut points must inlcuded in the opened loop:
%     test_dists_1=vecnorm(raw_loop(:,1)-[cut_x(1);cut_y(1)]);
%     test_dists_2=vecnorm(raw_loop(:,1)-[cut_x(2);cut_y(2)]);
%     if    test_dists_1<test_dists_2
%         opened_group{loop_num}=[[cut_x(1) cut_y(1)]' opened_group{loop_num} [cut_x(2) cut_y(2)]' ] ;
%     else 
%         opened_group{loop_num}=[[cut_x(2) cut_y(2)]' opened_group{loop_num} [cut_x(1) cut_y(1)]' ] ;
%     end


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