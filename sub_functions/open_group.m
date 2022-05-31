function opened_group=open_group(group_to_open,cut_shape)
%open the groups with the cut rectangle and interconnect them

raw_group=group_to_open.uv;


%find the group indice most near to the center of the cut_shape
%average point positions together with their neightbours




cut_center=mean(cut_shape,2);

folded_dists=   vecnorm(group_to_open.uv(:,circshift(1:size(group_to_open.uv,2),2))-cut_center)+...
                        vecnorm(group_to_open.uv(:,circshift(1:size(group_to_open.uv,2),1))-cut_center)+...
                        vecnorm(group_to_open.uv-cut_center)+...
                        vecnorm(group_to_open.uv(:,circshift(1:size(group_to_open.uv,2),1))-cut_center)+...
                        vecnorm(group_to_open.uv(:,circshift(1:size(group_to_open.uv,2),-2))-cut_center);

[~,nearest_to_cut_center_inds]=min(folded_dists);




%close the raw group
if raw_group(1,1)~=raw_group(1,end)&raw_group(2,1)~=raw_group(2,end)
raw_group=[raw_group raw_group(:,end)]; %close the raw group
end


%calculate the intersection points with the respective group segments
intersection_points=find_segment_intersections(raw_group,cut_shape); 
all_intersec_segments=arrayfun(@(x) intersection_points(x).segment_inds,1:numel(intersection_points),'UniformOutput',0);
all_intersec_positions=[intersection_points(:).uv];
num_intersec_points=sum(~isnan([all_intersec_segments{:}]));
[sorted_segments,sort_ind]=sort([all_intersec_segments{:}]);
sorted_positions=all_intersec_positions(:,sort_ind);
sorted_positions=sorted_positions(:,~isnan(sorted_segments));
sorted_segments=sorted_segments(~isnan(sorted_segments));
% cuts_per_side=arrayfun(@(x) all(~isnan(intersection_points(x).segment_inds)).*numel(intersection_points(x).segment_inds),1:numel(intersection_points));
% any_cut_flag=~isempty([intersection_points(:).uv]);

%Build cuts between the group and the cutshape
if num_intersec_points>0% Check wether there are any intersections


%check wether the group start&end are outside or inside the cutshape
test_intersection_points=find_segment_intersections(cut_shape,[cut_center raw_group(:,1)]);
if isnan(test_intersection_points(1).segment_inds)
    ends_within_cutshape=true;
else
    ends_within_cutshape=false;
end

%shift the cut pairs for 1; shift the order of the cutted segments for one position so that we have proper pairs
if ends_within_cutshape
sorted_segments=circshift(sorted_segments,-1);
sorted_positions=circshift(sorted_positions,-1,2);
end


ind_dists_to_cutshape=abs(sorted_segments-nearest_to_cut_center_inds);
pair_ind_dists_to_cutshape=zeros(1,numel(sorted_segments)/2);
for pair_ind=1:numel(sorted_segments)/2
pair_ind_dists_to_cutshape(pair_ind)=sum(ind_dists_to_cutshape(repelem(1:numel(sorted_segments)/2,2)==pair_ind));
end
[~,min_ind]=min(pair_ind_dists_to_cutshape);
final_pair_inds=find(repelem(1:numel(sorted_segments)/2,2)==min_ind);
final_segment_pair=sorted_segments(final_pair_inds);

%[final_segment_pair,sort_ind]=sort(final_segment_pair);
final_segment_positions=sorted_positions(:,final_pair_inds);
first_cut_position=final_segment_positions(:,1);
last_cut_position=final_segment_positions(:,2);
in_cut_rect_inds=(final_segment_pair(1)+1):final_segment_pair(2);


%delete the points within the cutshape
raw_group(:,in_cut_rect_inds)=[];

%shift the opening around the cutshape
if ~ends_within_cutshape & ~isempty(in_cut_rect_inds)
raw_group = circshift(raw_group,(in_cut_rect_inds(1)-1)*(-1),2);
% else %point indice overflow
% raw_group = circshift(raw_group,(in_cut_rect_inds(end)-1)*(-1),2);
elseif ~ends_within_cutshape & isempty(in_cut_rect_inds) % take care of the possibilty of no points within the cutshape
raw_group = circshift(raw_group,(sorted_segments(1))*(-1),2);
end

%add end points for a clean proper cut
%check which cutpoint belongs to which open end of the group
track1=[first_cut_position raw_group last_cut_position];
track2=[last_cut_position raw_group first_cut_position];
if sum(vecnorm(track1(:,2:end)-track1(:,1:end-1)))<sum(vecnorm(track2(:,2:end)-track2(:,1:end-1)))
opened_group=track1;
else
opened_group=track2;
end

else
error('Degenerate cut: Cutshape does not overlapp with the group; change the interconnection width');
end

% close all;
% figure; 
% plot(group_to_open.uv(1,:),group_to_open.uv(2,:),'b');
% hold on; 
% plot(opened_group(1,:),opened_group(2,:),'g');
% plot(rectangle_cut(1,:),rectangle_cut(2,:),'r');
% hold off;

end