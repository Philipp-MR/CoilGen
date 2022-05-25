function opened_loop=open_loop(loop_to_open,cut_shape)
%open the loops with the cut rectangle and interconnect them

raw_loop=loop_to_open.uv;

%close the raw loop
if raw_loop(1,1)~=raw_loop(1,end)&raw_loop(2,1)~=raw_loop(2,end)
raw_loop=[raw_loop raw_loop(:,end)]; %close the raw loop
end

%calculate the intersection points with the respective loop segments
intersection_points=find_segment_intersections(raw_loop,cut_shape); 
all_intersec_segments=arrayfun(@(x) intersection_points(x).segment_inds,1:numel(intersection_points),'UniformOutput',0);
all_intersec_positions=[intersection_points(:).uv];
num_intersec_points=sum(~isnan([all_intersec_segments{:}]));
[sorted_segments,sort_ind]=sort([all_intersec_segments{:}]);
sorted_positions=all_intersec_positions(:,sort_ind);
sorted_positions=sorted_positions(:,~isnan(sorted_segments));
sorted_segments=sorted_segments(~isnan(sorted_segments));
% cuts_per_side=arrayfun(@(x) all(~isnan(intersection_points(x).segment_inds)).*numel(intersection_points(x).segment_inds),1:numel(intersection_points));
% any_cut_flag=~isempty([intersection_points(:).uv]);
cut_center=mean(cut_shape,2);

%Build cuts between the loop and the cutshape
if num_intersec_points>0% Check wether there are any intersections


%check wether the loop start&end are outside or inside the cutshape
test_intersection_points=find_segment_intersections(cut_shape,[cut_center raw_loop(:,1)]);
if isnan(test_intersection_points(1).segment_inds)
    ends_within_cutshape=true;
else
    ends_within_cutshape=false;
end

%chose the cutted parts of the loop which has longest overlapp with the cutshape
if ~ends_within_cutshape %shift the cut pairs for 1 
for cut_pair_ind=1:num_intersec_points/2
seg_1=sorted_segments(1+2*(cut_pair_ind-1));
seg_2=sorted_segments(2+2*(cut_pair_ind-1));
cut_pair_lengths(cut_pair_ind)=sum(vecnorm(raw_loop(:,(seg_1+1):seg_2)-raw_loop(:,(seg_1):(seg_2-1))));
end
[~,min_ind]=max(cut_pair_lengths);
seg_1=sorted_segments(1+2*(min_ind-1));
seg_2=sorted_segments(2+2*(min_ind-1));
first_cut_position=sorted_positions(:,1+2*(min_ind-1));
last_cut_position=sorted_positions(:,2+2*(min_ind-1));
in_cut_rect_inds=(seg_1+1):seg_2;
else  %shift the order of the cutted segments for one position so that we have proper pairs
sorted_segments=circshift(sorted_segments,-1);
sorted_positions=circshift(sorted_positions,-1,2);
for cut_pair_ind=1:num_intersec_points/2
seg_1=sorted_segments(1+2*(cut_pair_ind-1));
seg_2=sorted_segments(2+2*(cut_pair_ind-1));
if cut_pair_ind~=num_intersec_points/2
in_cut_rect_inds=(seg_1+1):seg_2;
cut_pair_lengths(cut_pair_ind)=sum(vecnorm(raw_loop(:,(seg_1+1):seg_2)-raw_loop(:,(seg_1):(seg_2-1))));
else %the last pair will have a point indice overflow
in_cut_rect_inds=[(seg_1+1):size(raw_loop,2) 1:seg_2];
cut_pair_lengths(cut_pair_ind)=sum(vecnorm(raw_loop(:,[(seg_1+1):size(raw_loop,2) 1:seg_2])-raw_loop(:,[(seg_1):size(raw_loop,2) 1:seg_2-1])));
end
end
[~,min_ind]=max(cut_pair_lengths);
seg_1=sorted_segments(1+2*(min_ind-1));
seg_2=sorted_segments(2+2*(min_ind-1));
first_cut_position=sorted_positions(:,1+2*(min_ind-1));
last_cut_position=sorted_positions(:,2+2*(min_ind-1));
end

%delete the points within the cutshape
raw_loop(:,in_cut_rect_inds)=[];

%shift the opening around the cutshape
if ~ends_within_cutshape & ~isempty(in_cut_rect_inds)
raw_loop = circshift(raw_loop,(in_cut_rect_inds(1)-1)*(-1),2);
% else %point indice overflow
% raw_loop = circshift(raw_loop,(in_cut_rect_inds(end)-1)*(-1),2);
elseif ~ends_within_cutshape & isempty(in_cut_rect_inds) % take care of the possibilty of no points within the cutshape
raw_loop = circshift(raw_loop,(sorted_segments(1))*(-1),2);
end

%add end points for a clean proper cut
%check which cutpoint belongs to which open end of the loop
track1=[first_cut_position raw_loop last_cut_position];
track2=[last_cut_position raw_loop first_cut_position];
if sum(vecnorm(track1(:,2:end)-track1(:,1:end-1)))<sum(vecnorm(track2(:,2:end)-track2(:,1:end-1)))
opened_loop=track1;
else
opened_loop=track2;
end

else
error('Degenerate cut: Cutshape does not overlapp with the loop; change the interconnection width');
end

% close all;
% figure; 
% plot(loop_to_open.uv(1,:),loop_to_open.uv(2,:),'b');
% hold on; 
% plot(opened_loop(1,:),opened_loop(2,:),'g');
% plot(rectangle_cut(1,:),rectangle_cut(2,:),'r');
% hold off;

end