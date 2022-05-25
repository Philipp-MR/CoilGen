function intersection_points=find_segment_intersections(loop,test_polygon)
% Copyright@ Philipp Amrein, Uniklinik Freiburg , 2022

%find intersection points loop and a polyogon (2d)

intersection_points(size(test_polygon,2)-1).segment_inds=[];
intersection_points((size(test_polygon,2)-1)).uv=[];


%the test polygon should be closed if it is not a single edge

if size(test_polygon,2)>2
if test_polygon(1,1)~=test_polygon(1,end)&test_polygon(2,1)~=test_polygon(2,end) 
test_polygon=[test_polygon test_polygon(:,1)]; %open the raw loop
end
end


%iterate over the test segments
for seg_ind=1:(size(test_polygon,2)-1)

x_1=repmat(test_polygon(1,seg_ind),[1,size(loop,2)-1]);
x_2=repmat(test_polygon(1,seg_ind+1),[1,size(loop,2)-1]);
x_3=loop(1,1:end-1);
x_4=loop(1,2:end);

y_1=repmat(test_polygon(2,seg_ind),[1,size(loop,2)-1]);
y_2=repmat(test_polygon(2,seg_ind+1),[1,size(loop,2)-1]);
y_3=loop(2,1:end-1);
y_4=loop(2,2:end);

d1x = x_2 - x_1;    
d1y = y_2 - y_1;
d2x = x_4 - x_3;     
d2y = y_4 - y_3;

s = (-d1y.*(x_1-x_3)+d1x.*(y_1-y_3))./(-d2x.*d1y+d1x.*d2y);
t = ( d2x.*(y_1-y_3)-d2y.*(x_1-x_3))./(-d2x.*d1y+d1x.*d2y);

intersection_segment_inds=find((s >= 0) & (s <= 1) & (t >= 0) & (t <= 1));

if ~isempty(intersection_segment_inds)

x_out = x_1(intersection_segment_inds) + (t(intersection_segment_inds) .* d1x(intersection_segment_inds));
y_out = y_1(intersection_segment_inds) + (t(intersection_segment_inds) .* d1y(intersection_segment_inds));

intersection_points(seg_ind).segment_inds=[intersection_points(seg_ind).segment_inds intersection_segment_inds];
intersection_points(seg_ind).uv=[intersection_points(seg_ind).uv [x_out;y_out]];

else

intersection_points(seg_ind).segment_inds=nan;
intersection_points(seg_ind).uv=[intersection_points(seg_ind).uv [nan;nan]];
end

end

end