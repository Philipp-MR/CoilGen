function [intersec_point,cut_flag]=plane_line_intersect(plane_normal,plane_pos,point_0,point_1)
intersec_point=zeros(1,3);
line_vec = point_1-point_0;
diff_vec = point_0 - plane_pos;
D = dot(plane_normal,line_vec);
N = -dot(plane_normal,diff_vec);
cut_flag=0;
if abs(D) < 10^-7        % The segment is parallel to plane
if N == 0           % The segment lies in plane
cut_flag=2;
return
else
cut_flag=0;       %no intersection
return
end
end
%compute the intersection parameter
sI = N / D;
intersec_point = point_0+ sI.*line_vec;
if (sI < -0.0000001 || sI > 1.0000001)
    cut_flag= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    cut_flag=1;
end
end