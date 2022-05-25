function [I,check]=calc_plane_line_intersection(n,V0,P0,P1)
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