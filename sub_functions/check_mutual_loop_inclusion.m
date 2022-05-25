 function inside_flag=check_mutual_loop_inclusion(test_poly,target_poly)
%check if the testpolyon lies fully enclosed within the second polygon
% this check is done with the wnding number algorithm test for each vertex toward the
% second polygon
% Copyright: Philipp Amrein, Uniklinik Freiburg Mai 2022

winding_numbers=zeros(1,size(test_poly,2));

for point_ind=1:size(test_poly,2)

A=repmat(test_poly(:,point_ind),[1 size(target_poly,2)-1]);
B=target_poly(:,2:end);
C=target_poly(:,1:end-1);

vec1=C-A;
vec2=B-A;

angle = atan2( vec1(1,:).*vec2(2,:) - vec1(2,:).*vec2(1,:), vec1(1,:).*vec2(1,:) + vec1(2,:).*vec2(2,:) );

winding_numbers(point_ind)=round(abs(sum(angle)/(2*pi)));


end

if all(winding_numbers==1)
inside_flag=true;
else
inside_flag=false;
end



 end