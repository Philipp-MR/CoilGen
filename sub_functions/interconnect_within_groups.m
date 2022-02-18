function coil_parts=interconnect_within_groups(coil_parts,force_cut_selection)


coil_parts(numel(coil_parts)).connected_group=[];
coil_parts(numel(coil_parts)).return_paths=[];


for part_ind=1:numel(coil_parts)

%this options might not be needed anymore:
group_orientation=zeros(1,numel(coil_parts(part_ind).groups));



%choose between the high and low cut shapes
%choose the one which has significantly less points to remove
%avoid that the cuts cut muliple patches from the track

cut_selection=zeros(1,numel(coil_parts(part_ind).groups)); % 0 for low cut and 1 for the high cut
for group_num=1:numel(coil_parts(part_ind).groups)
    
    overlapp_points_high=0;
    overlapp_points_low=0;
    muliple_patches_high=[];
    muliple_patches_low=[];
    
    for loop_ind=1:numel(coil_parts(part_ind).groups(group_num).loops)
    in_points_high = inpolygon(coil_parts(part_ind).groups(group_num).loops(loop_ind).uv(1,:),coil_parts(part_ind).groups(group_num).loops(loop_ind).uv(2,:),coil_parts(part_ind).rectangle_cuts.high(group_num).uv(1,:),coil_parts(part_ind).rectangle_cuts.high(group_num).uv(2,:));
    in_points_low = inpolygon(coil_parts(part_ind).groups(group_num).loops(loop_ind).uv(1,:),coil_parts(part_ind).groups(group_num).loops(loop_ind).uv(2,:),coil_parts(part_ind).rectangle_cuts.low(group_num).uv(1,:),coil_parts(part_ind).rectangle_cuts.low(group_num).uv(2,:));
    overlapp_points_high=overlapp_points_high+sum(in_points_high);
    overlapp_points_low=overlapp_points_low+sum(in_points_low);
    muliple_patches_high=muliple_patches_high+sum(diff(in_points_high)<0)>1;
    muliple_patches_low=muliple_patches_low+sum(diff(in_points_low)<0)>1;
    end
    muliple_patches_high=muliple_patches_high~=0;
    muliple_patches_low=muliple_patches_low~=0;
    
    if  xor(muliple_patches_high,muliple_patches_high) % there is only one degenerate cut
            if muliple_patches_high
            cut_selection(group_num)=0;
            else
            cut_selection(group_num)=1;
            end
    else % both are or are not degenerate
            if overlapp_points_high >overlapp_points_low
            cut_selection(group_num)=0;
            else
            cut_selection(group_num)=1;
            end
    end   

    
end

% force cut selection if set
for group_num=1:numel(coil_parts(part_ind).groups)
    if ~isempty(force_cut_selection)
    if ~isnumeric(force_cut_selection{1})
        if strcmp(force_cut_selection{1},'high')
        cut_selection(group_num)=1;
        else
        cut_selection(group_num)=0;
        end   
    else
        if force_cut_selection(group_num)==1
        cut_selection(group_num)=1;
        else
        cut_selection(group_num)=0;
        end
    end
    end
end
%cut_selection(group_num)=coil_parts(part_ind).groups(group_num).loops(1).current_orientation;

%create the connected group
clear connected_group return_paths
connected_group(numel(coil_parts(part_ind).groups))=struct();
return_paths(numel(coil_parts(part_ind).groups))=struct();

planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices');


%if indicated flipp the interconnection order within the group ("inner first vs outer first")
coil_parts(part_ind).groups=flip_group_interconnection_order(coil_parts(part_ind).groups,group_orientation);


for group_num=1:numel(coil_parts(part_ind).groups)
 
if numel(coil_parts(part_ind).groups(group_num).loops)>1
if cut_selection(group_num)==1
rect_points=coil_parts(part_ind).rectangle_cuts.high(group_num).uv;
else
rect_points=coil_parts(part_ind).rectangle_cuts.low(group_num).uv;
end
[opened_loops,return_path_points,use_other_cut_rect_flag]=open_loop_group(coil_parts(part_ind).groups(group_num),rect_points);
%try the other opening cut, if the first one does not work
if use_other_cut_rect_flag==1
cut_selection(group_num)=~cut_selection(group_num);
if cut_selection(group_num)==1
rect_points=coil_parts(part_ind).rectangle_cuts.high(group_num).uv;
else
rect_points=coil_parts(part_ind).rectangle_cuts.low(group_num).uv;
end
[opened_loops,return_path_points,~]=open_loop_group(coil_parts(part_ind).groups(group_num),rect_points);
end
%build the interconnected group by adding the opened loops
connected_group(group_num).uv=[];
connected_group(group_num).v=[];
connected_group(group_num).spiral_in.uv=[];
connected_group(group_num).spiral_out.uv=[];
for loop_ind=1:numel(coil_parts(part_ind).groups(group_num).loops)
connected_group(group_num).uv=[connected_group(group_num).uv opened_loops{loop_ind}];
connected_group(group_num).spiral_in.uv=[connected_group(group_num).spiral_in.uv opened_loops{loop_ind}];
connected_group(group_num).spiral_out.uv=[connected_group(group_num).spiral_out.uv opened_loops{numel(coil_parts(part_ind).groups(group_num).loops)+1-loop_ind}];
end
%add the return path
connected_group(group_num).uv=[connected_group(group_num).uv return_path_points];
[connected_group(group_num).v,~]=uv_to_xyz(connected_group(group_num).uv,planary_mesh,curved_mesh);
return_paths(group_num).uv=return_path_points;
[return_paths(group_num).v,~]=uv_to_xyz(return_path_points,planary_mesh,curved_mesh);

[connected_group(group_num).spiral_in.v,~]=uv_to_xyz(connected_group(group_num).spiral_in.uv,planary_mesh,curved_mesh);
[connected_group(group_num).spiral_out.v,~]=uv_to_xyz(connected_group(group_num).spiral_out.uv,planary_mesh,curved_mesh);

else % group consists out of only one loop
   
connected_group(group_num).uv=coil_parts(part_ind).groups(group_num).loops(1).uv(:,[1:end-1]);
[connected_group(group_num).v,~]=uv_to_xyz(connected_group(group_num).uv,planary_mesh,curved_mesh);
    
connected_group(group_num).spiral_in.uv=connected_group(group_num).uv;
connected_group(group_num).spiral_out.uv=connected_group(group_num).uv;
connected_group(group_num).spiral_in.v=connected_group(group_num).v;
connected_group(group_num).spiral_in.v=connected_group(group_num).v;


return_paths(group_num).uv=mean(coil_parts(part_ind).groups(group_num).loops(1).uv(:,[1 end]),2);
[return_paths(group_num).v,~ ]=uv_to_xyz(return_paths(group_num).uv,planary_mesh,curved_mesh);
    
end

end
    
%Assign the outputs
coil_parts(part_ind).connected_group=connected_group;
coil_parts(part_ind).return_paths=return_paths;

end

function [opened_group,return_path,use_the_other_cut_rect_flag]=open_loop_group(group_to_open,rectangle_cut)
%open the loops with the cut rectangle and interconnect them

opened_group=cell(1,numel(group_to_open.loops));
return_path=[];
use_the_other_cut_rect_flag=0;
for loop_num=numel(group_to_open.loops):-1:1
%open the loops
raw_loop=group_to_open.loops(loop_num).uv(:,1:end-1);


[in_cut_rect_ind,~] = inpolygon(raw_loop(1,:),raw_loop(2,:),rectangle_cut(1,:),rectangle_cut(2,:)); %find points within cut rect
[buff_cut,~] = polyxpoly(raw_loop(1,:),raw_loop(2,:),rectangle_cut(1,:),rectangle_cut(2,:)); %find cutpoints of loop and cut rect

if any(in_cut_rect_ind) %there are points within cut rectanlge
    cut_points_inds=zeros(1,size(raw_loop,2));
    cut_points_inds(in_cut_rect_ind)=1;
    cut_points_inds=cut_points_inds+circshift(cut_points_inds,1)+circshift(cut_points_inds,-1);
    cut_points_inds(cut_points_inds>1)=1;
    first_half_inds=zeros(1,numel(cut_points_inds));
    second_half_inds=zeros(1,numel(cut_points_inds));
    first_half_inds(1:floor(numel(cut_points_inds)/2))=cut_points_inds(1:floor(numel(cut_points_inds)/2));
    second_half_inds(floor(numel(cut_points_inds)/2):end)=cut_points_inds( floor(numel(cut_points_inds)/2):end);
    if cut_points_inds(1)==1 && cut_points_inds(end)==1
    cut_points=raw_loop(:,[min(find(second_half_inds)):numel(cut_points_inds) 1:max(find(first_half_inds))]);
    else
    cut_points=raw_loop(:,find(cut_points_inds));
    end
    raw_loop=circshift(raw_loop,min(find(cut_points_inds))*(-1),2);
    [in_cut_rect_ind,~] = inpolygon(raw_loop(1,:),raw_loop(2,:),rectangle_cut(1,:),rectangle_cut(2,:)); 
    raw_loop(:,in_cut_rect_ind)=[];
    
    
    %check if still shifting needed (ugly workaround)
     [buff_cut_shift,~] = polyxpoly(raw_loop(1,:),raw_loop(2,:),rectangle_cut(1,:),rectangle_cut(2,:));
    if ~isempty(buff_cut_shift)
    %shift the points of the loop so that the opening cut is between the start
    %and end of the loop
    loop_to_shift=true;
    shift_ind_try_ind=0;
    while loop_to_shift %this solution for shifting works, but is slow..
    [buff_cut_shift,~] = polyxpoly(raw_loop(1,:),raw_loop(2,:),rectangle_cut(1,:),rectangle_cut(2,:));
    shift_ind_try_ind=shift_ind_try_ind+1;
    if ~isempty(buff_cut_shift)
    raw_loop=circshift(raw_loop,1,2);
    else 
    break;
    end
    if shift_ind_try_ind>2*size(raw_loop,2)
    use_the_other_cut_rect_flag=1;
    return;
    end
    end
    end
    
    %add the "cut"points for a clean cut
    [cut_x,cut_y] = polyxpoly(cut_points(1,:),cut_points(2,:),rectangle_cut(1,:),rectangle_cut(2,:));
    opened_group{loop_num}=raw_loop;
    %find out the order in which the cut points must inlcuded in the opened loop:
    test_dists_1=vecnorm(raw_loop(:,1)-[cut_x(1);cut_y(1)]);
    test_dists_2=vecnorm(raw_loop(:,1)-[cut_x(2);cut_y(2)]);
    if    test_dists_1<test_dists_2
        opened_group{loop_num}=[[cut_x(1) cut_y(1)]' opened_group{loop_num} [cut_x(2) cut_y(2)]' ] ;
    else 
        opened_group{loop_num}=[[cut_x(2) cut_y(2)]' opened_group{loop_num} [cut_x(1) cut_y(1)]' ] ;
    end
    
%     return_point=calc_3D_to_2D_mean([cut_x';cut_y'],parameterized_mesh);
%     %check if this point might be outside the mesh (might be the case for strong curved boundaries..)
%     [return_point_triangle,~] = pointLocation(triangulation(parameterized_mesh.f,parameterized_mesh.uv),return_point(1),return_point(2));
%     if isnan(return_point_triangle) %repair the return point 
%         %find the boundary which is cut by the retun path
%          boundary_cut_point=[1;1];
%          for boundary_test_ind=1:numel(parameterized_mesh.loops)
%         boundary_test_points=parameterized_mesh.uv(parameterized_mesh.loops{boundary_test_ind},:)';
%         [xi,yi] = polyxpoly(boundary_test_points(1,:),boundary_test_points(2,:),[return_path(1,:) return_point(1)],[return_path(2,:) return_point(2)]);
%          boundary_cut_point=[boundary_cut_point [xi; yi]];
%          end
%          boundary_cut_point(:,1)=[];
%          boundary_cut_point=return_path(:,end)+(boundary_cut_point-return_path(:,end)).*0.9; %shift the point inside the mesh
%          return_point=boundary_cut_point;
%     end
%  	return_path=[return_path return_point];
    
    longitudinal_vector=[cut_x(2)-cut_x(1); cut_y(2)-cut_y(1)];
    othorgonal_vector=[longitudinal_vector(2,:); longitudinal_vector(1,:)*(-1)];
    othorgonal_vector=othorgonal_vector./vecnorm(othorgonal_vector).*vecnorm([cut_x(2)-cut_x(1); cut_y(2)-cut_y(1)])./2;
    othorgonal_vector=othorgonal_vector./vecnorm(othorgonal_vector);
    [center_cut_x,center_cut_y]=polyxpoly(cut_points(1,:),cut_points(2,:),[mean(cut_x)+othorgonal_vector(1) mean(cut_x)-othorgonal_vector(1)],[mean(cut_y)+othorgonal_vector(2) mean(cut_y)-othorgonal_vector(2)]);
    return_path=[return_path [center_cut_x(1); center_cut_y(1)]];

else  %either no points within rect or loop already aligned with the endings around the cut rect

if   isempty(buff_cut) %the ends of the loops are already in the cut rect
    %add the "cut"points for a clean cut
    [cut_x,cut_y]=polyxpoly(raw_loop(1,[1:end 1]),raw_loop(2,[1:end 1]),rectangle_cut(1,:),rectangle_cut(2,:));
    opened_group{loop_num}=raw_loop;
    %find out the order in which the cut points must inlcuded in the opened loop:
    test_dists_1=vecnorm(raw_loop(:,1)-[cut_x(1);cut_y(1)]);
    test_dists_2=vecnorm(raw_loop(:,1)-[cut_x(2);cut_y(2)]);
    if    test_dists_1<test_dists_2
        opened_group{loop_num}=[[cut_x(1) cut_y(1)]' opened_group{loop_num} [cut_x(2) cut_y(2)]' ] ;
    else 
        opened_group{loop_num}=[[cut_x(2) cut_y(2)]' opened_group{loop_num} [cut_x(1) cut_y(1)]' ] ;
    end
    %return_path=[return_path calc_3D_to_2D_mean([cut_x';cut_y'],parameterized_mesh)];
    return_path=[return_path [mean(cut_x); mean(cut_y)]];
else % there are no points in the opening rectanlge but the also not the endings insde the opening rectangle
    
    %shift the points of the loop so that the opening cut is between the start
    %and end of the loop
    loop_to_shift=true;
    shift_ind_try_ind=0;
    while loop_to_shift %this solution for shifting works, but is slow..
    [buff_cut,~] = polyxpoly(raw_loop(1,:),raw_loop(2,:),rectangle_cut(1,:),rectangle_cut(2,:));
    shift_ind_try_ind=shift_ind_try_ind+1;
    if ~isempty(buff_cut)
    raw_loop=circshift(raw_loop,1,2);
    else 
    break;
    end
    if shift_ind_try_ind>2*size(raw_loop,2)
    use_the_other_cut_rect_flag=1;
    return;
    end
    end
    %add the "cut"points for a clean cut
    [cut_x,cut_y]=polyxpoly(raw_loop(1,[1:end 1]),raw_loop(2,[1:end 1]),rectangle_cut(1,:),rectangle_cut(2,:));
    opened_group{loop_num}=raw_loop;
    %find out the order in which the cut points must inlcuded in the opened loop:
    test_dists_1=vecnorm(raw_loop(:,1)-[cut_x(1);cut_y(1)]);
    test_dists_2=vecnorm(raw_loop(:,1)-[cut_x(2);cut_y(2)]);
    if    test_dists_1<test_dists_2
        opened_group{loop_num}=[[cut_x(1) cut_y(1)]' opened_group{loop_num} [cut_x(2) cut_y(2)]' ] ;
    else 
        opened_group{loop_num}=[[cut_x(2) cut_y(2)]' opened_group{loop_num} [cut_x(1) cut_y(1)]' ] ;
    end
    return_path=[return_path [mean(cut_x); mean(cut_y)]];
    %return_path=[return_path calc_3D_to_2D_mean([cut_x';cut_y'],parameterized_mesh)];
    
end


end



end

end


function group_container=flip_group_interconnection_order(group_container,group_orientation)
for group_ind=1:numel(group_container)
%flip the order of loops within the groups if indicated by group_ordering
if group_orientation(group_ind)==-1
v_cell=fliplr({group_container(group_ind).loops.point_coordinates});
uv_cell=fliplr({group_container(group_ind).loops.uv});
point_num_cell=fliplr({group_container(group_ind).loops.number_points});
potential_cell=fliplr({group_container(group_ind).loops.potential});
for vec_ind=1:numel(group_container(group_ind).loops)
group_container(group_ind).loops(vec_ind).point_coordinates=v_cell{vec_ind};
group_container(group_ind).loops(vec_ind).uv=uv_cell{vec_ind};
group_container(group_ind).loops(vec_ind).number_points=point_num_cell{vec_ind};
group_container(group_ind).loops(vec_ind).potential=potential_cell{vec_ind};
end
end
end
end



end


% function mean_out=calc_3D_to_2D_mean(in_points,parameterized_mesh)
% %build a mean in 3D by uv coords
% mean_out=mean(in_points,2);
% return;
% [target_triangles,bary_centric_coords] = pointLocation(triangulation(parameterized_mesh.f,parameterized_mesh.uv),in_points(1,:)',in_points(2,:)');
% if any(isnan(target_triangles))
% mean_out=mean(in_points,2);
% return;
% else
% points_3d = barycentricToCartesian(triangulation(parameterized_mesh.f,parameterized_mesh.v),target_triangles,bary_centric_coords)';
% %build the mean of these points
% mean_3d=mean(points_3d,2);
% %find the uv coords of this 3D mean point
% nearest_vertex_3D = nearestNeighbor(triangulation(parameterized_mesh.f,parameterized_mesh.v),mean_3d(1),mean_3d(2),mean_3d(3));
% %find all triangles with that mesh
% possible_triangles = vertexAttachments(triangulation(parameterized_mesh.f,parameterized_mesh.v),nearest_vertex_3D);
% possible_triangles_normals= faceNormal(triangulation(parameterized_mesh.f,parameterized_mesh.v),possible_triangles{1}');
% mean_3d_repmat=repmat(mean_3d',[numel(possible_triangles{1}) 1]);
% vert0=parameterized_mesh.v(parameterized_mesh.f(possible_triangles{1},1),:);
% vert1=parameterized_mesh.v(parameterized_mesh.f(possible_triangles{1},2),:);
% vert2=parameterized_mesh.v(parameterized_mesh.f(possible_triangles{1},3),:);
% [intersect, ~, u, v, ~] = TriangleRayIntersection (mean_3d_repmat, possible_triangles_normals, vert0, vert1, vert2, 'linetype','line');
% if any(intersect)
% intersect_tri=find(intersect);
% bary_u=u(intersect_tri(1));
% bary_v=v(intersect_tri(1));
% bary_w=1-bary_u-bary_v;
% %consider the possiblity that bary_u and bary_v might be swapped
% mean_out_1=barycentricToCartesian(triangulation(parameterized_mesh.f,parameterized_mesh.uv),possible_triangles{1}(intersect_tri(1)),[bary_u bary_v bary_w])';
% mean_out_2=barycentricToCartesian(triangulation(parameterized_mesh.f,parameterized_mesh.uv),possible_triangles{1}(intersect_tri(1)),[bary_v bary_u bary_w])';
% [~,min_ind_out]=min(vecnorm([mean_out_1 mean_out_2]-mean(in_points)'));
% if min_ind_out==1
% mean_out=mean_out_1;
% else
% mean_out=mean_out_2; 
% end
% else
% mean_out=mean(in_points,2); 
% end
% end
% end

% % out_canidates_2D=[];
% % out_canidates_3D=[];
% % for possible_tri_ind=1:numel(possible_triangles{1})
% % P0=parameterized_mesh.v(parameterized_mesh.f(possible_triangles{1}(possible_tri_ind),1),:)';
% % P1=parameterized_mesh.v(parameterized_mesh.f(possible_triangles{1}(possible_tri_ind),2),:)';
% % P2=parameterized_mesh.v(parameterized_mesh.f(possible_triangles{1}(possible_tri_ind),3),:)';
% % N=cross(P2-P1,P1-P0); % Normal to the plane of the triangle
% % N=N./vecnorm(N);
% % Q1=mean_3d-N;
% % Q2=mean_3d+N;
% % x = Q1 + (Q2-Q1).*dot(P0-Q1,N)./dot(Q2-Q1,N); % The point of intersection
% % v_n=cross(P2-P1,P1-P0);
% % A=vecnorm(v_n);
% % n_vec=v_n./A;
% % bary_u=sum(cross(P2-P1,x-P1).*n_vec)./A;
% % bary_v=sum(cross(P0-P2,x-P2).*n_vec)./A;
% % bary_w=1-bary_u-bary_v;
% % if (1>bary_u && bary_u>0) && (1>bary_v && bary_v>0) && (1>bary_w && bary_w>0)
% % out_canidates_2D=[out_canidates_2D barycentricToCartesian(triangulation(parameterized_mesh.f,parameterized_mesh.uv),possible_triangles{1}(possible_tri_ind),[bary_u bary_v bary_w])'];
% % out_canidates_3D=[out_canidates_3D x];
% % end
% % end
% % if isempty(out_canidates_2D)
% % mean_out=mean(in_points,2);
% % return;
% % else %take the points nearest to the actual 3D mean point
% % [~,min_ind_out]=min(vecnorm(out_canidates_3D-mean_3d));
% % mean_out=out_canidates_2D(:,min_ind_out);
% % end




% %open the loops with the cut rectangle and interconnect them
% for loop_num=1:numel(group_container(group_num).loops)
% %open the loops
% raw_loop=group_container(group_num).loops(loop_num).uv;
% %find loop points within the cut rect and points on intersection of
% %cut_rect and loops
% [cut_x,cut_y] = polyxpoly(raw_loop(1,:),raw_loop(2,:),rect_points(1,:),rect_points(2,:));
% longitudinal_vector=[cut_x(2)-cut_x(1); cut_y(2)-cut_y(1)];
% othorgonal_vector=[longitudinal_vector(2,:); longitudinal_vector(1,:)*(-1)];
% othorgonal_vector=othorgonal_vector./vecnorm(othorgonal_vector).*vecnorm([cut_x(2)-cut_x(1); cut_y(2)-cut_y(1)])./2;
% cut_center=[mean(cut_x); mean(cut_y)];
% %find the position of the cut within the loop
% [~,~,cut_segment_inds] = polyxpoly(raw_loop(1,:),raw_loop(2,:),[cut_center(1)+othorgonal_vector(1) cut_center(1)-othorgonal_vector(1)],[cut_center(2)+othorgonal_vector(2) cut_center(2)-othorgonal_vector(2)]);
% %raw_loop(:,end)=[]; %open the close end
% raw_loop=circshift(raw_loop,(-1)*(cut_segment_inds(1)-1),2);
% 
% raw_loop(:,size(raw_loop,2)+(-1)*(cut_segment_inds(1)-1))=[]; %open the closed end
% 
% [in_cut_rect_ind,~] = inpolygon(raw_loop(1,:),raw_loop(2,:),rect_points(1,:),rect_points(2,:));
% if ~isempty(in_cut_rect_ind)
% raw_loop(:,in_cut_rect_ind)=[];
% end
% opened_loops{loop_num}=raw_loop;
% %add the "cut"points for a clean cut
% if dot([cut_x(1) cut_y(1)]'-[cut_center(1) cut_center(2)]',opened_loops{loop_num}(:,1)-[cut_center(1) cut_center(2)]')>0
% opened_loops{loop_num}=[[cut_x(1) cut_y(1)]' opened_loops{loop_num} [cut_x(2) cut_y(2)]' ] ;
% else
% opened_loops{loop_num}=[[cut_x(2) cut_y(2)]' opened_loops{loop_num} [cut_x(1) cut_y(1)]' ] ;
% end 
% end



% for loop_num=1:numel(group_container(group_num).loops)
% 
%         %open the loops
%         [in_cut_rect_ind,~] = inpolygon(group_container(group_num).loops(loop_num).uv(1,:),group_container(group_num).loops(loop_num).uv(2,:),rect_points(1,:),rect_points(2,:));
%         in_cut_rect_ind=find(in_cut_rect_ind==1);
%         [cut_x,cut_y] = polyxpoly(group_container(group_num).loops(loop_num).uv(1,:),group_container(group_num).loops(loop_num).uv(2,:),rect_points(1,:),rect_points(2,:));
%         % delete the "inside" points and replace them with the "cut" points
%         %delete the "double" point at the loop end
%         in_cut_rect_ind(in_cut_rect_ind==size( group_container(group_num).loops(loop_num).uv,2))=[];
%         group_container(group_num).loops(loop_num).uv(:,end)=[]; 
%         %consider the possibility that the current loop opening is inside the
%         %opening rect
%         if ~isempty(intersect([size(group_container(group_num).loops(loop_num).uv,2) 1],in_cut_rect_ind))
%         group_container(group_num).loops(loop_num).uv(:,in_cut_rect_ind)=[];
%         opened_loops{loop_num}=group_container(group_num).loops(loop_num).uv;
%         else
%         group_container(group_num).loops(loop_num).uv(:,in_cut_rect_ind)=[];
%         opened_loops{loop_num}=circshift(group_container(group_num).loops(loop_num).uv,(-1)*(min(in_cut_rect_ind)-1),2);
%         end
%         %add the "cut"points for a clean cut
%         if vecnorm([cut_x(1) cut_y(1)]'-opened_loops{loop_num}(:,1))<vecnorm([cut_x(2) cut_y(2)]'-opened_loops{loop_num}(:,1))
%         opened_loops{loop_num}=[[cut_x(1) cut_y(1)]' opened_loops{loop_num} [cut_x(2) cut_y(2)]' ] ;
%         else
%         opened_loops{loop_num}=[[cut_x(2) cut_y(2)]' opened_loops{loop_num} [cut_x(1) cut_y(1)]' ] ;
%         end 
% end