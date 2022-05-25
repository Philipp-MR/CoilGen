function coil_parts=calc_contours_by_triangular_potential_cuts(coil_parts)
%center the stream function potential around zero and add zeros around the
%periphery

part(numel(coil_parts)).raw=[]; %initialize container

for part_ind=1:numel(coil_parts)

coil_parts(part_ind).coil_mesh.vertices=coil_parts(part_ind).coil_mesh.vertices';
coil_parts(part_ind).coil_mesh.uv=coil_parts(part_ind).coil_mesh.uv';
coil_parts(part_ind).coil_mesh.faces=coil_parts(part_ind).coil_mesh.faces';

%build the edges of the mesh and for each edge the attached triangle
edge_nodes = edges(triangulation(coil_parts(part_ind).coil_mesh.faces,coil_parts(part_ind).coil_mesh.uv));
edge_attached_triangles = edgeAttachments(triangulation(coil_parts(part_ind).coil_mesh.faces,coil_parts(part_ind).coil_mesh.uv),edge_nodes);
num_attached_tris=cellfun(@numel,edge_attached_triangles);
edge_nodes=edge_nodes(num_attached_tris==2,:);
edge_attached_triangles_inds= [cellfun(@(x) x(1),edge_attached_triangles(num_attached_tris==2)) cellfun(@(x) x(2),edge_attached_triangles(num_attached_tris==2))];



%find the node indices of triangles for the corresponding edges
edge_attached_triangles=cell(size(edge_attached_triangles_inds));
for x_ind=1:size(edge_attached_triangles_inds,1)
for y_ind=1:size(edge_attached_triangles_inds,2)
edge_attached_triangles{x_ind,y_ind}=coil_parts(part_ind).coil_mesh.faces(edge_attached_triangles_inds(x_ind,y_ind),:);
end
end

%take only the edge oppoising nodes of these trianglesa
edge_opposed_nodes=zeros(size(edge_attached_triangles));
for x_ind=1:size(edge_attached_triangles_inds,1)
edge_opposed_nodes(x_ind,1)=setdiff(edge_attached_triangles{x_ind,1},edge_attached_triangles{x_ind,2});
edge_opposed_nodes(x_ind,2)=setdiff(edge_attached_triangles{x_ind,2},edge_attached_triangles{x_ind,1});
end


%test for all edges wether they cut one of the potential levels
rep_contour_level_list=repmat(coil_parts(part_ind).potential_level_list,[size(edge_nodes,1),1]);
edge_node_potentials=coil_parts(part_ind).stream_function(edge_nodes);
min_edge_potentials=min(edge_node_potentials,[],2) ;
max_edge_potentials=max(edge_node_potentials,[],2);
min_edge_potentials=repmat(min_edge_potentials,1,numel(coil_parts(part_ind).potential_level_list));
max_edge_potentials=repmat(max_edge_potentials,1,numel(coil_parts(part_ind).potential_level_list));
%
tri_below_pot_step=max_edge_potentials>rep_contour_level_list;
tri_above_pot_step=min_edge_potentials<rep_contour_level_list;
potential_cut_criteria=tri_below_pot_step&tri_above_pot_step;
potential_cut_criteria=double(potential_cut_criteria);
potential_cut_criteria(potential_cut_criteria == 0) = NaN;

%calculate for each edge which cuts a potential the uv coordinates
edge_lengths=vecnorm(coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,2),:)-coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,1),:),2,2);
edge_lengths=repmat(edge_lengths,1,numel(coil_parts(part_ind).potential_level_list));
edge_potential_span=edge_node_potentials(:,2)-edge_node_potentials(:,1);
edge_potential_span=repmat(edge_potential_span,1,numel(coil_parts(part_ind).potential_level_list));

pot_dist_to_step=rep_contour_level_list-repmat(edge_node_potentials(:,1),1,numel(coil_parts(part_ind).potential_level_list));
cut_point_distance_to_edge_node_1=abs(pot_dist_to_step./edge_potential_span.*edge_lengths);

u_component_edge_vectors=coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,2),1)-coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,1),1);
u_component_edge_vectors=repmat(u_component_edge_vectors,1,numel(coil_parts(part_ind).potential_level_list));
v_component_edge_vectors=coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,2),2)-coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,1),2);
v_component_edge_vectors=repmat(v_component_edge_vectors,1,numel(coil_parts(part_ind).potential_level_list));

first_edge_node_u=repmat(coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,1),1),1,numel(coil_parts(part_ind).potential_level_list));
first_edge_node_v=repmat(coil_parts(part_ind).coil_mesh.uv(edge_nodes(:,1),2),1,numel(coil_parts(part_ind).potential_level_list));
u_cut_point=potential_cut_criteria.*(first_edge_node_u+u_component_edge_vectors.*cut_point_distance_to_edge_node_1./edge_lengths);
v_cut_point=potential_cut_criteria.*(first_edge_node_v+v_component_edge_vectors.*cut_point_distance_to_edge_node_1./edge_lengths);


%create cell by sorting the cut points to the corresponding potential
%levels
potential_sorted_cut_points=cell(1,numel(coil_parts(part_ind).potential_level_list));

for pot_ind=1:numel(coil_parts(part_ind).potential_level_list)
potential_sorted_cut_points{pot_ind}=[u_cut_point(:,pot_ind) v_cut_point(:,pot_ind)  [1:size(edge_nodes,1)]'];
potential_sorted_cut_points{pot_ind}(isnan(potential_sorted_cut_points{pot_ind}(:,1)),:)=[];
end

%creat a struct with the unsorted points
empty_potential_groups=cellfun(@isempty,potential_sorted_cut_points);
part(part_ind).raw.unsorted_points(sum(~empty_potential_groups))=struct();
running_ind=1;
for struct_ind=find(~empty_potential_groups)
part(part_ind).raw.unsorted_points(running_ind).potential=coil_parts(part_ind).potential_level_list(struct_ind);
part(part_ind).raw.unsorted_points(running_ind).edge_ind=potential_sorted_cut_points{struct_ind}(:,3);
part(part_ind).raw.unsorted_points(running_ind).uv=potential_sorted_cut_points{struct_ind}(:,[1 2]);
running_ind=running_ind+1;
end


%seperate loops within potential groups
%building loops by edge connecitvity information



%select the one of two possible triangles which has not been used yet
part(part_ind).raw.unarranged_loops.loop.edge_inds=[];
part(part_ind).raw.unarranged_loops.loop.uv=[];


%create loops
for potential_group=1:numel(part(part_ind).raw.unsorted_points)
    
    
all_current_edges=edge_nodes(part(part_ind).raw.unsorted_points(potential_group).edge_ind,:);
all_current_opposed_nodes=edge_opposed_nodes(part(part_ind).raw.unsorted_points(potential_group).edge_ind,:);
all_current_uv_coords=part(part_ind).raw.unsorted_points(potential_group).uv;
% all_opposed_triangles=edge_attached_triangles_array(part(part_ind).raw.unsorted_points(potential_group).edge_ind,:);
% all_current_triangles=edge_attached_triangles(part(part_ind).raw.unsorted_points(potential_group).edge_ind,:);   

set_new_start=1;    
num_build_loops=0;
edge_already_used=zeros(size(all_current_edges,1),1);    


    %beginn to connect
    while ~all(edge_already_used)
        
        if set_new_start==1 %set new start if loop is closed within one potential group or line of segments has ended
        %initialize the starting position
        num_build_loops=num_build_loops+1;
        starting_edge=min(find(edge_already_used==0));
        part(part_ind).raw.unarranged_loops(potential_group).loop(num_build_loops).uv=[all_current_uv_coords(starting_edge,:)];
        part(part_ind).raw.unarranged_loops(potential_group).loop(num_build_loops).edge_inds=all_current_edges(starting_edge,:);
        %mark the start postion as "included"
        edge_already_used(starting_edge)=1;
        current_edge=starting_edge;
        %find the next edge
        %find the edges which contains the opposed nodes of the current edge as
        %well as one of the nodes of the current edge
        current_edge_nodes=all_current_edges(current_edge,:);
        neighbouring_free_next_edges=find(any(all_current_edges==all_current_opposed_nodes(current_edge,1) | ...
        all_current_edges==all_current_opposed_nodes(current_edge,2),2) & any(all_current_edges==current_edge_nodes(1) | all_current_edges==current_edge_nodes(2),2));
    
        if isempty(neighbouring_free_next_edges)
        break;
        elseif numel(neighbouring_free_next_edges)==1
        next_edge=neighbouring_free_next_edges;
        else
        %select as a stariting direction one of the two possible nodes
        if ~edge_already_used(neighbouring_free_next_edges(1))
        next_edge=neighbouring_free_next_edges(1);
        else
        next_edge=neighbouring_free_next_edges(2);
        end
        set_new_start=0;
        end
        end
        
        while ~(next_edge==starting_edge)
            %include the next point
            edge_already_used(next_edge)=1;
            part(part_ind).raw.unarranged_loops(potential_group).loop(num_build_loops).uv=[part(part_ind).raw.unarranged_loops(potential_group).loop(num_build_loops).uv; all_current_uv_coords(next_edge,:)];
            part(part_ind).raw.unarranged_loops(potential_group).loop(num_build_loops).edge_inds=[part(part_ind).raw.unarranged_loops(potential_group).loop(num_build_loops).edge_inds; all_current_edges(next_edge,:)];
            current_edge=next_edge;
            %find the next edge
            %find the edges which contains the opposed nodes of the current edge as
            %well as one of the nodes of the current edge
            current_edge_nodes=all_current_edges(current_edge,:);
            possible_next_edges=find(any(all_current_edges==all_current_opposed_nodes(current_edge,1) | ...
            all_current_edges==all_current_opposed_nodes(current_edge,2),2) & any(all_current_edges==current_edge_nodes(1) | all_current_edges==current_edge_nodes(2),2));
            possible_next_edges=setdiff(possible_next_edges,find(edge_already_used==1));
            %check if the starting edge is under the possible next edges
            if isempty(possible_next_edges)
            break;
            else
            if  numel(possible_next_edges)==1
            next_edge=possible_next_edges;
            else
            %select as a stariting direction one of the two possible nodes
            if ~edge_already_used(possible_next_edges(1))
            next_edge=possible_next_edges(1);
            else
            next_edge=possible_next_edges(2);
            end
            end
            end
        end
        set_new_start=1;
    end
end


%evalute for each loop the current orientation
part(part_ind).raw.unarranged_loops(numel(part(part_ind).raw.unsorted_points)).loop(1).current_orientation=[];
for pot_ind=1:numel(part(part_ind).raw.unsorted_points)
    center_segment_potential=part(part_ind).raw.unsorted_points(pot_ind).potential;
for loop_ind=1:numel(part(part_ind).raw.unarranged_loops(pot_ind).loop)
if numel(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).edge_inds)>2
test_edge=floor(size(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).edge_inds,1)/2);
first_edge=part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).edge_inds(test_edge,:);
second_edge=part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).edge_inds(test_edge+1,:);
node_1=intersect(first_edge,second_edge);
node_2=setdiff(first_edge,node_1);
node_3=setdiff(second_edge,node_1);
node_1_uv=coil_parts(part_ind).coil_mesh.uv(node_1,:);
node_2_uv=coil_parts(part_ind).coil_mesh.uv(node_2,:);
node_3_uv=coil_parts(part_ind).coil_mesh.uv(node_3,:);
node_1_pot=coil_parts(part_ind).stream_function(node_1);
node_2_pot=coil_parts(part_ind).stream_function(node_2);
node_3_pot=coil_parts(part_ind).stream_function(node_3);
%calculate the 2D gradient of the triangle
center_segment_positon=(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(test_edge,:)+part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(test_edge+1,:))./2;
vec_center_node_1=node_1_uv-center_segment_positon;
vec_center_node_2=node_2_uv-center_segment_positon;
vec_center_node_3=node_3_uv-center_segment_positon;
pot_diff_center_node_1=node_1_pot-center_segment_potential;
pot_diff_center_node_2=node_2_pot-center_segment_potential;
pot_diff_center_node_3=node_3_pot-center_segment_potential;
pot_gradient_vec=vec_center_node_1*pot_diff_center_node_1+vec_center_node_2*pot_diff_center_node_2+vec_center_node_3*pot_diff_center_node_3;
%test the chirality of the segment on the potential gradient on that
%segment
segment_vec=part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(test_edge+1,:)-part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(test_edge,:);
cross_vec=cross([segment_vec 0],[pot_gradient_vec 0]);
part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).current_orientation=sign(cross_vec(3));

if part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).current_orientation==-1
    part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv=flipud(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv);
    part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).edge_inds=flipud(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).edge_inds);
end
else
error('Some loops are to small and contain only 2 points, therefore ill-defined');
end
end
end



%build the contour lines
coil_parts(part_ind).contour_lines.uv=[];
coil_parts(part_ind).contour_lines.potential=[];
coil_parts(part_ind).contour_lines.current_orientation=[];
build_ind=1;
for pot_ind=1:numel(part(part_ind).raw.unsorted_points)
for loop_ind=1:numel(part(part_ind).raw.unarranged_loops(pot_ind).loop)
%coil_parts(part_ind).contour_lines(build_ind).uv=[part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv; part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(1,:)]';
coil_parts(part_ind).contour_lines(build_ind).uv=part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv';
coil_parts(part_ind).contour_lines(build_ind).potential=part(part_ind).raw.unsorted_points(pot_ind).potential;
%find the current orientation (for comparisson with other loops)
uv_center=mean(coil_parts(part_ind).contour_lines(build_ind).uv,2);
uv_to_center_vecs=coil_parts(part_ind).contour_lines(build_ind).uv(:,1:end-1)-uv_center;
uv_to_center_vecs=[uv_to_center_vecs; zeros(1,size(uv_to_center_vecs,2))];
uv_vecs=coil_parts(part_ind).contour_lines(build_ind).uv(:,2:end)-coil_parts(part_ind).contour_lines(build_ind).uv(:,1:end-1);
uv_vecs=[uv_vecs; zeros(1,size(uv_vecs,2))];
rot_vecs=cross(uv_to_center_vecs,uv_vecs);
track_orientation=sign(sum(rot_vecs(3,:)));
coil_parts(part_ind).contour_lines(build_ind).current_orientation=track_orientation;
build_ind=build_ind+1;
end
end

coil_parts(part_ind).coil_mesh.vertices=coil_parts(part_ind).coil_mesh.vertices';
coil_parts(part_ind).coil_mesh.uv=coil_parts(part_ind).coil_mesh.uv';
coil_parts(part_ind).coil_mesh.faces=coil_parts(part_ind).coil_mesh.faces';

end


% %plot points
% figure; hold on; axis equal;
% triplot(triangulation(parameterized_mesh.f,parameterized_mesh.uv),'color',[0.9 0.9 0.9]);
% for pot_ind=1:numel(potential_level_list)
% if ~isempty(potential_sorted_cut_points{pot_ind})
% scatter(potential_sorted_cut_points{pot_ind}(:,1),potential_sorted_cut_points{pot_ind}(:,2));
% end
% end
% hold off

% %plot build loops
% my_color=colorcube(3*numel(part(part_ind).raw.unarranged_loops));
% my_color(find([[my_color(:,1)==my_color(:,2)].*[my_color(:,1)==my_color(:,3)]]'),:)=[]; %#ok
% color_ind=1;
% figure; hold on; axis equal;
% %triplot(triangulation(parameterized_mesh.f,parameterized_mesh.uv),'color',[0.9 0.9 0.9]);
% for pot_ind=1:numel(part(part_ind).raw.unarranged_loops)
% for loop_ind=1:numel(part(part_ind).raw.unarranged_loops(pot_ind).loop)
% if ~isempty(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv)
% plot(part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(:,1),part(part_ind).raw.unarranged_loops(pot_ind).loop(loop_ind).uv(:,2),'color',my_color(color_ind,:));
% end
% end
% color_ind=color_ind+1;
% end
% hold off


        
 
        
    

end