function coil_parts= shift_return_paths(coil_parts,input)
    
vertical_separation=input.normal_shift_length;


%local variables
vec_normal_local_smoothing_length=2;
shift_length=2;
smoothing_length=5;
up_sample_factor=1;
input.smooth_factor=1;

coil_parts(numel(coil_parts)).shift_array=[];
coil_parts(numel(coil_parts)).points_to_shift=[];

for part_ind=1:numel(coil_parts)

planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.vertices');
wire_path_out=coil_parts(part_ind).wire_path;


%equilze the point to point distances
% mean_dist=mean(vecnorm(wire_path_out.uv(:,1:end-1)-wire_path_out.uv(:,2:end)));
% wire_path_out.uv = equilze_point_distances(wire_path_out.uv,mean_dist);

% %delete superimposing points
% points_to_delete=[1 vecnorm(curved_wire_path_out.uv(:,2:end)-curved_wire_path_out.uv(:,1:end-1)) 1]==0;
% curved_wire_path_out.uv(:,points_to_delete)=[];
% normal_vectors_wire_path_out.uv(:,points_to_delete)=[];
% is_crosspoint(:,points_to_delete)=[];


% %re-equiize the point to point distance


if input.smooth_flag
wire_path_out.uv = smooth_track_by_folding(wire_path_out.uv ,input.smooth_factor);
%wire_path_out.uv = equilize_point_distances_spline(wire_path_out.uv,input.smooth_factor);
[wire_path_out.v,wire_path_out.uv]=uv_to_xyz(wire_path_out.uv,planary_mesh,curved_mesh);
end


% wire_path_out = upsample_loop(wire_path_out,1);
% %downsample the wire
% %if size(wire_path_out.uv,2)>max_point_number+1
% %[wire_buff,~,~] = interparc([0:1/max_point_number:1 1],wire_path_out.uv(1,:),wire_path_out.uv(2,:));
% if up_sample_factor~=1
% wire_path_out.uv = equilize_point_distances_spline(wire_path_out.uv,up_sample_factor);
% end

%close the final track
if all(wire_path_out.uv(:,end)==wire_path_out.uv(:,1))
wire_path_out.uv=[wire_path_out.uv wire_path_out.uv(:,1)];
wire_path_out.v=[wire_path_out.v wire_path_out.v(:,1)];
end


%convert the tracks to 3d
%[wire_path_out.v,wire_path_out.uv]=uv_to_xyz(wire_path_out.uv,planary_mesh,curved_mesh);




% % %convert the tracks to 3d
% % [wire_path_out.uv,wire_path_out.v,~]=uv_to_xyz_loc(wire_path_out.uv,zeros(1,size(wire_path_out.uv,2)),planary_mesh,curved_mesh);

%detect wire crossings
%[cross_segments,~]=self_intersect_iterative(wire_path_out.uv);

%detect wire crossings
[cross_points,cross_segments] = InterX(wire_path_out.uv);

% %Add the cross points to the wire track
% cross_pair_inds=zeros(size(cross_segments));
% for cross_ind=1:size(cross_points,2)
% first_ind=cross_segments(cross_ind,1);
% second_ind=cross_segments(cross_ind,2);
% if first_ind<second_ind
% wire_path_out.uv=[ wire_path_out.uv(:,1:first_ind) cross_points(:,cross_ind) wire_path_out.uv(:,first_ind+1:second_ind) cross_points(:,cross_ind) wire_path_out.uv(:,second_ind+1:end)];
% cross_pair_inds(cross_ind,:)= [first_ind+1; second_ind+2];
% %update the cross_segment since additional points have been added
% cross_segments(cross_segments>first_ind)=cross_segments(cross_segments>first_ind)+1;
% cross_segments(cross_segments>second_ind)=cross_segments(cross_segments>second_ind)+1;
% else
% wire_path_out.uv=[ wire_path_out.uv(:,1:second_ind) cross_points(:,cross_ind) wire_path_out.uv(:,second_ind+1:first_ind) cross_points(:,cross_ind) wire_path_out.uv(:,first_ind+1:end)];
% cross_pair_inds(cross_ind,:)=[second_ind+1; first_ind+2];
% %update the cross_segment since additional points have been added
% cross_segments(cross_segments>second_ind)=cross_segments(cross_segments>second_ind)+1;
% cross_segments(cross_segments>first_ind)=cross_segments(cross_segments>first_ind)+1;
% end
% end
    


%Group the segments with crossings

if ~isempty(cross_segments)

sorted_crossed_segments=sort([cross_segments(:,1)' cross_segments(:,2)']);
neighours_weight=zeros(size(sorted_crossed_segments));
scale_factor=sum(sorted_crossed_segments.^4);
%Mark segments with have crossed segments in its vicinity
for crossed_seg_ind=1:numel(sorted_crossed_segments)
check_dists=abs(sorted_crossed_segments(crossed_seg_ind)-sorted_crossed_segments);
check_dists(check_dists==0)=[];
neighours_weight(crossed_seg_ind)=min(check_dists);
end

segment_to_shift=false(1,size(wire_path_out.uv,2)-1);
%For the pairs of crossed segments select which one will be shifted ( by the segments neighbour information)
for crossed_pair_ind=1:size(cross_segments,1)
ind1=find(sorted_crossed_segments==cross_segments(crossed_pair_ind,1));
ind2=find(sorted_crossed_segments==cross_segments(crossed_pair_ind,2));
if numel(ind1)>1 | numel(ind2)>1 %take care of the possibility that a segment has multible crossings
if numel(ind1)>numel(ind2)
segment_to_shift(cross_segments(crossed_pair_ind,1))=true;
else
segment_to_shift(cross_segments(crossed_pair_ind,2))=true;
end
else
if neighours_weight(ind1)<neighours_weight(ind2)
segment_to_shift(cross_segments(crossed_pair_ind,1))=true;
else
segment_to_shift(cross_segments(crossed_pair_ind,2))=true;
end
end
end

%Find the point indices of those segments
points_to_shift=zeros(1,size(wire_path_out.uv,2));
points_to_shift(unique([find(segment_to_shift) find(segment_to_shift)+1]))=1;

%Calculate vertex normals for later
vertex_normal=vertexNormal(curved_mesh);
%make sure that they are pointing to the outisde of the surface
if mean(dot([vertexNormal(curved_mesh)]',[curved_mesh.Points-mean(curved_mesh.Points)]'))<0
vertex_normal=vertex_normal.*(-1);
end

%calculate the normal vectors along the wire path
normal_vectors_wire_path=zeros(3,size(wire_path_out.uv,2));
for point_ind=1:size(wire_path_out.uv,2)
[target_triangle_normal,~] = pointLocation(planary_mesh,wire_path_out.uv(1,point_ind),wire_path_out.uv(2,point_ind));
if ~isnan(target_triangle_normal)
nodes_target_triangle=curved_mesh(target_triangle_normal,:);
node_normals_target_triangle=vertex_normal(nodes_target_triangle,:);
normal_vectors_wire_path(:,point_ind)=mean(node_normals_target_triangle)';
else
normal_vectors_wire_path(:,point_ind)=normal_vectors_wire_path(:,point_ind-1);
end
end



%shift the cross points along the normal to avoid wire intersection
shift_array=conv(points_to_shift,[1:smoothing_length ones(1,shift_length).*smoothing_length (smoothing_length-1):-1:1]./smoothing_length,'same');
shift_array(shift_array>1)=1;
%make sure that shifting the points does not lead to a intersection with
%the surface (this might happen if the surface is has extensiv curved edges)
for point_ind=1:size(wire_path_out.uv,2)
if point_ind>vec_normal_local_smoothing_length && point_ind<(size(wire_path_out.uv,2)-vec_normal_local_smoothing_length)  
shift_vec=mean(normal_vectors_wire_path(:,(point_ind-floor(vec_normal_local_smoothing_length/2)):(point_ind+floor(vec_normal_local_smoothing_length/2))),2).*shift_array(point_ind).*vertical_separation;
else
shift_vec=normal_vectors_wire_path(:,point_ind).*shift_array(point_ind).*vertical_separation;
end
wire_path_out.v(:,point_ind)=wire_path_out.v(:,point_ind)+shift_vec;
end


else
    
shift_array=zeros(1,size(wire_path_out.uv,2));
points_to_shift=zeros(1,size(wire_path_out.uv,2));
    
    
end

coil_parts(part_ind).shift_array=shift_array;
coil_parts(part_ind).points_to_shift=points_to_shift;
coil_parts(part_ind).wire_path.uv=wire_path_out.uv;
coil_parts(part_ind).wire_path.v=wire_path_out.v;
  
end

%%%%%  subfunctions  %%%%%%
function points_in = equilze_point_distances(points_in,maximal_distance)
%limit the point to point distane in the 2d representation
zzzz=0;
for rrrr=2:size(points_in,2)  
point_dist=vecnorm(points_in(:,rrrr+zzzz)-points_in(:,rrrr+zzzz-1));
if point_dist>maximal_distance
%find the number and positions of the points to insert
num_ins_points=floor(point_dist/maximal_distance);
shift_vec_loc=(points_in(:,rrrr+zzzz)-points_in(:,rrrr+zzzz-1))./point_dist;
shift_vec_loc=shift_vec_loc.*maximal_distance;
points_to_insert=points_in(:,rrrr+zzzz-1)+repmat(shift_vec_loc,1,num_ins_points).*repmat(1:num_ins_points,2,1);
points_in=[points_in(:,1:(rrrr+zzzz-1)) points_to_insert points_in(:,(rrrr+zzzz):end)];
zzzz=zzzz+num_ins_points;   
end       
end
end




function [P,intersect_edge_inds] = InterX(L1,varargin)
    if nargin == 1
        L2 = L1;    hF = @lt;   %...Avoid the inclusion of common points
    else
        L2 = varargin{1}; hF = @le;
    end
       
    %...Preliminary stuff
    x1  = L1(1,:)';  x2 = L2(1,:);
    y1  = L1(2,:)';  y2 = L2(2,:);
    dx1 = diff(x1); dy1 = diff(y1);
    dx2 = diff(x2); dy2 = diff(y2);
    
    %...Determine 'signed distances'   
    S1 = dx1.*y1(1:end-1) - dy1.*x1(1:end-1);
    S2 = dx2.*y2(1:end-1) - dy2.*x2(1:end-1);
    
    C1 = feval(hF,D(bsxfun(@times,dx1,y2)-bsxfun(@times,dy1,x2),S1),0);
    C2 = feval(hF,D((bsxfun(@times,y1,dx2)-bsxfun(@times,x1,dy2))',S2'),0)';
    %...Obtain the segments where an intersection is expected
    [i,j] = find(C1 & C2); 
    if isempty(i)
        P = zeros(2,0);
        intersect_edge_inds=[];
        return; 
    end
    
    %...Transpose and prepare for output
    i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
    L = dy2(j).*dx1(i) - dy1(i).*dx2(j);
    i = i(L~=0); j=j(L~=0); L=L(L~=0);  %...Avoid divisions by 0
    
    intersect_edge_inds=unique(sort([i' j],2),'rows');
    %intersect_edge_inds=[i' j];
    
    %...Solve system of eqs to get the common points
    P = unique([dx2(j).*S1(i) - dx1(i).*S2(j), ...
                dy2(j).*S1(i) - dy1(i).*S2(j)]./[L L],'rows')';
            
    function u = D(x,y)
        u = bsxfun(@minus,x(:,1:end-1),y).*bsxfun(@minus,x(:,2:end),y);
    end
end


% % % function [crossing_segments,cross_points_to_add]=self_intersect_iterative(points)
% % % %Calculate     
% % % %to avoid memory problems splitt the curve into several parts
% % % %Philipp Amrein, 4.2021
% % % track_part_length=1000;
% % % 
% % % 
% % % 
% % % num_points=size(points,2);
% % % points_inds=1:num_points;
% % % 
% % % 
% % % if size(points,2)>track_part_length
% % % track_part_inds=[1:track_part_length:num_points num_points];
% % % if track_part_inds(end-1)==track_part_inds(end)
% % % track_part_inds(end)=[];
% % % end
% % % points_part(numel(track_part_inds)-1).inds=[];
% % % for part_ind=1:numel(track_part_inds)-1
% % % part_inds=track_part_inds(part_ind):(track_part_inds(part_ind+1));
% % % points_part(part_ind).coords=points(:,part_inds);
% % % points_part(part_ind).inds=part_inds;
% % % end
% % % else
% % % points_part.coords=points;
% % % points_part.inds=points_inds;
% % % end
% % % 
% % % crossing_segments=[];
% % % cross_points_to_add=[];
% % % for part_ind=1:numel(points_part)
% % % part_start_ind=points_part(part_ind).inds(1);
% % % part_end_ind=points_part(part_ind).inds(end);
% % % 
% % % %check for  crossings within this prt
% % % [within_x,within_y,part_intersect_segments]=InterX(points_part(part_ind).coords(1,:),points_part(part_ind).coords(2,:));
% % % cross_points_to_add=[cross_points_to_add; [within_x within_y]];
% % % part_intersect_segments=part_intersect_segments+(part_start_ind-1);
% % % crossing_segments=[crossing_segments; part_intersect_segments];
% % % 
% % % %check for crossings with all following parts  
% % % [rest_x,rest_y,rest_intersects] = polyxpoly(points_part(part_ind).coords(1,:),points_part(part_ind).coords(2,:),points(1,part_end_ind+1:num_points),points(2,part_end_ind+1:num_points));
% % % cross_points_to_add=[cross_points_to_add; [rest_x rest_y]];
% % % rest_intersects(:,1)=rest_intersects(:,1)+(part_start_ind-1);
% % % rest_intersects(:,2)=rest_intersects(:,2)+part_end_ind;
% % % crossing_segments=[crossing_segments; rest_intersects];
% % % 
% % % end
% % %     
% % % cross_points_to_add=cross_points_to_add';
% % % 
% % % 
% % % end



function points_out = upsample_loop(points_in,upsample_factor)
%limit the point to point distane in the 2d representation
points_out=points_in;
seg_lengths=vecnorm(points_out.v(:,2:end)-points_out.v(:,1:end-1));
target_segment_length=mean(seg_lengths)/upsample_factor;
%for each segment calculate the of points to insert
point_num_to_insert=floor(seg_lengths./target_segment_length)-1;
segment_cells_v=arrayfun(@(x) [points_in.v(:,x-1) points_in.v(:,x)],2:size(points_in.v,2),'UniformOutput',false);
segment_cells_uv=arrayfun(@(x) [points_in.uv(:,x-1) points_in.uv(:,x)],2:size(points_in.uv,2),'UniformOutput',false);
for seg_ind=1:numel(seg_lengths)
seg_vec_v=[segment_cells_v{seg_ind}(:,2)-segment_cells_v{seg_ind}(:,1)]./(point_num_to_insert(seg_ind)+1);
points_to_insert_v=segment_cells_v{seg_ind}(:,1)+seg_vec_v.*(1:point_num_to_insert(seg_ind));
segment_cells_v{seg_ind}=[segment_cells_v{seg_ind}(:,1) points_to_insert_v];
seg_vec_uv=[segment_cells_uv{seg_ind}(:,2)-segment_cells_uv{seg_ind}(:,1)]./(point_num_to_insert(seg_ind)+1);
points_to_insert_uv=segment_cells_uv{seg_ind}(:,1)+seg_vec_uv.*(1:point_num_to_insert(seg_ind));
segment_cells_uv{seg_ind}=[segment_cells_uv{seg_ind}(:,1) points_to_insert_uv];
end
points_out.v=[segment_cells_v{:} points_in.v(:,end)];
points_out.uv=[segment_cells_uv{:} points_in.uv(:,end)];
end


function sampled_points = equilize_point_distances_spline(unsampled_points,upsampling_factor)
%Equilize the point the distances, keeping the same number of points
number_points=size(unsampled_points,2);
t=0:(1/(number_points*upsampling_factor)):1;
t = t(:);
nt = numel(t);
ndim = 2;
%delete non-unique points
ind_point_x=diff(unsampled_points(1,:))==0;
ind_point_y=diff(unsampled_points(2,:))==0;
unique_points=~(ind_point_x|ind_point_y);
unsampled_points=unsampled_points(:,unique_points')';
pt = NaN(nt,ndim);
chordlen = sqrt(sum(diff(unsampled_points,[],1).^2,2));
chordlen = chordlen/sum(chordlen);
cumarc = [0;cumsum(chordlen)];
spline_interpol = cell(1,ndim);
spld = spline_interpol;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for spl_dim = 1:ndim
spline_interpol{spl_dim} = spline(cumarc,unsampled_points(:,spl_dim));
nc = numel(spline_interpol{spl_dim}.coefs);
if nc < 4
spline_interpol{spl_dim}.coefs = [zeros(1,4-nc),spline_interpol{spl_dim}.coefs];
spline_interpol{spl_dim}.order = 4;
end
xp = spline_interpol{spl_dim};
xp.coefs = xp.coefs*diffarray;
xp.order = 3;
spld{spl_dim} = xp;
end
cumarc = spline_interpol{1}.breaks;
n = numel(cumarc);
polyarray = zeros(ndim,3);
segment_length = zeros(n-1,1);
opts = odeset('reltol',1.e-9);
for iii = 1:spline_interpol{1}.pieces
for jjj = 1:ndim
polyarray(jjj,:) = spld{jjj}.coefs(iii,:);
end
[tout,yout] = ode45(@(t,y) segkernel(t,y,ndim,polyarray),[0,chordlen(iii)],0,opts); %#ok
segment_length(iii) = yout(end);
end
totalsplinelength = sum(segment_length);
cumulative_segment_length = [0;cumsum(segment_length)];
[junk,tbins] = histc(t*totalsplinelength,cumulative_segment_length); %#ok
tbins((tbins <= 0) | (t <= 0)) = 1;
tbins((tbins >= n) | (t >= 1)) = n - 1;
s = totalsplinelength*t;
opts = odeset('reltol',1.e-9,'events',@ode_events);
ti = t;
for ii = 1:nt
si = s(ii) - cumulative_segment_length(tbins(ii));
for jj = 1:ndim
polyarray(jj,:) = spld{jj}.coefs(tbins(ii),:);
end
[tout,yout,te,ye] = ode45(@(t,y) segkernel(t,y,ndim,polyarray),[0,chordlen(tbins(ii))],-si,opts); %#ok
if ~isempty(te)
ti(ii) = te(1) + cumarc(tbins(ii));
else
if abs(yout(1)) < abs(yout(end))
ti(ii) = tout(1) + cumarc(tbins(ii));
else
ti(ii) = tout(end) + cumarc(tbins(ii));
end
end
end
for L = 1:ndim
pt(:,L) = ppval(spline_interpol{L},ti);
end
dudt = zeros(nt,ndim);
for L = 1:ndim
dudt(:,L) = ppval(spld{L},ti);
end
function val = segkernel(t,y,dim_num,poly_array) %#ok
val = zeros(size(t));
for k = 1:dim_num
val = val + polyval(poly_array(k,:),t).^2;
end
val = sqrt(val);
end
function [value,isterminal,direction] = ode_events(t,y) %#ok
value = y;
isterminal = ones(size(y));
direction = ones(size(y));
end
sampled_points=transpose(pt);
end


end
