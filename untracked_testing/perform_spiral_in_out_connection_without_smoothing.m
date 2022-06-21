function layout_out = perform_spiral_in_out_connection_without_smoothing(coil_track,track_width,layer_speration_width,parameterized_mesh,group_center,pcb_end_extention_factor,pcb_end_center_factor)
%Create for a layout a 2D sweep with defined with
%Mirror the output to create a second layer to have a spiral in spiral out
%design
%this works only for a single loop group!
%the mirroring works only if the cut path is along one of the axis


layer_1_out(numel(coil_track)).v=[];
layer_2_out(numel(coil_track)).v=[];

left_track(numel(coil_track)).spiral_out=[];
left_track(numel(coil_track)).spiral_in=[];
right_track(numel(coil_track)).spiral_out=[];
right_track(numel(coil_track)).spiral_in=[];

layer_1_tri(numel(coil_track)).uv=[];
layer_1_tri(numel(coil_track)).v=[];
layer_2_tri(numel(coil_track)).uv=[];
layer_2_tri(numel(coil_track)).v=[];

layer_1(numel(coil_track)).uv=[];
layer_1(numel(coil_track)).v=[];
layer_2(numel(coil_track)).uv=[];
layer_2(numel(coil_track)).v=[];

for group_ind=1:numel(coil_track)
% % %%Add a common postions between the later two layers
%swap_point_1=group_center.uv(:,group_ind);
swap_point_1=(coil_track(group_ind).spiral_in.uv(:,end)+coil_track(group_ind).spiral_out.uv(:,1))./2;
swap_point_2=(coil_track(group_ind).spiral_in.uv(:,1)+coil_track(group_ind).spiral_out.uv(:,end))./2;
swap_point_1=swap_point_1+(swap_point_1-group_center.uv(:,group_ind))./vecnorm((swap_point_1-group_center.uv(:,group_ind))).*pcb_end_center_factor;
swap_point_2=swap_point_2-(swap_point_2-group_center.uv(:,group_ind))./vecnorm((swap_point_2-group_center.uv(:,group_ind))).*pcb_end_center_factor;


coil_track(group_ind).spiral_out.uv=[swap_point_1 coil_track(group_ind).spiral_out.uv swap_point_2];

%Build the 2D sweeped track with right and left side border
[left_track(group_ind).spiral_out,right_track(group_ind).spiral_out]=build_2D_track(coil_track(group_ind).spiral_out.uv,track_width,parameterized_mesh);
% left_track(group_ind).spiral_out=fliplr(left_track(group_ind).spiral_out);
% right_track(group_ind).spiral_out=fliplr(right_track(group_ind).spiral_out);
% coil_track(group_ind).spiral_out.uv=fliplr(coil_track(group_ind).spiral_out.uv);



coil_track(group_ind).spiral_in.uv=[swap_point_1 fliplr(coil_track(group_ind).spiral_in.uv) swap_point_2 ];
[left_track(group_ind).spiral_in,right_track(group_ind).spiral_in]=build_2D_track(coil_track(group_ind).spiral_in.uv,track_width,parameterized_mesh);
left_track(group_ind).spiral_in=fliplr(left_track(group_ind).spiral_in);
right_track(group_ind).spiral_in=fliplr(right_track(group_ind).spiral_in);
coil_track(group_ind).spiral_in.uv=fliplr(coil_track(group_ind).spiral_in.uv);


coil_track(group_ind).spiral_in.v=uv_to_xyz_loc(coil_track(group_ind).spiral_in.uv,triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v));
coil_track(group_ind).spiral_out.v=uv_to_xyz_loc(coil_track(group_ind).spiral_out.uv,triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v));




%Stretch the ends a little to simplify connections between the layers
vec_add1=(left_track(group_ind).spiral_in(:,2)-left_track(group_ind).spiral_in(:,1))./vecnorm((left_track(group_ind).spiral_in(:,2)-left_track(group_ind).spiral_in(:,1))).*track_width.*pcb_end_extention_factor;
vec_add2=(left_track(group_ind).spiral_in(:,end)-left_track(group_ind).spiral_in(:,end-1))./vecnorm((left_track(group_ind).spiral_in(:,end)-left_track(group_ind).spiral_in(:,end-1))).*track_width.*pcb_end_extention_factor;
left_track(group_ind).spiral_in=[left_track(group_ind).spiral_in(:,1)-vec_add1 left_track(group_ind).spiral_in left_track(group_ind).spiral_in(:,end)+vec_add2];
right_track(group_ind).spiral_in=[right_track(group_ind).spiral_in(:,1)-vec_add1 right_track(group_ind).spiral_in right_track(group_ind).spiral_in(:,end)+vec_add2];
left_track(group_ind).spiral_out=[left_track(group_ind).spiral_out(:,1)-vec_add1 left_track(group_ind).spiral_out left_track(group_ind).spiral_out(:,end)+vec_add2];
right_track(group_ind).spiral_out=[right_track(group_ind).spiral_out(:,1)-vec_add1 right_track(group_ind).spiral_out right_track(group_ind).spiral_out(:,end)+vec_add2];

% left_track(group_ind).spiral_in=[left_track(group_ind).spiral_in(:,1)-(left_track(group_ind).spiral_in(:,2)-left_track(group_ind).spiral_in(:,1)) left_track(group_ind).spiral_in left_track(group_ind).spiral_in(:,end) left_track(group_ind).spiral_in(:,end)+(left_track(group_ind).spiral_in(:,end)-left_track(group_ind).spiral_in(:,end-1))];
% left_track(group_ind).spiral_out=[left_track(group_ind).spiral_out(:,1)-(left_track(group_ind).spiral_out(:,2)-left_track(group_ind).spiral_out(:,1)) left_track(group_ind).spiral_out left_track(group_ind).spiral_out(:,end) left_track(group_ind).spiral_out(:,end)+(left_track(group_ind).spiral_out(:,end)-left_track(group_ind).spiral_out(:,end-1))];
% right_track(group_ind).spiral_in=[right_track(group_ind).spiral_in(:,1)-(right_track(group_ind).spiral_in(:,2)-right_track(group_ind).spiral_in(:,1)) right_track(group_ind).spiral_in right_track(group_ind).spiral_in(:,end) right_track(group_ind).spiral_in(:,end)+(right_track(group_ind).spiral_in(:,end)-right_track(group_ind).spiral_in(:,end-1))];
% right_track(group_ind).spiral_out=[right_track(group_ind).spiral_out(:,1)-(right_track(group_ind).spiral_out(:,2)-right_track(group_ind).spiral_out(:,1)) right_track(group_ind).spiral_out right_track(group_ind).spiral_out(:,end) right_track(group_ind).spiral_out(:,end)+(right_track(group_ind).spiral_out(:,end)-right_track(group_ind).spiral_out(:,end-1))];


%build the conductor shape
layer_1(group_ind).uv=[left_track(group_ind).spiral_out fliplr(right_track(group_ind).spiral_out)];
layer_2(group_ind).uv=[left_track(group_ind).spiral_in fliplr(right_track(group_ind).spiral_in)];

layer_1(group_ind).v=uv_to_xyz_vec(layer_1(group_ind).uv,triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v));
layer_2(group_ind).v=uv_to_xyz_vec(layer_2(group_ind).uv,triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v));



warning('off','all');
poly1=polyshape(layer_1(group_ind).uv(1,:),layer_1(group_ind).uv(2,:));
poly2=polyshape(layer_2(group_ind).uv(1,:),layer_2(group_ind).uv(2,:));
warning('on','all');

layer_1_tri(group_ind).uv = triangulation(poly1);
layer_2_tri(group_ind) .uv= triangulation(poly2);

layer_1_tri(group_ind).v = triangulation(layer_1_tri(group_ind).uv.ConnectivityList,uv_to_xyz_vec(layer_1_tri(group_ind).uv.Points',triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v))');
layer_2_tri(group_ind).v = triangulation(layer_2_tri(group_ind).uv.ConnectivityList,uv_to_xyz_vec(layer_2_tri(group_ind).uv.Points',triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v))');

%Seperate the two layers along the normal
normals_layer_1=find_normals(layer_1_tri(group_ind).uv.Points,parameterized_mesh);
normals_layer_2=find_normals(layer_2_tri(group_ind).uv.Points,parameterized_mesh);

layer_1_out(group_ind).v=layer_1_tri(group_ind).v.Points+(normals_layer_1)'.*layer_speration_width/2;
layer_1_out(group_ind).f=layer_1_tri(group_ind).v.ConnectivityList;

layer_2_out(group_ind).v=layer_2_tri(group_ind).v.Points-(normals_layer_2)'.*layer_speration_width/2;
layer_2_out(group_ind).f=layer_2_tri(group_ind).v.ConnectivityList;

end

%Assign output
layout_out.layer_1_seperated=layer_1_out;
layout_out.layer_2_seperated=layer_2_out;
layout_out.layer_1_tri=layer_1_tri;
layout_out.layer_2_tri=layer_2_tri;
layout_out.layer_1=layer_1;
layout_out.layer_2=layer_2;
layout_out.coil_track=coil_track;

% % % % % %Plot
% % % % % if plot_pcb_layout
% % % % % figure; 
% % % % % hold on;
% % % % % trisurf(triangulation(parameterized_mesh.f,parameterized_mesh.v),'EdgeColor',[0.8 0.8 0.8],'faceAlpha',0);
% % % % % for plot_ind=1:numel(layer_1_out)
% % % % % trisurf(triangulation(layer_1_out(plot_ind).f,layer_1_out(plot_ind).v),'edgealpha',0,'facecolor','r');
% % % % % trisurf(triangulation(layer_2_out(plot_ind).f,layer_2_out(plot_ind).v),'edgealpha',0,'facecolor','b');
% % % % % end
% % % % % axis equal;
% % % % % hold off;
% % % % % 
% % % % % figure; 
% % % % % hold on;
% % % % % triplot(triangulation(parameterized_mesh.f,parameterized_mesh.uv),'color',[0.8 0.8 0.8]);
% % % % % for plot_ind=1:numel(layer_1_out)
% % % % % triplot(layer_1_tri(plot_ind).uv,'color','r');
% % % % % triplot(layer_2_tri(plot_ind).uv,'color','b');
% % % % % end
% % % % % axis equal;
% % % % % scatter([swap_point_1(1) swap_point_2(1)],[swap_point_1(2) swap_point_2(2)],'filled');
% % % % % hold off;
% % % % % 
% % % % % figure; 
% % % % % hold on;
% % % % % for plot_ind=1:numel(layer_1_out)
% % % % % plot(left_track(plot_ind).spiral_out(1,:),left_track(plot_ind).spiral_out(2,:));
% % % % % plot(left_track(plot_ind).spiral_in(1,:),left_track(plot_ind).spiral_in(2,:));
% % % % % plot(right_track(plot_ind).spiral_out(1,:),right_track(plot_ind).spiral_out(2,:));
% % % % % plot(right_track(plot_ind).spiral_in(1,:),right_track(plot_ind).spiral_in(2,:));
% % % % % scatter([swap_point_1(1) swap_point_2(1)],[swap_point_1(2) swap_point_2(2)],'filled');
% % % % % end
% % % % % hold off;
% % % % % 
% % % % % figure;
% % % % % hold on;
% % % % % for plot_ind=1:numel(layer_1_out)
% % % % % plot(coil_track(plot_ind).spiral_in.uv(1,:),coil_track(plot_ind).spiral_in.uv(2,:));
% % % % % plot(coil_track(plot_ind).spiral_out.uv(1,:),coil_track(plot_ind).spiral_out.uv(2,:));
% % % % % end
% % % % % scatter([swap_point_1(1) swap_point_2(1)],[swap_point_1(2) swap_point_2(2)],'filled');
% % % % % hold off;
% % % % % 
% % % % % 
% % % % % end


%assign output

% layout_out.tri_layer1_2d=tri1_2d;
% layout_out.tri_layer2_2d=tri2_2d;
% layout_out.tri_layer1_3d=tri1_3d;
% layout_out.tri_layer2_3d=tri2_3d;
% 
% layout_out.track1=coil_track_2d;
% layout_out.track2=fliplr([coil_track_2d(1,:).*(-1); coil_track_2d(2,:)]);
% 
% layout_out.layer1_points=layer_1;
% layout_out.layer2_points=layer_2;


function sampled_points = equilize_point_distances_spline(unsampled_points,upsampling_factor)
%Equilize the point the distances, keeping the same number of points
%delete non-unique points
ind_point_x=diff(unsampled_points(1,:))==0;
ind_point_y=diff(unsampled_points(2,:))==0;
unique_points=~(ind_point_x|ind_point_y);
unsampled_points=unsampled_points(:,unique_points')';
number_points=size(unsampled_points,1);
t=0:(1/((number_points-1)*upsampling_factor)):1;
t = t(:);
nt = numel(t);
ndim = 2;
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

function [x0,y0,segments]=selfintersect(x,y)
% Two similar curves are firstly created.
x1=x; x2=x;
y1=y; y2=y;
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);
% Compute number of line segments in each curve and some differences we'll
% need later.
n1 = length(x1) - 1;
n2 = length(x2) - 1;
dxy1 = diff([x1 y1]);
dxy2 = diff([x2 y2]);
[i,j] = find(repmat(min(x1(1:end-1),x1(2:end)),1,n2) <= ...
repmat(max(x2(1:end-1),x2(2:end)).',n1,1) & ...
repmat(max(x1(1:end-1),x1(2:end)),1,n2) >= ...
repmat(min(x2(1:end-1),x2(2:end)).',n1,1) & ...
repmat(min(y1(1:end-1),y1(2:end)),1,n2) <= ...
repmat(max(y2(1:end-1),y2(2:end)).',n1,1) & ...
repmat(max(y1(1:end-1),y1(2:end)),1,n2) >= ...
repmat(min(y2(1:end-1),y2(2:end)).',n1,1));
% Removing coincident and adjacent segments.
remove=find(abs(i-j)<2);
i(remove)=[];
j(remove)=[];
% Removing duplicate combinations of segments.
remove=[];
for ii=1:size(i,1)
ind=find((i(ii)==j(ii:end))&(j(ii)==i(ii:end)));
remove=[remove;ii-1+ind];
end
i(remove)=[];
j(remove)=[];
remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
i(remove) = [];
j(remove) = [];
% need anyway) and addition.
remove = isnan(sum(dxy1(i,:) + dxy2(j,:),2));
i(remove) = [];
j(remove) = [];
n = length(i);
T = zeros(4,n);
A = sparse([1 2 3 4],[3 3 4 4],-1,4,4,8);
B = -[x1(i) x2(j) y1(i) y2(j)].';
index_dxy1 = [1 3];  %  A(1) = A(1,1), A(3) = A(3,1)
index_dxy2 = [6 8];  %  A(6) = A(2,2), A(8) = A(4,2)
% Loop through possibilities.  Set warning not to trigger for anomalous
% results (i.e., when A is singular).
warning_state = warning('off','MATLAB:singularMatrix');
try
for k = 1:n
A(index_dxy1) = dxy1(i(k),:);
A(index_dxy2) = dxy2(j(k),:);
T(:,k) = A\B(:,k);
end
warning(warning_state)
catch
warning(warning_state)
rethrow(lasterror)
end
in_range = T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1;
anomalous = any(isnan(T));
if any(anomalous)
ia = i(anomalous);
ja = j(anomalous);
% set x0 and y0 to middle of overlapping region.
T(3,anomalous) = (max(min(x1(ia),x1(ia+1)),min(x2(ja),x2(ja+1))) + ...
min(max(x1(ia),x1(ia+1)),max(x2(ja),x2(ja+1))))/2;
T(4,anomalous) = (max(min(y1(ia),y1(ia+1)),min(y2(ja),y2(ja+1))) + ...
min(max(y1(ia),y1(ia+1)),max(y2(ja),y2(ja+1))))/2;
x0 = T(3,in_range | anomalous).';
y0 = T(4,in_range | anomalous).';
i=i(in_range | anomalous);
j=j(in_range | anomalous);
else
x0 = T(3,in_range).';
y0 = T(4,in_range).';
i=i(in_range);
j=j(in_range);
end
segments=sort([i,j],2);
end


function path_angles=calculate_track_angles(track_in)
%angles along a point track
path_vecs=track_in(:,2:end)-track_in(:,1:end-1);
path_vecs=[path_vecs path_vecs(:,end)];
vecs_a=path_vecs(:,2:end);
vecs_b=path_vecs(:,1:end-1);
path_angles = atan2(vecnorm(cross(vecs_a,vecs_b)),dot(vecs_a,vecs_b));
path_angles=path_angles./pi*180;
path_angles=[0 path_angles];
end




function [left_track,right_track]=build_2D_track(coil_track_2d,track_width_xyz,parameterized_mesh)

%coil_track_2d = equilize_point_distances_spline(coil_track_2d,1);

coil_track_2d =interparc([0:1/(numel(coil_track_2d(1,:))-1):1],coil_track_2d(1,:),coil_track_2d(2,:),'spline')';

%prepare the track direction
path_directions=coil_track_2d(:,2:end)-coil_track_2d(:,1:end-1); %edge vectors
path_directions=[path_directions path_directions(:,end)]; %add a repetition at the end for the last face
path_directions=path_directions./repmat(vecnorm(path_directions),[2 1]); % normalize
path_directions(:,2:end-1)=(path_directions(:,1:end-2)+path_directions(:,3:end))./2; % average the vector between its precessor and suczsecor
path_directions=path_directions./repmat(vecnorm(path_directions),[2 1]); % normalize again
othorgonal_vecs=[path_directions(2,:); path_directions(1,:)*(-1)];
othorgonal_vecs=othorgonal_vecs./repmat(vecnorm(othorgonal_vecs),[2 1]); % normalize again


%find the track width in the uv representation along the track
local_track_width=calc_local_track_width(parameterized_mesh,coil_track_2d,othorgonal_vecs,track_width_xyz);


%build the left and right edge track
left_track=coil_track_2d+othorgonal_vecs.*repmat(local_track_width,[2 1]);
right_track=coil_track_2d-othorgonal_vecs.*repmat(local_track_width,[2 1]);

%remove overly sharp edges
right_angles=calculate_track_angles([right_track; zeros(1,size(right_track,2))]);
left_angles=calculate_track_angles([left_track; zeros(1,size(left_track,2))]);
right_track=right_track(:,right_angles<35);
left_track=left_track(:,left_angles<35);
%smoth both edge tracks
left_track =interparc([0:1/(numel(left_track(1,:))-1):1],left_track(1,:),left_track(2,:),'spline')';
right_track =interparc([0:1/(numel(right_track(1,:))-1):1],right_track(1,:),right_track(2,:),'spline')';

% left_track(group_ind) = equilize_point_distances_spline(left_track(group_ind),1);
% right_track(group_ind) = equilize_point_distances_spline(right_track(group_ind),1);

end


function local_track_width=calc_local_track_width(parameterized_mesh,track_points,othorgonal_vector,opening_gab)
    
local_track_width=zeros(1,size(track_points,2));

opening_gab=opening_gab/(1000000);

for point_ind=1:size(track_points,2)

point_aa=track_points(:,point_ind)+opening_gab*othorgonal_vector(:,point_ind)/2;
point_bb=track_points(:,point_ind)-opening_gab*othorgonal_vector(:,point_ind)/2;
point_a_curved=uv_to_xyz(point_aa,triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v));
point_b_curved=uv_to_xyz(point_bb,triangulation(parameterized_mesh.f,parameterized_mesh.uv),triangulation(parameterized_mesh.f,parameterized_mesh.v));
curved_length=vecnorm(point_b_curved-point_a_curved);
local_track_width(:,point_ind)=opening_gab*opening_gab/curved_length;
local_track_width(:,point_ind)=local_track_width(:,point_ind)*(1000000);

end
end

function xyz=uv_to_xyz(uv,planary_mesh,curved_mesh)
% from the 2D surface coordinates to the xyz coordinates of the 3D coil
% surface
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,uv(1),uv(2));
xyz = barycentricToCartesian(curved_mesh,target_triangle,bary_centric_coord)';
end

function xyz=uv_to_xyz_vec(uv,planary_mesh,curved_mesh)
% from the 2D surface coordinates to the xyz coordinates of the 3D coil
% surface
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,uv(1,:)',uv(2,:)');
xyz = barycentricToCartesian(curved_mesh,target_triangle,bary_centric_coord)';
end


function normals_out=find_normals(path_uv,parameterized_mesh)
%calculate the normal vectors along the wire path
normals_out=zeros(3,size(path_uv,1));
planary_mesh=triangulation(parameterized_mesh.f,parameterized_mesh.uv);
for point_ind=1:size(path_uv,1)
[target_triangle_normal,~] = pointLocation(planary_mesh,path_uv(point_ind,1),path_uv(point_ind,2));
normals_out(:,point_ind)=parameterized_mesh.fn(target_triangle_normal,:)';
end
end


function [pt,dudt,fofthandle] = interparc(t,px,py,varargin)
% interparc: interpolate points along a curve in 2 or more dimensions
% usage: pt = interparc(t,px,py)    % a 2-d curve
% usage: pt = interparc(t,px,py,pz) % a 3-d curve
% usage: pt = interparc(t,px,py,pz,pw,...) % a 4-d or higher dimensional curve
% usage: pt = interparc(t,px,py,method) % a 2-d curve, method is specified
% usage: [pt,dudt,fofthandle] = interparc(t,px,py,...) % also returns derivatives, and a function handle
%
% Interpolates new points at any fractional point along
% the curve defined by a list of points in 2 or more
% dimensions. The curve may be defined by any sequence
% of non-replicated points.
%
% arguments: (input)
%  t   - vector of numbers, 0 <= t <= 1, that define
%        the fractional distance along the curve to
%        interpolate the curve at. t = 0 will generate
%        the very first point in the point list, and
%        t = 1 yields the last point in that list.
%        Similarly, t = 0.5 will yield the mid-point
%        on the curve in terms of arc length as the
%        curve is interpolated by a parametric spline.
%
%        If t is a scalar integer, at least 2, then
%        it specifies the number of equally spaced
%        points in arclength to be generated along
%        the curve.
%
%  px, py, pz, ... - vectors of length n, defining
%        points along the curve. n must be at least 2.
%        Exact Replicate points should not be present
%        in the curve, although there is no constraint
%        that the curve has replicate independent
%        variables.
%
%  method - (OPTIONAL) string flag - denotes the method
%        used to compute the points along the curve.
%
%        method may be any of 'linear', 'spline', or 'pchip',
%        or any simple contraction thereof, such as 'lin',
%        'sp', or even 'p'.
%        
%        method == 'linear' --> Uses a linear chordal
%               approximation to interpolate the curve.
%               This method is the most efficient.
%
%        method == 'pchip' --> Uses a parametric pchip
%               approximation for the interpolation
%               in arc length.
%
%        method == 'spline' --> Uses a parametric spline
%               approximation for the interpolation in
%               arc length. Generally for a smooth curve,
%               this method may be most accurate.
%
%        method = 'csape' --> if available, this tool will
%               allow a periodic spline fit for closed curves.
%               ONLY use this method if your points should
%               represent a closed curve.
%               
%               If the last point is NOT the same as the
%               first point on the curve, then the curve
%               will be forced to be periodic by this option.
%               That is, the first point will be replicated
%               onto the end.
%
%               If csape is not present in your matlab release,
%               then an error will result.
%
%        DEFAULT: 'spline'
%
%
% arguments: (output)
%  pt - Interpolated points at the specified fractional
%        distance (in arc length) along the curve.
%
%  dudt - when a second return argument is required,
%       interparc will return the parametric derivatives
%       (dx/dt, dy/dt, dz/dt, ...) as an array.
%
%  fofthandle - a function handle, taking numbers in the interval [0,1]
%       and evaluating the function at those points.
%
%       Extrapolation will not be permitted by this call.
%       Any values of t that lie outside of the interval [0,1]
%       will be clipped to the endpoints of the curve.
%
% Example:
% % Interpolate a set of unequally spaced points around
% % the perimeter of a unit circle, generating equally
% % spaced points around the perimeter.
% theta = sort(rand(15,1))*2*pi;
% theta(end+1) = theta(1);
% px = cos(theta);
% py = sin(theta);
%
% % interpolate using parametric splines
% pt = interparc(100,px,py,'spline');
%
% % Plot the result
% plot(px,py,'r*',pt(:,1),pt(:,2),'b-o')
% axis([-1.1 1.1 -1.1 1.1])
% axis equal
% grid on
% xlabel X
% ylabel Y
% title 'Points in blue are uniform in arclength around the circle'
%
%
% Example:
% % For the previous set of points, generate exactly 6
% % points around the parametric splines, verifying
% % the uniformity of the arc length interpolant.
% pt = interparc(6,px,py,'spline');
%
% % Convert back to polar form. See that the radius
% % is indeed 1, quite accurately.
% [TH,R] = cart2pol(pt(:,1),pt(:,2))
% % TH =
% %       0.86005
% %        2.1141
% %       -2.9117
% %        -1.654
% %      -0.39649
% %       0.86005
% % R =
% %             1
% %        0.9997
% %        0.9998
% %       0.99999
% %        1.0001
% %             1
%
% % Unwrap the polar angles, and difference them.
% diff(unwrap(TH))
% % ans =
% %        1.2541
% %        1.2573
% %        1.2577
% %        1.2575
% %        1.2565
%
% % Six points around the circle should be separated by
% % 2*pi/5 radians, if they were perfectly uniform. The
% % slight differences are due to the imperfect accuracy
% % of the parametric splines.
% 2*pi/5
% % ans =
% %        1.2566
%
%
% See also: arclength, spline, pchip, interp1
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/15/2010
% unpack the arguments and check for errors
if nargin < 3
  error('ARCLENGTH:insufficientarguments', ...
    'at least t, px, and py must be supplied')
end
t = t(:);
if (numel(t) == 1) && (t > 1) && (rem(t,1) == 0)
  % t specifies the number of points to be generated
  % equally spaced in arclength
  t = linspace(0,1,t)';
elseif any(t < 0) || any(t > 1)
  error('ARCLENGTH:impropert', ...
    'All elements of t must be 0 <= t <= 1')
end
% how many points will be interpolated?
nt = numel(t);
% the number of points on the curve itself
px = px(:);
py = py(:);
n = numel(px);
% are px and py both vectors of the same length?
if ~isvector(px) || ~isvector(py) || (length(py) ~= n)
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of the same length')
elseif n < 2
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of length at least 2')
end
% compose px and py into a single array. this way,
% if more dimensions are provided, the extension
% is trivial.
pxy = [px,py];
ndim = 2;
% the default method is 'linear'
method = 'spline';
% are there any other arguments?
if nargin > 3
  % there are. check the last argument. Is it a string?
  if ischar(varargin{end})
    method = varargin{end};
    varargin(end) = [];
    
    % method may be any of {'linear', 'pchip', 'spline', 'csape'.}
    % any any simple contraction thereof.
    valid = {'linear', 'pchip', 'spline', 'csape'};
    [method,errstr] = validstring(method,valid);
    if ~isempty(errstr)
      error('INTERPARC:incorrectmethod',errstr)
    end
  end
  
  % anything that remains in varargin must add
  % an additional dimension on the curve/polygon
  for i = 1:numel(varargin)
    pz = varargin{i};
    pz = pz(:);
    if numel(pz) ~= n
      error('ARCLENGTH:improperpxorpy', ...
        'pz must be of the same size as px and py')
    end
    pxy = [pxy,pz]; %#ok
  end
  
  % the final number of dimensions provided
  ndim = size(pxy,2);
end
% if csape, then make sure the first point is replicated at the end.
% also test to see if csape is available
if method(1) == 'c'
  if exist('csape','file') == 0
    error('CSAPE was requested, but you lack the necessary toolbox.')
  end
  
  p1 = pxy(1,:);
  pend = pxy(end,:);
  
  % get a tolerance on whether the first point is replicated.
  if norm(p1 - pend) > 10*eps(norm(max(abs(pxy),[],1)))
    % the two end points were not identical, so wrap the curve
    pxy(end+1,:) = p1;
    nt = nt + 1;
  end
end
% preallocate the result, pt
pt = NaN(nt,ndim);
% Compute the chordal (linear) arclength
% of each segment. This will be needed for
% any of the methods.
chordlen = sqrt(sum(diff(pxy,[],1).^2,2));
% Normalize the arclengths to a unit total
chordlen = chordlen/sum(chordlen);
% cumulative arclength
cumarc = [0;cumsum(chordlen)];
% The linear interpolant is trivial. do it as a special case
if method(1) == 'l'
  % The linear method.
  
  % which interval did each point fall in, in
  % terms of t?
  [junk,tbins] = histc(t,cumarc); %#ok
  
  % catch any problems at the ends
  tbins((tbins <= 0) | (t <= 0)) = 1;
  tbins((tbins >= n) | (t >= 1)) = n - 1;
  
  % interpolate
  s = (t - cumarc(tbins))./chordlen(tbins);
  % be nice, and allow the code to work on older releases
  % that don't have bsxfun
  pt = pxy(tbins,:) + (pxy(tbins+1,:) - pxy(tbins,:)).*repmat(s,1,ndim);
  
  % do we need to compute derivatives here?
  if nargout > 1
    dudt = (pxy(tbins+1,:) - pxy(tbins,:))./repmat(chordlen(tbins),1,ndim);
  end
  
  % do we need to create the spline as a piecewise linear function?
  if nargout > 2
    spl = cell(1,ndim);
    for i = 1:ndim
      coefs = [diff(pxy(:,i))./diff(cumarc),pxy(1:(end-1),i)];
      spl{i} = mkpp(cumarc.',coefs);
    end
    
    %create a function handle for evaluation, passing in the splines
    fofthandle = @(t) foft(t,spl);
  end
  
  % we are done at this point
  return
end
% If we drop down to here, we have either a spline
% or csape or pchip interpolant to work with.
% compute parametric splines
spl = cell(1,ndim);
spld = spl;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for i = 1:ndim
  switch method
    case 'pchip'
      spl{i} = pchip(cumarc,pxy(:,i));
    case 'spline'
      spl{i} = spline(cumarc,pxy(:,i));
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
    case 'csape'
      % csape was specified, so the curve is presumed closed,
      % therefore periodic
      spl{i} = csape(cumarc,pxy(:,i),'periodic');
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
  end
  
  % and now differentiate them
  xp = spl{i};
  xp.coefs = xp.coefs*diffarray;
  xp.order = 3;
  spld{i} = xp;
end
% catch the case where there were exactly three points
% in the curve, and spline was used to generate the
% interpolant. In this case, spline creates a curve with
% only one piece, not two.
if (numel(cumarc) == 3) && (method(1) == 's')
  cumarc = spl{1}.breaks;
  n = numel(cumarc);
  chordlen = sum(chordlen);
end
% Generate the total arclength along the curve
% by integrating each segment and summing the
% results. The integration scheme does its job
% using an ode solver.
% polyarray here contains the derivative polynomials
% for each spline in a given segment
polyarray = zeros(ndim,3);
seglen = zeros(n-1,1);
% options for ode45
opts = odeset('reltol',1.e-9);
for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:ndim
    polyarray(j,:) = spld{j}.coefs(i,:);
  end
  
  % integrate the arclength for the i'th segment
  % using ode45 for the integral. I could have
  % done this part with quad too, but then it
  % would not have been perfectly (numerically)
  % consistent with the next operation in this tool.
  [tout,yout] = ode45(@(t,y) segkernel(t,y),[0,chordlen(i)],0,opts); %#ok
  seglen(i) = yout(end);
end
% and normalize the segments to have unit total length
totalsplinelength = sum(seglen);
cumseglen = [0;cumsum(seglen)];
% which interval did each point fall into, in
% terms of t, but relative to the cumulative
% arc lengths along the parametric spline?
[junk,tbins] = histc(t*totalsplinelength,cumseglen); %#ok
% catch any problems at the ends
tbins((tbins <= 0) | (t <= 0)) = 1;
tbins((tbins >= n) | (t >= 1)) = n - 1;
% Do the fractional integration within each segment
% for the interpolated points. t is the parameter
% used to define the splines. It is defined in terms
% of a linear chordal arclength. This works nicely when
% a linear piecewise interpolant was used. However,
% what is asked for is an arclength interpolation
% in terms of arclength of the spline itself. Call s
% the arclength traveled along the spline.
s = totalsplinelength*t;
% the ode45 options will now include an events property
% so we can catch zero crossings.
opts = odeset('reltol',1.e-9,'events',@ode_events);
ti = t;
for i = 1:nt
  % si is the piece of arc length that we will look
  % for in this spline segment.
  si = s(i) - cumseglen(tbins(i));
  
  % extract polynomials for the derivatives
  % in the interval the point lies in
  for j = 1:ndim
    polyarray(j,:) = spld{j}.coefs(tbins(i),:);
  end
  
  % we need to integrate in t, until the integral
  % crosses the specified value of si. Because we
  % have defined totalsplinelength, the lengths will
  % be normalized at this point to a unit length.
  %
  % Start the ode solver at -si, so we will just
  % look for an event where y crosses zero.
  [tout,yout,te,ye] = ode45(@(t,y) segkernel(t,y),[0,chordlen(tbins(i))],-si,opts); %#ok
  
  % we only need that point where a zero crossing occurred
  % if no crossing was found, then we can look at each end.
  if ~isempty(te)
    ti(i) = te(1) + cumarc(tbins(i));
  else
    % a crossing must have happened at the very
    % beginning or the end, and the ode solver
    % missed it, not trapping that event.
    if abs(yout(1)) < abs(yout(end))
      % the event must have been at the start.
      ti(i) = tout(1) + cumarc(tbins(i));
    else
      % the event must have been at the end.
      ti(i) = tout(end) + cumarc(tbins(i));
    end
  end
end
% Interpolate the parametric splines at ti to get
% our interpolated value.
for L = 1:ndim
  pt(:,L) = ppval(spl{L},ti);
end
% do we need to compute first derivatives here at each point?
if nargout > 1
  dudt = zeros(nt,ndim);
  for L = 1:ndim
    dudt(:,L) = ppval(spld{L},ti);
  end
end
% create a function handle for evaluation, passing in the splines
if nargout > 2
  fofthandle = @(t) foft(t,spl);
end
% ===============================================
%  nested function for the integration kernel
% ===============================================
  function val = segkernel(t,y) %#ok
    % sqrt((dx/dt)^2 + (dy/dt)^2 + ...)
    val = zeros(size(t));
    for k = 1:ndim
      val = val + polyval(polyarray(k,:),t).^2;
    end
    val = sqrt(val);
    
  end % function segkernel
% ===============================================
%  nested function for ode45 integration events
% ===============================================
  function [value,isterminal,direction] = ode_events(t,y) %#ok
    % ode event trap, looking for zero crossings of y.
    value = y;
    isterminal = ones(size(y));
    direction = ones(size(y));
  end % function ode_events

% ===============================================
%       end mainline - interparc
% ===============================================
%       begin subfunctions
% ===============================================
% ===============================================
%  subfunction for evaluation at any point externally
% ===============================================
function f_t = foft(t,spl)
% tool allowing the user to evaluate the interpolant at any given point for any values t in [0,1]
pdim = numel(spl);
f_t = zeros(numel(t),pdim);
% convert t to a column vector, clipping it to [0,1] as we do.
t = max(0,min(1,t(:)));
% just loop over the splines in the cell array of splines
for i = 1:pdim
  f_t(:,i) = ppval(spl{i},t);
end
end % function foft
function [str,errorclass] = validstring(arg,valid)
% validstring: compares a string against a set of valid options
% usage: [str,errorclass] = validstring(arg,valid)
%
% If a direct hit, or any unambiguous shortening is found, that
% string is returned. Capitalization is ignored.
%
% arguments: (input)
%  arg - character string, to be tested against a list
%        of valid choices. Capitalization is ignored.
%
%  valid - cellstring array of alternative choices
%
% Arguments: (output)
%  str - string - resulting choice resolved from the
%        list of valid arguments. If no unambiguous
%        choice can be resolved, then str will be empty.
%
%  errorclass - string - A string argument that explains
%        the error. It will be one of the following
%        possibilities:
%
%        ''  --> No error. An unambiguous match for arg
%                was found among the choices.
%
%        'No match found' --> No match was found among 
%                the choices provided in valid.
%
%        'Ambiguous argument' --> At least two ambiguous
%                matches were found among those provided
%                in valid.
%        
%
% Example:
%  valid = {'off' 'on' 'The sky is falling'}
%  
%
% See also: parse_pv_pairs, strmatch, strcmpi
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/25/2010
ind = find(strncmpi(lower(arg),valid,numel(arg)));
if isempty(ind)
  % No hit found
  errorclass = 'No match found';
  str = '';
elseif (length(ind) > 1)
  % Ambiguous arg, hitting more than one of the valid options
  errorclass = 'Ambiguous argument';
  str = '';
  return
else
  errorclass = '';
  str = valid{ind};
end
end % function validstring


end % mainline - interparc


function points_out_3d=uv_to_xyz_loc(points_in_2d,planary_mesh,curved_mesh)
% from the 2D surface coordinates to the xyz coordinates of the 3D coil
% surface
points_out_3d=zeros(3,size(points_in_2d,2));
for rrrr=1:size(points_in_2d,2)
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,points_in_2d(1,rrrr),points_in_2d(2,rrrr));
while isnan(target_triangle) % in case the point is directly between two triangles and no directly triangle indice can be found
[target_triangle,bary_centric_coord] = pointLocation(planary_mesh,points_in_2d(1,rrrr)+ eps('double')*10*(1-2*rand(1)),points_in_2d(2,rrrr)+ eps('double')*10*(1-2*rand(1)));
end
points_out_3d(:,rrrr) = barycentricToCartesian(curved_mesh,target_triangle,bary_centric_coord)';
end
end


end

