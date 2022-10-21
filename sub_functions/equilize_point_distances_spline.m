function sampled_points = equilize_point_distances_spline(unsampled_points,upsampling_factor)
%Equilize the point the distances, keeping the same number of points
unsampled_points=[unsampled_points(1,:) unsampled_points(1,end)+(unsampled_points(1,end)-unsampled_points(1,end-1));unsampled_points(2,:) unsampled_points(2,end)+(unsampled_points(2,end)-unsampled_points(2,end-1))];
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