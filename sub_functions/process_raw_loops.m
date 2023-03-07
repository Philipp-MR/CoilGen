function coil_parts= process_raw_loops(coil_parts,input,target_field)
%take care of loops crossing boundaries or off boundary loops

%Smooth the contours
if input.smooth_flag
for part_ind=1:numel(coil_parts)
for loop_num=1:numel(coil_parts(part_ind).contour_lines)
coil_parts(part_ind).contour_lines(loop_num).uv= smooth_track_by_folding(coil_parts(part_ind).contour_lines(loop_num).uv,input.smooth_factor);
end
end
end

%Generate the curved coordinates
for part_ind=1:numel(coil_parts)
planary_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.uv');
curved_mesh=triangulation(coil_parts(part_ind).coil_mesh.faces',coil_parts(part_ind).coil_mesh.v);
for loop_num=1:numel(coil_parts(part_ind).contour_lines)
[coil_parts(part_ind).contour_lines(loop_num).v,coil_parts(part_ind).contour_lines(loop_num).uv]=uv_to_xyz(coil_parts(part_ind).contour_lines(loop_num).uv,planary_mesh,curved_mesh);
end
end


%Remove loops that do no contribute enough to the target field
% if input.min_loop_signifcance~=0
coil_parts=evaluate_loop_significance(coil_parts,target_field);
for part_ind=1:numel(coil_parts)
loops_to_delete=coil_parts(part_ind).loop_signficance<input.min_loop_signifcance;
coil_parts(part_ind).contour_lines(loops_to_delete)=[];
end
% else
% coil_parts(part_ind).combined_loop_field=zeros(size(target_field.b));
% coil_parts(part_ind).loop_signficance=ones(1,numel(coil_parts(part_ind).contour_lines));
% end

%Close the loops
for part_ind=1:numel(coil_parts)
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
if coil_parts(part_ind).contour_lines(loop_ind).uv(1,end)~=coil_parts(part_ind).contour_lines(loop_ind).uv(1,1) & coil_parts(part_ind).contour_lines(loop_ind).uv(2,end)~=coil_parts(part_ind).contour_lines(loop_ind).uv(2,1) 
coil_parts(part_ind).contour_lines(loop_ind).uv=[coil_parts(part_ind).contour_lines(loop_ind).uv coil_parts(part_ind).contour_lines(loop_ind).uv(:,1)]; %close the loops
coil_parts(part_ind).contour_lines(loop_ind).v=[coil_parts(part_ind).contour_lines(loop_ind).v coil_parts(part_ind).contour_lines(loop_ind).v(:,1)]; %close the loops
end
end
end


%Calculate the combinded wire length for the unconnected loops
coil_parts(numel(coil_parts)).combined_loop_length=0;
for part_ind=1:numel(coil_parts)
coil_parts(part_ind).combined_loop_length=sum(arrayfun(@(x) sum(vecnorm(coil_parts(part_ind).contour_lines(x).v(:,2:end)-coil_parts(part_ind).contour_lines(x).v(:,1:end-1))),1:numel(coil_parts(part_ind).contour_lines)));
end

function coil_parts=evaluate_loop_significance(coil_parts,target_field)
%Calcaltue relative errors between the different input and
%output fields
%Calcuate the combined field of the unconnected contours

coil_parts(numel(coil_parts)).combined_loop_field=[];
coil_parts(numel(coil_parts)).loop_signficance=[];

coil_parts(numel(coil_parts)).field_by_loops=[];
for part_inds=1:numel(coil_parts)
coil_parts(part_inds).field_by_loops=zeros(3,size(target_field.b,2),numel(coil_parts(part_inds).contour_lines));
for loop_inds=1:numel(coil_parts(part_inds).contour_lines)
coil_parts(part_inds).field_by_loops(:,:,loop_inds)=biot_savart_calc_b(coil_parts(part_inds).contour_lines(loop_inds).v,target_field).*coil_parts(part_inds).contour_step;
end
combined_loop_field=squeeze(sum(coil_parts(part_inds).field_by_loops,3));
loop_signficance=zeros(1,numel(coil_parts(part_inds).contour_lines));
for loop_inds=1:numel(coil_parts(part_inds).contour_lines)
%loop_signficance(loop_inds)=max(abs(coil_parts(part_inds).field_by_loops(3,:,loop_inds)))./(mean(abs(combined_loop_field(3,:)))/numel(coil_parts(part_inds).contour_lines)).*100;
loop_signficance(loop_inds)=max(abs(coil_parts(part_inds).field_by_loops(3,:,loop_inds)))./(mean(abs(combined_loop_field(3,:)))).*100;
end
coil_parts(part_inds).combined_loop_field=combined_loop_field;
coil_parts(part_inds).loop_signficance=loop_signficance;
end


function b_field=biot_savart_calc_b(wire_elements,target_f)
%Calculate b field with biot savarts law  by wire elements given as a sequence of coordinate points
num_tp=size(target_f.b,2);
%to avoid memory problems splitt the curve into several parts
track_part_length=1000;
if size(wire_elements,2)>track_part_length
track_part_inds=[1:track_part_length:size(wire_elements,2) size(wire_elements,2)];
if track_part_inds(end-1)==track_part_inds(end)
track_part_inds(end)=[];
end
for parts_ind=1:numel(track_part_inds)-1
wire_part(parts_ind).coord=wire_elements(:,track_part_inds(parts_ind):track_part_inds(parts_ind+1));
wire_part(parts_ind).seg_coords=(wire_part(parts_ind).coord(:,1:end-1)+wire_part(parts_ind).coord(:,2:end))./2;
wire_part(parts_ind).currents=wire_part(parts_ind).coord(:,2:end)-wire_part(parts_ind).coord(:,1:end-1);
%wire_part(part_ind).seg_length=wire_elements(:,track_part_inds(part_ind):track_part_inds(part_ind+1));
end
else
wire_part.coord=wire_elements;
wire_part.seg_coords=(wire_part.coord(:,1:end-1)+wire_part.coord(:,2:end))./2;
wire_part.currents=wire_part.coord(:,2:end)-wire_part.coord(:,1:end-1);
end

%Calculate the magnetic field with Biot-Savarts law
b_field=zeros(3,num_tp);
for parts_ind=1:numel(wire_part)
target_p=repmat(target_f.coords,[1 1 size(wire_part(parts_ind).seg_coords,2)]); % copies target point coordinates
target_p=permute(target_p,[1 3 2]);
cur_pos= repmat(wire_part(parts_ind).seg_coords(:,:),[1 1 num_tp]); % copies vector of current position
cur_dir= repmat(wire_part(parts_ind).currents(:,:),[1 1 num_tp]); % copies vector of current direction
R = target_p-cur_pos;  % distance from current path to each point
%R(:,sqrt(dot(R,R,1))<min_distance_to_wire)=10^12; % if distance is to small the distance is blown up to ignore large values
len_R=(1./vecnorm(R,2,1)).^3; % distance factor of biot savart law
len_R=repmat(len_R,[3 1 1]);
dB=10^(-7)*cross(cur_dir,R,1).*len_R;
b_field(:,:)=b_field(:,:)+squeeze(sum(dB,2));
end
end
end


end

