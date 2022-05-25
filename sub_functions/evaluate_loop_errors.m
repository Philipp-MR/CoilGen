function [coil_parts,combined_field_layout,combined_field_loops,combined_field_layout_per1Amp,combined_field_loops_per1Amp,field_error_vals,opt_current_layout,layout_gradient]=evaluate_loop_errors(coil_parts,target_field,sf_b_field)
%Calcaltue relative errors between the different input and
%output fields


coil_parts(numel(coil_parts)).field_by_loops=[];

for part_ind=1:numel(coil_parts)
%Calcuate the combined field of the unconnected contours
coil_parts(part_ind).field_by_loops=zeros(3,size(target_field.b,2));
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
loop_field=biot_savart_calc_b(coil_parts(part_ind).contour_lines(loop_ind).v,target_field);
coil_parts(part_ind).field_by_loops=coil_parts(part_ind).field_by_loops+loop_field;
end
coil_parts(part_ind).field_by_loops=coil_parts(part_ind).field_by_loops.*coil_parts(part_ind).contour_step;
%Calculate the field of connected, final layouts
% Find the ideal current strength for the connected layout to match the target field
%coil_parts(part_ind).opt_current_layout=abs(mean(target_field.b(3,:)./coil_parts(part_ind).field_by_layout(3,:)));
end


%find the current polarity for the different coil parts regarding the
%target field
possible_polarities=cellfun(@(x) x=='1',num2cell(dec2bin(0:(2^(numel(coil_parts))-1))));
possible_polarities=double(possible_polarities);
possible_polarities(possible_polarities==0)=-1;


%Combine the coil_part fields scaled with all possible polarietes and chose the one
%which is closest to the target field
combined_field_loops=zeros([size(possible_polarities,1),size(coil_parts(part_ind).field_by_loops)]);
pol_projections_loops=zeros(1,size(possible_polarities,1));
for pol_ind=1:size(possible_polarities,1)
for part_ind=1:numel(coil_parts)
combined_field_loops(pol_ind,:,:)=squeeze(combined_field_loops(pol_ind,:,:))+possible_polarities(pol_ind,part_ind).*coil_parts(part_ind).field_by_loops;
end
%Project the combined field onto the target field
% pol_projections_layout(pol_ind)=sum(sum(squeeze(combined_field_layout(pol_ind,:,:)).*target_field.b,1));
% pol_projections_loops(pol_ind)=sum(sum(squeeze(combined_field_loops(pol_ind,:,:)).*target_field.b,1));
pol_projections_loops(pol_ind)=sum(vecnorm(squeeze(combined_field_loops(pol_ind,:,:))-target_field.b));
end

%Chose the best combination
[~,best_dir_loops]=min(pol_projections_loops);
combined_field_loops=squeeze(combined_field_loops(best_dir_loops,:,:));


%adjust the current direction for the loops (in case of wrong direction)
for part_ind=1:numel(coil_parts)
if possible_polarities(best_dir_loops,part_ind)~=1
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
coil_parts(part_ind).contour_lines(loop_ind).v=fliplr(coil_parts(part_ind).contour_lines(loop_ind).v);
coil_parts(part_ind).contour_lines(loop_ind).uv=fliplr(coil_parts(part_ind).contour_lines(loop_ind).uv);
end
end
end



%Only look at the z-component
target_z=target_field.b(3,:);
sf_z=sf_b_field(3,:); %field of stream function
loop_z=combined_field_loops(3,:);
field_error_vals.max_rel_error_layout_vs_target=0;
field_error_vals.mean_rel_error_layout_vs_target=0;
field_error_vals.max_rel_error_unconnected_contours_vs_target=max(abs((loop_z-target_z)./max(abs(target_z))),[],'all').*100;
field_error_vals.mean_rel_error_unconnected_contours_vs_target=mean(abs((loop_z-target_z)./max(abs(target_z))),'all').*100;
field_error_vals.max_rel_error_layout_vs_stream_function_field=0;
field_error_vals.mean_rel_error_layout_vs_stream_function_field=0;
field_error_vals.max_rel_error_unconnected_contours_vs_stream_function_field=max(abs((loop_z-sf_z)./max(abs(sf_z))),[],'all').*100;
field_error_vals.mean_rel_error_unconnected_contours_vs_stream_function_field=mean(abs((loop_z-sf_z)./max(abs(sf_z))),'all').*100;


%Go back to the fields for 1 Ampere (Unit Current)
combined_field_loops_per1Amp=zeros(size(combined_field_loops));
for part_ind=1:numel(coil_parts)
combined_field_loops_per1Amp=combined_field_loops_per1Amp+coil_parts(part_ind).field_by_loops./max(arrayfun(@(x) coil_parts(x).contour_step,1:numel(coil_parts))).*possible_polarities(best_dir_loops,part_ind);
end


%calculate for each point the average of the local gradients to its neighbours
layout_gradient.local_gx=zeros(1,size(target_field.b,2));
layout_gradient.local_gy=zeros(1,size(target_field.b,2));
layout_gradient.local_gz=zeros(1,size(target_field.b,2));
layout_gradient.local_gx_loops=zeros(1,size(target_field.b,2));
layout_gradient.local_gy_loops=zeros(1,size(target_field.b,2));
layout_gradient.local_gz_loops=zeros(1,size(target_field.b,2));
target_Delaunay = delaunayTriangulation(target_field.coords(1,:)',target_field.coords(2,:)',target_field.coords(3,:)');
%target_neighbors=nearestNeighbor(target_Delaunay,target_Delaunay.Points');
target_neighbors=edges(target_Delaunay);
for point_ind=1:size(target_field.coords,2)
    neighbour_inds=target_neighbors(target_neighbors(:,1)==point_ind,2);
    not_same_x=~(target_field.coords(1,point_ind)-target_field.coords(1,neighbour_inds)==0);
    not_same_y=~(target_field.coords(2,point_ind)-target_field.coords(2,neighbour_inds)==0);
    not_same_z=~(target_field.coords(3,point_ind)-target_field.coords(3,neighbour_inds)==0);
    layout_gradient.local_gx_loops(point_ind)=mean((combined_field_loops_per1Amp(3,neighbour_inds(not_same_x))-combined_field_loops_per1Amp(3,point_ind))./(target_field.coords(1,neighbour_inds(not_same_x))-target_field.coords(1,point_ind)));
    layout_gradient.local_gy_loops(point_ind)=mean((combined_field_loops_per1Amp(3,neighbour_inds(not_same_y))-combined_field_loops_per1Amp(3,point_ind))./(target_field.coords(2,neighbour_inds(not_same_y))-target_field.coords(2,point_ind)));
    layout_gradient.local_gz_loops(point_ind)=mean((combined_field_loops_per1Amp(3,neighbour_inds(not_same_z))-combined_field_loops_per1Amp(3,point_ind))./(target_field.coords(3,neighbour_inds(not_same_z))-target_field.coords(3,point_ind)));
end
% %calculate gradient according the center system of the coordinate system
layout_gradient.gx=zeros(1,size(target_field.b,2));
layout_gradient.gy=zeros(1,size(target_field.b,2));
layout_gradient.gz=zeros(1,size(target_field.b,2));
layout_gradient.gx_loops=combined_field_loops_per1Amp(3,:)./target_field.coords(1,:); %~dBz/dx
layout_gradient.gy_loops=combined_field_loops_per1Amp(3,:)./target_field.coords(2,:); %~dBz/dy
layout_gradient.gz_loops=combined_field_loops_per1Amp(3,:)./target_field.coords(3,:); %~dBz/dz
% Remove Inf values
layout_gradient.gx_loops(~isfinite(layout_gradient.gx_loops))=nan;
layout_gradient.gy_loops(~isfinite(layout_gradient.gy_loops))=nan;
layout_gradient.gz_loops(~isfinite(layout_gradient.gz_loops))=nan;
layout_gradient.local_gx_loops(~isfinite(layout_gradient.local_gx_loops))=nan;
layout_gradient.local_gy_loops(~isfinite(layout_gradient.local_gy_loops))=nan;
layout_gradient.local_gz_loops(~isfinite(layout_gradient.local_gz_loops))=nan;
% go to the gradient unit mT/m/A
layout_gradient.gx=layout_gradient.gx.*1000;
layout_gradient.gy=layout_gradient.gy.*1000;
layout_gradient.gz=layout_gradient.gz.*1000;
layout_gradient.gx_loops=layout_gradient.gx_loops.*1000;
layout_gradient.gy_loops=layout_gradient.gy_loops.*1000;
layout_gradient.gz_loops=layout_gradient.gz_loops.*1000;
layout_gradient.local_gx=layout_gradient.local_gx.*1000;
layout_gradient.local_gy=layout_gradient.local_gy.*1000;
layout_gradient.local_gz=layout_gradient.local_gz.*1000;
layout_gradient.local_gx_loops=layout_gradient.local_gx_loops.*1000;
layout_gradient.local_gy_loops=layout_gradient.local_gy_loops.*1000;
layout_gradient.local_gz_loops=layout_gradient.local_gz_loops.*1000;
%Calculate mean and standart deviation
layout_gradient.mean_gx=0;
layout_gradient.mean_gy=0;
layout_gradient.mean_gz=0;
layout_gradient.mean_gx_loops=mean(layout_gradient.gx_loops,'omitnan');
layout_gradient.mean_gy_loops=mean(layout_gradient.gy_loops,'omitnan');
layout_gradient.mean_gz_loops=mean(layout_gradient.gz_loops,'omitnan');
layout_gradient.std_gx=0;
layout_gradient.std_gy=0;
layout_gradient.std_gz=0;
layout_gradient.std_gx_loops=std(layout_gradient.gx_loops,'omitnan');
layout_gradient.std_gy_loops=std(layout_gradient.gy_loops,'omitnan');
layout_gradient.std_gz_loops=std(layout_gradient.gz_loops,'omitnan');
layout_gradient.mean_local_gx=0;
layout_gradient.mean_local_gy=0;
layout_gradient.mean_local_gz=0;
layout_gradient.mean_local_gx_loops=mean(layout_gradient.local_gx_loops,'omitnan');
layout_gradient.mean_local_gy_loops=mean(layout_gradient.local_gy_loops,'omitnan');
layout_gradient.mean_local_gz_loops=mean(layout_gradient.local_gz_loops,'omitnan');
layout_gradient.std_local_gx=0;
layout_gradient.std_local_gy=0;
layout_gradient.std_local_gz=0;
layout_gradient.std_local_gx_loops=std(layout_gradient.local_gx_loops,'omitnan');
layout_gradient.std_local_gy_loops=std(layout_gradient.local_gy_loops,'omitnan');
layout_gradient.std_local_gz_loops=std(layout_gradient.local_gz_loops,'omitnan');


combined_field_layout=zeros(size(combined_field_loops));
combined_field_layout_per1Amp=zeros(size(combined_field_loops));

% Find the ideal current strength for the connected layout to match the target field
opt_current_layout=0;

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