function [coil_parts,combined_field_layout,combined_field_loops,combined_field_layout_per1Amp,combined_field_loops_per1Amp,field_error_vals,opt_current_layout]=evaluate_field_errors(coil_parts,input,target_field,sf_b_field)
%Calcaltue relative errors between the different input and
%output fields




coil_parts(numel(coil_parts)).opt_current_layout=[];
coil_parts(numel(coil_parts)).field_by_loops=[];
coil_parts(numel(coil_parts)).field_by_layout=[];


for part_ind=1:numel(coil_parts)
%Calcuate the combined field of the unconnected contours
coil_parts(part_ind).field_by_loops=zeros(3,size(target_field.b,2));
for loop_ind=1:numel(coil_parts(part_ind).contour_lines)
loop_field=biot_savart_calc_b(coil_parts(part_ind).contour_lines(loop_ind).v,target_field);
coil_parts(part_ind).field_by_loops=coil_parts(part_ind).field_by_loops+loop_field;
end
coil_parts(part_ind).field_by_loops=coil_parts(part_ind).field_by_loops.*coil_parts(part_ind).contour_step;
%Calculate the field of connected, final layouts
if ~input.skip_postprocessing
coil_parts(part_ind).field_by_layout=biot_savart_calc_b(coil_parts(part_ind).wire_path.v,target_field);
coil_parts(part_ind).field_by_layout=coil_parts(part_ind).field_by_layout.*coil_parts(part_ind).contour_step; %scaled with the current of the discretization
else
coil_parts(part_ind).field_by_layout=coil_parts(part_ind).field_by_loops;
end
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
combined_field_layout=zeros([size(possible_polarities,1),size(coil_parts(part_ind).field_by_layout)]);
combined_field_loops=zeros([size(possible_polarities,1),size(coil_parts(part_ind).field_by_loops)]);
pol_projections_loops=zeros(1,size(possible_polarities,1));
pol_projections_layout=zeros(1,size(possible_polarities,1));
for pol_ind=1:size(possible_polarities,1)
for part_ind=1:numel(coil_parts)
combined_field_layout(pol_ind,:,:)=squeeze(combined_field_layout(pol_ind,:,:))+possible_polarities(pol_ind,part_ind).*coil_parts(part_ind).field_by_layout;
combined_field_loops(pol_ind,:,:)=squeeze(combined_field_loops(pol_ind,:,:))+possible_polarities(pol_ind,part_ind).*coil_parts(part_ind).field_by_loops;
end
%Project the combined field onto the target field
% pol_projections_layout(pol_ind)=sum(sum(squeeze(combined_field_layout(pol_ind,:,:)).*target_field.b,1));
% pol_projections_loops(pol_ind)=sum(sum(squeeze(combined_field_loops(pol_ind,:,:)).*target_field.b,1));
pol_projections_layout(pol_ind)=sum(vecnorm(squeeze(combined_field_layout(pol_ind,:,:))-target_field.b));
pol_projections_loops(pol_ind)=sum(vecnorm(squeeze(combined_field_loops(pol_ind,:,:))-target_field.b));
end

%Chose the best combination
[~,best_dir_layout]=min(pol_projections_layout);
[~,best_dir_loops]=min(pol_projections_loops);
combined_field_layout=squeeze(combined_field_layout(best_dir_layout,:,:));
combined_field_loops=squeeze(combined_field_loops(best_dir_loops,:,:));


%adjust the current direction for the layout (in case of wrong direction)
if ~input.skip_postprocessing
for part_ind=1:numel(coil_parts)
if possible_polarities(best_dir_layout,part_ind)~=1
coil_parts(part_ind).wire_path.v=fliplr(coil_parts(part_ind).wire_path.v);
coil_parts(part_ind).wire_path.uv=fliplr(coil_parts(part_ind).wire_path.uv);
end
end
end

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
layout_z=combined_field_layout(3,:);
loop_z=combined_field_loops(3,:);
field_error_vals.max_rel_error_layout_vs_target=max(abs((layout_z-target_z)./max(abs(target_z))),[],'all').*100;
field_error_vals.mean_rel_error_layout_vs_target=mean(abs((layout_z-target_z)./max(abs(target_z))),'all').*100;
field_error_vals.max_rel_error_unconnected_contours_vs_target=max(abs((loop_z-target_z)./max(abs(target_z))),[],'all').*100;
field_error_vals.mean_rel_error_unconnected_contours_vs_target=mean(abs((loop_z-target_z)./max(abs(target_z))),'all').*100;
field_error_vals.max_rel_error_layout_vs_stream_function_field=max(abs((layout_z-sf_z)./max(abs(sf_z))),[],'all').*100;
field_error_vals.mean_rel_error_layout_vs_stream_function_field=mean(abs((layout_z-sf_z)./max(abs(sf_z))),'all').*100;
field_error_vals.max_rel_error_unconnected_contours_vs_stream_function_field=max(abs((loop_z-sf_z)./max(abs(sf_z))),[],'all').*100;
field_error_vals.mean_rel_error_unconnected_contours_vs_stream_function_field=mean(abs((loop_z-sf_z)./max(abs(sf_z))),'all').*100;

%Go back to the fields for 1 Ampere (Unit Current)
combined_field_layout_per1Amp=zeros(size(combined_field_layout));
combined_field_loops_per1Amp=zeros(size(combined_field_layout));
for part_ind=1:numel(coil_parts)
combined_field_layout_per1Amp=combined_field_layout_per1Amp+coil_parts(part_ind).field_by_layout./max(arrayfun(@(x) coil_parts(x).contour_step,1:numel(coil_parts))).*possible_polarities(best_dir_layout,part_ind);
combined_field_loops_per1Amp=combined_field_loops_per1Amp+coil_parts(part_ind).field_by_loops./max(arrayfun(@(x) coil_parts(x).contour_step,1:numel(coil_parts))).*possible_polarities(best_dir_loops,part_ind);
end




% Find the ideal current strength for the connected layout to match the target field
opt_current_layout=abs(mean(target_field.b(3,:)./combined_field_layout(3,:)));


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