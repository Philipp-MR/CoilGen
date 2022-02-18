function [field_by_loops,relative_max_field_error,relative_mean_field_error]=evaluate_resulting_field_pcb(pcb_layout,target_field,potential_step,plot_flag)



target_field.Bz=target_field.Bz';
    

field_by_loops=zeros(3,size(target_field,2));
for track_ind=1:numel(pcb_layout.coil_track)
loop_field1=biot_savart_calc_b(pcb_layout.coil_track(track_ind).spiral_in.v,target_field.coords);
loop_field2=biot_savart_calc_b(pcb_layout.coil_track(track_ind).spiral_out.v,target_field.coords);
field_by_loops=field_by_loops+loop_field1+loop_field2;
end

% layout_c=field_by_loops(3,:)./mean(abs(field_by_loops(3,:))); %Normalize for comparison
% target_c=target_field.bz./mean(abs(target_field.bz)); %Normalize for comparison

% layout_c=field_by_loops(3,:)./max(abs(field_by_loops(3,:))); %Normalize for comparison
% target_c=target_field.bz./max(abs(target_field.bz)); %Normalize for comparison

field_by_loops=field_by_loops.*potential_step/2; %assign the current strength by the ptential step

layout_c=field_by_loops(3,:);
target_c=target_field.Bz;


if size(layout_c)~=size(target_c)
target_c=target_c';
end

relative_max_field_error=max(abs((layout_c-target_c)./max(abs(target_c))),[],'all').*100;
relative_mean_field_error=mean(abs((layout_c-target_c)./max(abs(target_c)))).*100;


%adjust the current direction for the layout
if relative_max_field_error>100
field_by_loops=field_by_loops.*(-1);
layout_c=field_by_loops(3,:)./max(abs(field_by_loops(3,:)));
target_c=target_field.Bz./max(abs(target_field.Bz));
if size(layout_c)~=size(target_c)
target_c=target_c';
end
relative_max_field_error=max(abs((layout_c-target_c)./max(abs(target_c))),[],'all').*100;
relative_mean_field_error=mean(abs((layout_c-target_c)./max(abs(target_c)))).*100;
end



function b_field=biot_savart_calc_b(wire_elements,target_coords)
    
%to avoid memory problems splitt the curve into several parts
track_part_length=1000;
    
    
    
    num_tp=size(target_coords,2);

        %to avoid memory problems splitt the curve into several parts
        
        if size(wire_elements,2)<track_part_length
            wire_part.coord=wire_elements;
            wire_part.seg_coords=(wire_part.coord(:,1:end-1)+wire_part.coord(:,2:end))./2;
            wire_part.currents=wire_part.coord(:,2:end)-wire_part.coord(:,1:end-1);
        else
        track_part_inds=[1:track_part_length:size(wire_elements,2) size(wire_elements,2)];
        for part_ind=1:numel(track_part_inds)-1
            wire_part(part_ind).coord=wire_elements(:,track_part_inds(part_ind):track_part_inds(part_ind+1));
            wire_part(part_ind).seg_coords=(wire_part(part_ind).coord(:,1:end-1)+wire_part(part_ind).coord(:,2:end))./2;
            wire_part(part_ind).currents=wire_part(part_ind).coord(:,2:end)-wire_part(part_ind).coord(:,1:end-1);
        end
        end
        
        %Calculate the magnetic field with Biot-Savarts law
        b_field=zeros(3,num_tp);
        for part_ind=1:numel(wire_part)
                target_p=repmat(target_coords,[1 1 size(wire_part(part_ind).seg_coords,2)]); % copies target point coordinates
                target_p=permute(target_p,[1 3 2]);
                cur_pos= repmat(wire_part(part_ind).seg_coords(:,:),[1 1 num_tp]); % copies vector of current position
                cur_dir= repmat(wire_part(part_ind).currents(:,:),[1 1 num_tp]); % copies vector of current direction
                R = target_p-cur_pos;  % distance from current path to each point
                %R(:,sqrt(dot(R,R,1))<min_distance_to_wire)=10^12; % if distance is to small the distance is blown up to ignore large values
                len_R=(1./vecnorm(R,2,1)).^3; % distance factor of biot savart law
                len_R=repmat(len_R,[3 1 1]);
                dB=10^(-7)*cross(cur_dir,R,1).*len_R;
                b_field(:,:)=b_field(:,:)+squeeze(sum(dB,2));
        end

end

% % % %Plot the field results
% % % if plot_flag
% % % 
% % % 
% % % figure;
% % % t = tiledlayout('flow');
% % % title(t,"Max. Relative Error:"+" "+num2str(relative_max_field_error  ,3)+"%"+", "+"Mean Relative Error:"+" "+num2str(relative_mean_field_error,3) +"%" );
% % % nexttile;
% % % hold on; 
% % % title('B-Field Layout');
% % % axis equal;
% % % for track_ind=1:numel(pcb_layout.coil_track)
% % % plot3(pcb_layout.coil_track(track_ind).spiral_in.v(1,:),pcb_layout.coil_track(track_ind).spiral_in.v(2,:),pcb_layout.coil_track(track_ind).spiral_in.v(3,:));
% % % plot3(pcb_layout.coil_track(track_ind).spiral_out.v(1,:),pcb_layout.coil_track(track_ind).spiral_out.v(2,:),pcb_layout.coil_track(track_ind).spiral_out.v(3,:));
% % % end
% % % 
% % % scatter3(target_field.coords(1,:),target_field.coords(2,:),target_field.coords(3,:),[],layout_c);
% % % %trisurf(triangulation(result_coil(plot_ind).out.parameterized_mesh .f,result_coil(plot_ind).out.parameterized_mesh .v),'facecolor',[1 1 1],'facealpha',0.5,'edgecolor',[0.8 0.8 0.8],'edgealpha',0.5);
% % % view(45,45);
% % % xlabel('x[m]');
% % % ylabel('y[m]');
% % % zlabel('z[m]');
% % % hold off;
% % % nexttile;
% % % hold on; 
% % % title('Layout field Histo');
% % % histogram(layout_c);
% % % hold off;
% % % nexttile;
% % % hold on; axis equal; title('Target field');
% % % scatter3(target_field.coords(1,:),target_field.coords(2,:),target_field.coords(3,:),[],target_c);
% % % view(45,45);
% % % xlabel('x[m]');
% % % ylabel('y[m]');
% % % zlabel('z[m]');
% % % hold off;
% % % nexttile;
% % % hold on; 
% % % title('Target Histo');
% % % histogram(target_c);
% % % hold off;
% % % 
% % % 
% % % set(gcf,'color','w');
% % % set(gcf, 'Position',  [1000, 100, 1000, 1000])
% % % 
% % % 
% % % end



end