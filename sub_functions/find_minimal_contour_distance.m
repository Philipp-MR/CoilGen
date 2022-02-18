function coil_parts= find_minimal_contour_distance(coil_parts,pcb_track_width_factor)

%find the minimal distance in the xyz domian between contours to assign a proper conductor width later

coil_parts(numel(coil_parts)).pcb_track_width=[];

for part_ind=1:numel(coil_parts)


min_vals=[];

for ind_1=1:numel(coil_parts(part_ind).contour_lines)
    for ind_2=1:numel(coil_parts(part_ind).contour_lines)
        
        if ind_1~=ind_2
        
        loop_a=coil_parts(part_ind).contour_lines(ind_1).point_coordinates;
        loop_b=coil_parts(part_ind).contour_lines(ind_2).point_coordinates;
        
        loop_a_mat=repmat(loop_a,[1 1 size(loop_b,2)]);
        loop_b_mat=repmat(loop_b,[1 1 size(loop_a,2)]);
        
        dist_mat=((squeeze(loop_a_mat(1,:,:))'-squeeze(loop_b_mat(1,:,:))).^2+(squeeze(loop_a_mat(2,:,:))'-squeeze(loop_b_mat(2,:,:))).^2+(squeeze(loop_a_mat(3,:,:))'-squeeze(loop_b_mat(3,:,:))).^2).^(1/2);
        
        min_vals=[min_vals min(dist_mat(:))];
        
        end
        
    end
end

coil_parts(part_ind).pcb_track_width=pcb_track_width_factor*min(min_vals)/2;


end

end

