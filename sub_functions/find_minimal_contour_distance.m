function coil_parts= find_minimal_contour_distance(coil_parts,input)
%find the minimal distance in the xyz domian between contours to assign a proper conductor width later


coil_parts(numel(coil_parts)).pcb_track_width=[];

for part_ind=1:numel(coil_parts)
if ~input.skip_calculation_min_winding_distance

min_vals=[];

for ind_1=1:numel(coil_parts(part_ind).contour_lines)
    for ind_2=1:numel(coil_parts(part_ind).contour_lines)
        if ind_1~=ind_2
        [min_dist,~,~,~,~]=find_min_mutual_loop_distance(coil_parts(part_ind).contour_lines(ind_1),coil_parts(part_ind).contour_lines(ind_2),false);
        min_vals=[min_vals min_dist];
        end
    end
end
coil_parts(part_ind).pcb_track_width=min(min_vals);

else
coil_parts(part_ind).pcb_track_width=0;
end

end

end

