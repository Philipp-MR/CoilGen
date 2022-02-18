function single_ind_to_plot= find_even_leveled_solution(coil_layouts)
%Chose a even leveled solution for plotting
valid_layouts=arrayfun(@(x) ~isempty(coil_layouts(x).out),1:numel(coil_layouts)); 
levels=arrayfun(@(x) coil_layouts(x).out.input_data.levels,1:numel(coil_layouts));
offsets=arrayfun(@(x) coil_layouts(x).out.input_data.pot_offset_factor,1:numel(coil_layouts));
if any(~mod(levels,2).*valid_layouts)
single_ind_to_plot=find(~mod(levels,2).*valid_layouts);
single_ind_to_plot=single_ind_to_plot(ceil(numel(find(~mod(levels,2).*valid_layouts))/2));
else
single_ind_to_plot=find(valid_layouts);
single_ind_to_plot=single_ind_to_plot(ceil(numel(find(valid_layouts))/2));
end
end

