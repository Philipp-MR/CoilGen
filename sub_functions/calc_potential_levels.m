function [coil_parts,primary_surface_ind] = calc_potential_levels(coil_parts,combined_mesh,input)
%center the stream function potential around zero and add zeros around the
%periphery

num_levels=input.levels;
level_offset=input.pot_offset_factor;
level_set_method=input.level_set_method;


switch level_set_method
    
    case "primary" %select the stream function of the part with the highest current density for the finding level width
        
%for the general case of multiple independent surfaces; the expression
%"num_levels" defines the number of contour steps within the surface part
%of the biggest range of stream function values i.e. the most dense current
%density
sf_range_per_part=zeros(1,numel(coil_parts));
for part_ind=1:numel(coil_parts)
sf_range_per_part(part_ind)=max(coil_parts(part_ind).stream_function)-min(coil_parts(part_ind).stream_function);
end
[~,primary_surface_ind]=max(sf_range_per_part);
%calc the contour step width (the later current strength accordingly)
contour_step=sf_range_per_part(primary_surface_ind)./(num_levels-1+2*level_offset);
%Set for each coil part the set of potentials with the same overal step
%width and set them within their potential boundaries
coil_parts(numel(coil_parts)).potential_level_list=[];
coil_parts(numel(coil_parts)).contour_step=[];
for part_ind=1:numel(coil_parts)
%create the potential level list by number of steps and offset ratio
coil_parts(part_ind).contour_step=contour_step;
if part_ind==primary_surface_ind % use the precise level offset for the primary part
%num_pot_steps=floor(sf_range_per_part(part_ind)/contour_step); %Calculate for each part the number of turns;
%coil_parts(part_ind).potential_level_list=([1:num_pot_steps(part_ind)]-(1-level_offset))*contour_step+min(coil_parts(part_ind).stream_function);
coil_parts(part_ind).potential_level_list=[0:num_levels-1].*contour_step+(min(coil_parts(part_ind).stream_function)+level_offset*contour_step);
else % center the contours regularly for the other parts around its mean potential
pot_residual=sf_range_per_part(part_ind)-2*level_offset*contour_step;
num_pot_steps=floor(pot_residual/contour_step); %Calculate for each part the number of turns;
coil_parts(part_ind).potential_level_list=[0:num_pot_steps-1].*contour_step+(min(coil_parts(part_ind).stream_function)+level_offset*contour_step);
%Center the levels for even dist to the potential maximum and minimum
dist_to_pot_max=max(coil_parts(part_ind).stream_function)-coil_parts(part_ind).potential_level_list(end);
dist_to_pot_min=coil_parts(part_ind).potential_level_list(1)-min(coil_parts(part_ind).stream_function);
if dist_to_pot_max<dist_to_pot_min
   coil_parts(part_ind).potential_level_list=coil_parts(part_ind).potential_level_list-(dist_to_pot_max-dist_to_pot_min)/2;
else
   coil_parts(part_ind).potential_level_list=coil_parts(part_ind).potential_level_list+(dist_to_pot_max-dist_to_pot_min)/2;
end
end
end
        
    case "combined" %select combined stream function finding level set
        
for part_ind=1:numel(coil_parts)
coil_parts(part_ind).contour_step=(max(combined_mesh.stream_function)-min(combined_mesh.stream_function))./(num_levels-1+2*level_offset);
coil_parts(part_ind).potential_level_list=([1:num_levels]-(1-level_offset))*coil_parts(part_ind).contour_step+min(combined_mesh.stream_function);
end
primary_surface_ind=1;

    case "independent" % calculate the level set inepentendtly for each part
        
%Calculate for each part the set of levels independently
if level_set_method
for part_ind=1:numel(coil_parts)
coil_parts(part_ind).contour_step=(max(coil_parts(part_ind).stream_function)-min(coil_parts(part_ind).stream_function))./(num_levels-1+2*level_offset);
coil_parts(part_ind).potential_level_list=([1:num_levels]-(1-level_offset))*coil_parts(part_ind).contour_step+min(coil_parts(part_ind).stream_function  );
end
end
primary_surface_ind=1;
        
end




% %Delete level which is directly one the overall potential mean can lead
% %to non-symmetrical layouts
% if remove_level_of_avg_potential
% for part_ind=1:numel(coil_parts)
% mean_pot=(max(coil_parts(part_ind).stream_function)+min(coil_parts(part_ind).stream_function))/2;
% is_avg_pot_level=abs(coil_parts(part_ind).potential_level_list-mean_pot)./max(coil_parts(part_ind).stream_function)<0.01;
% coil_parts(part_ind).potential_level_list(is_avg_pot_level)=[];
% end
% end


% %Plot the levels
% figure;
% hold on;
% color_mat=jet(numel(coil_parts));
% for part_ind=1:numel(coil_parts)
% plot([0 1],ones(1,2).*min(coil_parts(part_ind).stream_function),'color',color_mat(part_ind,:),'linewidth',5);
% plot([0 1],ones(1,2).*max(coil_parts(part_ind).stream_function),'color',color_mat(part_ind,:),'linewidth',5);
% plot([0 1],ones(1,2).*(max(coil_parts(part_ind).stream_function)+min(coil_parts(part_ind).stream_function))/2,'-','color',color_mat(part_ind,:),'linewidth',5);
% end
% for part_ind=1:numel(coil_parts)
%     for level_ind=1:numel(coil_parts(part_ind).potential_level_list)
%         plot([0 1],ones(1,2).*coil_parts(part_ind).potential_level_list(level_ind),'color',color_mat(part_ind,:));
%     end
% end

% %create the potential level list by number of steps and offset ratio
% contour_step=(max(combined_mesh.stream_function)-min(combined_mesh.stream_function))./(num_levels-1+2*level_offset);
% for part_ind=1:numel(coil_parts)
% coil_parts(part_ind).contour_step=contour_step;
% coil_parts(part_ind).potential_level_list=([1:num_levels]-(1-level_offset))*coil_parts(part_ind).contour_step+min(coil_parts(part_ind).stream_function  );
% end

end