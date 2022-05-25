function cut_rectangle=build_cut_rectangle(loop,center_point,segment_ind,cut_witdh,cut_height_ratio)
%build a rectangular cut shape


cut_points_left=loop.uv(:,segment_ind+1);
cut_points_right=loop.uv(:,segment_ind);

longitudinal_vector=cut_points_left-cut_points_right;
longitudinal_vector=longitudinal_vector./vecnorm(longitudinal_vector);
othorgonal_vector=[longitudinal_vector(2,:); longitudinal_vector(1,:)*(-1)];
othorgonal_vector=othorgonal_vector./(vecnorm(othorgonal_vector));

%Scale the points to the targeted width
othorgonal_vector=othorgonal_vector.*(cut_witdh.*cut_height_ratio);
longitudinal_vector=longitudinal_vector.*cut_witdh;

%Create the rectangular points
cut_rectangle=[center_point+longitudinal_vector/2+othorgonal_vector/2 ...
                                center_point+longitudinal_vector/2-othorgonal_vector/2 ...
                                center_point-longitudinal_vector/2-othorgonal_vector/2 ...
                                center_point-longitudinal_vector/2+othorgonal_vector/2 ...
                                center_point+longitudinal_vector/2+othorgonal_vector/2];

end