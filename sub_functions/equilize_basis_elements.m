function equilized_basis_elements=equilize_basis_elements(basis_elements)
%fill the basis elements with zeros so that every elements has the same
%number of triangles and further calulations can be vectorized

equilized_basis_elements=basis_elements;
equilized_basis_elements(1).is_real_triangle=[];

highst_triangle_count_per_node=max(arrayfun(@(x) numel(basis_elements(x).area),1:numel(basis_elements)));

for hhhh=1:numel(basis_elements)
    
equilized_basis_elements(hhhh).is_real_triangle=zeros(1,highst_triangle_count_per_node);
equilized_basis_elements(hhhh).is_real_triangle(1:numel(basis_elements(hhhh).area))=1;

equilized_basis_elements(hhhh).triangle_points_ABC=zeros(highst_triangle_count_per_node,3,3);
equilized_basis_elements(hhhh).triangle_points_ABC(1:numel(basis_elements(hhhh).area),:,:)=basis_elements(hhhh).triangle_points_ABC;

equilized_basis_elements(hhhh).face_normal=zeros(highst_triangle_count_per_node,3);
equilized_basis_elements(hhhh).face_normal(1:numel(basis_elements(hhhh).area),:)=basis_elements(hhhh).face_normal;

equilized_basis_elements(hhhh).current=zeros(highst_triangle_count_per_node,3);
equilized_basis_elements(hhhh).current(1:numel(basis_elements(hhhh).area),:)=basis_elements(hhhh).current;

equilized_basis_elements(hhhh).area=zeros(highst_triangle_count_per_node,1);
equilized_basis_elements(hhhh).area(1:numel(basis_elements(hhhh).area))=basis_elements(hhhh).area;

end


end