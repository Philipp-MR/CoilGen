function local_opening_gab=calc_local_opening_gab(coil_mesh,middle_point,othorgonal_vector,opening_gab,vol_diagonal)

%Calulate the local opening width within the flat 2D domian


opening_gab=opening_gab/(1000000*vol_diagonal);
point_aa=middle_point+opening_gab*othorgonal_vector/2;
point_bb=middle_point-opening_gab*othorgonal_vector/2;
[point_a_curved,~]=uv_to_xyz(point_aa,triangulation(coil_mesh.faces',coil_mesh.uv'),triangulation(coil_mesh.faces',coil_mesh.vertices'));
[point_b_curved,~]=uv_to_xyz(point_bb,triangulation(coil_mesh.faces',coil_mesh.uv'),triangulation(coil_mesh.faces',coil_mesh.vertices'));
curved_length=vecnorm(point_b_curved-point_a_curved);
local_opening_gab=opening_gab*opening_gab/curved_length;
local_opening_gab=local_opening_gab*(1000000*vol_diagonal);


end