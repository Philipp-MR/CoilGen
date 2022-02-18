function layout_gradient= calculate_gradient(coil_parts,target_field,combined_field_layout_per1Amp,combined_field_loops_per1Amp,combined_mesh)
%CALCULATE_LOCAL_GRADIENT
%calculate for each point the average of the local gradients to its neighbours


mesh_diag_length=norm(combined_mesh.bounding_box(:,2)-combined_mesh.bounding_box(:,1));


layout_gradient.local_target_gx=zeros(1,size(target_field.b,2));
layout_gradient.local_target_gy=zeros(1,size(target_field.b,2));
layout_gradient.local_target_gz=zeros(1,size(target_field.b,2));
layout_gradient.local_gx=zeros(1,size(target_field.b,2));
layout_gradient.local_gy=zeros(1,size(target_field.b,2));
layout_gradient.local_gz=zeros(1,size(target_field.b,2));
layout_gradient.local_gx_loops=zeros(1,size(target_field.b,2));
layout_gradient.local_gy_loops=zeros(1,size(target_field.b,2));
layout_gradient.local_gz_loops=zeros(1,size(target_field.b,2));


% %define the neighbourhood in terms of a delaunay triangulation;
% target_Delaunay = delaunayTriangulation(target_field.coords(1,:)',target_field.coords(2,:)',target_field.coords(3,:)');
% %target_neighbors=nearestNeighbor(target_Delaunay,target_Delaunay.Points');
% target_neighbors=edges(target_Delaunay);
% target_neighbors_dists=vecnorm(target_Delaunay.Points(target_neighbors(:,2),:)-target_Delaunay.Points(target_neighbors(:,1),:),2,2);
% target_neighbors(target_neighbors_dists<mean(target_neighbors_dists))=[];

full_neighbourhood_points=[];
for target_group_ind=1:max(target_field.target_field_group_inds)
is_member=find(target_field.target_field_group_inds==target_group_ind);
%define the neighbourhood in terms of a delaunay triangulation;
target_Delaunay = delaunayTriangulation(target_field.coords(1,is_member)',target_field.coords(2,is_member)',target_field.coords(3,is_member)');
%target_neighbors=nearestNeighbor(target_Delaunay,target_Delaunay.Points');
target_neighbors=edges(target_Delaunay);
target_neighbors=[is_member(target_neighbors(:,1));is_member(target_neighbors(:,2))]';
full_neighbourhood_points=[full_neighbourhood_points; target_neighbors];
end
%target_neighbors_dists=vecnorm(target_field.coords(:,full_neighbourhood_points(:,2))-target_field.coords(:,full_neighbourhood_points(:,1)),2,1);
%target_neighbors(target_neighbors_dists<mean(target_neighbors_dists))=[];


for point_ind=1:size(target_field.coords,2)
    neighbour_inds=full_neighbourhood_points(full_neighbourhood_points(:,1)==point_ind,2);
    not_same_x=~(target_field.coords(1,point_ind)-target_field.coords(1,neighbour_inds)==0);
    not_same_y=~(target_field.coords(2,point_ind)-target_field.coords(2,neighbour_inds)==0);
    not_same_z=~(target_field.coords(3,point_ind)-target_field.coords(3,neighbour_inds)==0);
    
    denominator_x=(target_field.coords(1,neighbour_inds(not_same_x))-target_field.coords(1,point_ind));
    denominator_y=(target_field.coords(2,neighbour_inds(not_same_y))-target_field.coords(2,point_ind));
    denominator_z=(target_field.coords(3,neighbour_inds(not_same_z))-target_field.coords(3,point_ind));
    
    denominator_x(abs(denominator_x)<mesh_diag_length/100000)=nan; % ignore point pairs for which the xyz are very similar
    denominator_y(abs(denominator_y)<mesh_diag_length/100000)=nan;
    denominator_z(abs(denominator_z)<mesh_diag_length/100000)=nan;
    
	layout_gradient.local_target_gx(point_ind)=mean((target_field.b(3,neighbour_inds(not_same_x))-target_field.b(3,point_ind))./denominator_x);
	layout_gradient.local_target_gy(point_ind)=mean((target_field.b(3,neighbour_inds(not_same_y))-target_field.b(3,point_ind))./denominator_y);
	layout_gradient.local_target_gz(point_ind)=mean((target_field.b(3,neighbour_inds(not_same_z))-target_field.b(3,point_ind))./denominator_z);
	layout_gradient.local_gx(point_ind)=mean((combined_field_layout_per1Amp(3,neighbour_inds(not_same_x))-combined_field_layout_per1Amp(3,point_ind))./denominator_x);
	layout_gradient.local_gy(point_ind)=mean((combined_field_layout_per1Amp(3,neighbour_inds(not_same_y))-combined_field_layout_per1Amp(3,point_ind))./denominator_y);
	layout_gradient.local_gz(point_ind)=mean((combined_field_layout_per1Amp(3,neighbour_inds(not_same_z))-combined_field_layout_per1Amp(3,point_ind))./denominator_z);
    layout_gradient.local_gx_loops(point_ind)=mean((combined_field_loops_per1Amp(3,neighbour_inds(not_same_x))-combined_field_loops_per1Amp(3,point_ind))./denominator_x);
    layout_gradient.local_gy_loops(point_ind)=mean((combined_field_loops_per1Amp(3,neighbour_inds(not_same_y))-combined_field_loops_per1Amp(3,point_ind))./denominator_y);
    layout_gradient.local_gz_loops(point_ind)=mean((combined_field_loops_per1Amp(3,neighbour_inds(not_same_z))-combined_field_loops_per1Amp(3,point_ind))./denominator_z);
end
% %calculate gradient according the center system of the coordinate system
layout_gradient.target_gx=target_field.b(3,:)./target_field.coords(1,:); %~dBz/dx
layout_gradient.target_gy=target_field.b(3,:)./target_field.coords(2,:); %~dBz/dy
layout_gradient.target_gz=target_field.b(3,:)./target_field.coords(3,:); %~dBz/dz
layout_gradient.gx=combined_field_layout_per1Amp(3,:)./target_field.coords(1,:); %~dBz/dx
layout_gradient.gy=combined_field_layout_per1Amp(3,:)./target_field.coords(2,:); %~dBz/dy
layout_gradient.gz=combined_field_layout_per1Amp(3,:)./target_field.coords(3,:); %~dBz/dz
layout_gradient.gx_loops=combined_field_loops_per1Amp(3,:)./target_field.coords(1,:); %~dBz/dx
layout_gradient.gy_loops=combined_field_loops_per1Amp(3,:)./target_field.coords(2,:); %~dBz/dy
layout_gradient.gz_loops=combined_field_loops_per1Amp(3,:)./target_field.coords(3,:); %~dBz/dz
% Remove Inf values
layout_gradient.gx(~isfinite(layout_gradient.gx))=nan;
layout_gradient.gy(~isfinite(layout_gradient.gy))=nan;
layout_gradient.gz(~isfinite(layout_gradient.gz))=nan;
layout_gradient.gx_loops(~isfinite(layout_gradient.gx_loops))=nan;
layout_gradient.gy_loops(~isfinite(layout_gradient.gy_loops))=nan;
layout_gradient.gz_loops(~isfinite(layout_gradient.gz_loops))=nan;
layout_gradient.local_gx_loops(~isfinite(layout_gradient.local_gx_loops))=nan;
layout_gradient.local_gy_loops(~isfinite(layout_gradient.local_gy_loops))=nan;
layout_gradient.local_gz_loops(~isfinite(layout_gradient.local_gz_loops))=nan;
% go to the gradient unit mT/m/A
layout_gradient.target_gx=layout_gradient.target_gx.*1000;
layout_gradient.target_gy=layout_gradient.target_gy.*1000;
layout_gradient.target_gz=layout_gradient.target_gz.*1000;
layout_gradient.local_target_gx=layout_gradient.local_target_gx.*1000;
layout_gradient.local_target_gy=layout_gradient.local_target_gy.*1000;
layout_gradient.local_target_gz=layout_gradient.local_target_gz.*1000;
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
layout_gradient.mean_gx=mean(layout_gradient.gx,'omitnan');
layout_gradient.mean_gy=mean(layout_gradient.gy,'omitnan');
layout_gradient.mean_gz=mean(layout_gradient.gz,'omitnan');
layout_gradient.mean_gx_loops=mean(layout_gradient.gx_loops,'omitnan');
layout_gradient.mean_gy_loops=mean(layout_gradient.gy_loops,'omitnan');
layout_gradient.mean_gz_loops=mean(layout_gradient.gz_loops,'omitnan');
layout_gradient.std_gx=std(layout_gradient.gx,'omitnan');
layout_gradient.std_gy=std(layout_gradient.gy,'omitnan');
layout_gradient.std_gz=std(layout_gradient.gz,'omitnan');
layout_gradient.std_gx_loops=std(layout_gradient.gx_loops,'omitnan');
layout_gradient.std_gy_loops=std(layout_gradient.gy_loops,'omitnan');
layout_gradient.std_gz_loops=std(layout_gradient.gz_loops,'omitnan');
layout_gradient.mean_local_gx=mean(layout_gradient.local_gx,'omitnan');
layout_gradient.mean_local_gy=mean(layout_gradient.local_gy,'omitnan');
layout_gradient.mean_local_gz=mean(layout_gradient.local_gz,'omitnan');
layout_gradient.mean_local_gx_loops=mean(layout_gradient.local_gx_loops,'omitnan');
layout_gradient.mean_local_gy_loops=mean(layout_gradient.local_gy_loops,'omitnan');
layout_gradient.mean_local_gz_loops=mean(layout_gradient.local_gz_loops,'omitnan');
layout_gradient.std_local_gx=std(layout_gradient.local_gx,'omitnan');
layout_gradient.std_local_gy=std(layout_gradient.local_gy,'omitnan');
layout_gradient.std_local_gz=std(layout_gradient.local_gz,'omitnan');
layout_gradient.std_local_gx_loops=std(layout_gradient.local_gx_loops,'omitnan');
layout_gradient.std_local_gy_loops=std(layout_gradient.local_gy_loops,'omitnan');
layout_gradient.std_local_gz_loops=std(layout_gradient.local_gz_loops,'omitnan');


%Another idea would be to select the unique x,y,z values of the point cloud
%and build a meshgrid on it for which the gradient can then be regularly
%calculated

end

