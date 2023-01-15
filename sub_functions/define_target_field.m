function [target_field_out,is_supressed_point]= define_target_field(coil_parts,target_mesh,secondary_target_mesh,input)



%Define the target field
if ~strcmp(input.target_field_definition_file,'none')

if ispc
loaded_target_field=load(cd+"\"+"target_fields"+"\"+input.target_field_definition_file);
else
loaded_target_field=load(cd+"/"+"target_fields"+"/"+input.target_field_definition_file);
end
struct_name=fieldnames(loaded_target_field);
loaded_target_field=getfield(loaded_target_field,struct_name{1});
if isfield(loaded_target_field,input.target_field_definition_field_name)
loaded_field=getfield(loaded_target_field,input.target_field_definition_field_name);
if size(loaded_field,1)==1
target_field_out.b=[zeros(size(loaded_field));zeros(size(loaded_field));loaded_field];
else
target_field_out.b=loaded_field;
end
is_supressed_point=zeros(1,size(target_field_out.b,2));
target_field_out.coords=loaded_target_field.coords;
target_field_out.weights=ones(size(target_field_out.b));
target_field_out.target_field_group_inds=ones(1,size(target_field_out.b,2));
else
error("The target field with name"+" "+input.target_field_definition_file+" "+"does not exist, in the provided file");
end


else

if ~isempty(target_mesh) %create evenly distributed points within the suface of the "target mesh"
    
    
if ~input.use_only_target_mesh_verts
    
target_mesh_x_bounds=[min(target_mesh.vertices(:,1)) max(target_mesh.vertices(:,1))];
target_mesh_y_bounds=[min(target_mesh.vertices(:,2)) max(target_mesh.vertices(:,2))];
target_mesh_z_bounds=[min(target_mesh.vertices(:,3)) max(target_mesh.vertices(:,3))];

x_size=target_mesh_x_bounds(2)-target_mesh_x_bounds(1);
y_size=target_mesh_y_bounds(2)-target_mesh_y_bounds(1);
z_size=target_mesh_z_bounds(2)-target_mesh_z_bounds(1);

% num_points_per_x_dim=ceil(x_size/target_point_resolution);
% num_points_per_y_dim=ceil(y_size/target_point_resolution);
% num_points_per_z_dim=ceil(z_size/target_point_resolution);

num_points_per_x_dim=input.target_region_resolution;
num_points_per_y_dim=input.target_region_resolution;
num_points_per_z_dim=input.target_region_resolution;

target_x_coords=[target_mesh_x_bounds(1):(x_size/(num_points_per_x_dim-1)):target_mesh_x_bounds(2)];
target_y_coords=[target_mesh_y_bounds(1):(y_size/(num_points_per_y_dim-1)):target_mesh_y_bounds(2)];
target_z_coords=[target_mesh_z_bounds(1):(z_size/(num_points_per_z_dim-1)):target_mesh_z_bounds(2)];

[target_grid_x,target_grid_y,target_grid_z]=meshgrid(target_x_coords,target_y_coords,target_z_coords);
target_points=[target_grid_x(:) target_grid_y(:)  target_grid_z(:)]';

%Remove the points which are not inside the target surface
in = intriangulation((target_mesh.vertices-mean(target_mesh.vertices)).*0.9+mean(target_mesh.vertices),target_mesh.faces,target_points');

target_points=target_points(:,in);
target_points=[target_mesh.vertices' target_points]; %Add the surface vertices from the target mesh


else
    
target_points=target_mesh.vertices';  
end


else %Define the target point coordintes as points inside a sphere of a given radius
%num_points_per_dim=ceil(2*input.target_region_radius./target_point_resolution);
num_points_per_dim=input.target_region_resolution;
target_x_coords=[-2:(1/(num_points_per_dim-1)):2].*input.target_region_radius;
target_y_coords=[-2:(1/(num_points_per_dim-1)):2].*input.target_region_radius;
target_z_coords=[-2:(1/(num_points_per_dim-1)):2].*input.target_region_radius;
[target_grid_x,target_grid_y,target_grid_z]=meshgrid(target_x_coords,target_y_coords,target_z_coords);
target_points=[target_grid_x(:) target_grid_y(:)  target_grid_z(:)]';
%select the points that are inside a sphere 
target_points(:,sqrt(target_points(1,:).^2+target_points(2,:).^2+target_points(3,:).^2)>input.target_region_radius)=[];
all_verts=[];
for part_ind=1:numel(coil_parts)
all_verts=[all_verts coil_parts(part_ind).coil_mesh.vertices];
end
if input.set_roi_into_mesh_center
target_points=target_points-mean(all_verts,2);
end

end

%Remove identical points
[~,unique_inds,~] = unique(target_points','stable','rows');
target_points=target_points(:,unique_inds);


%define the target field shape
field_func=str2func("@(x,y,z)"+" "+input.field_shape_function);
target_field=zeros(size(target_points));
target_field(3,:)=field_func(target_points(1,:),target_points(2,:),target_points(3,:));



%Add points where the magnetic field should be supressed (=>0) 
if ~isempty(secondary_target_mesh)
num_supressed_points=size(secondary_target_mesh.vertices',2);
target_points=[target_points secondary_target_mesh.vertices'];
target_field=[target_field zeros(size(secondary_target_mesh.vertices'))];
is_supressed_point=false(1,size(target_points,2));
is_supressed_point((size(target_points,2)-num_supressed_points):size(target_points,2))=true;
else
is_supressed_point=false(1,size(target_points,2));
end

%Scale the fields to a targeted strength 
max_field_point_ind=find(target_field(3,:)==max(target_field(3,:)));
min_field_point_ind=find(target_field(3,:)==min(target_field(3,:)));
max_target_distance=vecnorm(target_points(:,max_field_point_ind(1))-target_points(:,min_field_point_ind(1)));
max_field_difference=max(target_field(3,:))-min(target_field(3,:));
if abs(max_field_difference)>10^(-10)
target_field=target_field./(max_field_difference/max_target_distance(1)).*input.target_gradient_strength;
else
target_field=target_field.*input.target_gradient_strength;
end

%Define weigthins from 0 to 1 that weights the significant of target points
target_field_weighting=ones(1,size(target_field,2));
target_field_weighting(is_supressed_point)=input.secondary_target_weight;
target_field_group_inds=ones(1,size(target_field,2));
target_field_group_inds(is_supressed_point)=2;

% %Permute also the target field
% perm_inds=randperm(size(target_points,2));
% target_points=target_points(:,perm_inds);
% target_field=target_field(:,perm_inds);
% target_field_weighting=target_field_weighting(perm_inds);
% target_field_group_inds=target_field_group_inds(perm_inds);

%Calcuate the gradients from the symbolic definition of the target field
[target_dbzbx,target_dbzby,target_dbzbz]=symbolic_calucation_of_gradient(input,target_field);

target_field_out.b=target_field;
target_field_out.coords=target_points;
target_field_out.weights=target_field_weighting;
target_field_out.target_field_group_inds=target_field_group_inds;
target_field_out.target_gradient_dbdxyz=[target_dbzbx; target_dbzby; target_dbzbz];



end

end

function [target_dbzbx,target_dbzby,target_dbzbz]=symbolic_calucation_of_gradient(input,target_field)
%Calcuate the gradients from the symbolic definition of the target field
try
syms x y z
dbzdx_fun =  string(diff(str2sym(input.field_shape_function),x));
dbzdy_fun  = string(diff(str2sym(input.field_shape_function),y));
dbzdz_fun  = string(diff(str2sym(input.field_shape_function),z));
%change the function handle for array-wise 
dbzdx_fun=replace(dbzdx_fun,"/","./");
dbzdx_fun=replace(dbzdx_fun,"^",".^");
dbzdx_fun=replace(dbzdx_fun,"*",".*");
dbzdy_fun=replace(dbzdy_fun,"/","./");
dbzdy_fun=replace(dbzdy_fun,"^",".^");
dbzdy_fun=replace(dbzdy_fun,"*",".*");
dbzdz_fun=replace(dbzdz_fun,"/","./");
dbzdz_fun=replace(dbzdz_fun,"^",".^");
dbzdz_fun=replace(dbzdz_fun,"*",".*");
dbzdx_fun=str2func("@(x,y,z)"+" "+dbzdx_fun);
dbzdy_fun=str2func("@(x,y,z)"+" "+dbzdy_fun);
dbzdz_fun=str2func("@(x,y,z)"+" "+dbzdz_fun);
target_dbzbx=dbzdx_fun(target_field(1,:),target_field(2,:),target_field(3,:));
target_dbzby=dbzdy_fun(target_field(1,:),target_field(2,:),target_field(3,:));
target_dbzbz=dbzdz_fun(target_field(1,:),target_field(2,:),target_field(3,:));
if numel(target_dbzbx)==1
target_dbzbx=repelem(target_dbzbx,size(target_field,2));
end
if numel(target_dbzby)==1
target_dbzby=repelem(target_dbzby,size(target_field,2));
end
if numel(target_dbzbz)==1
target_dbzbz=repelem(target_dbzbz,size(target_field,2));
end
catch
target_dbzbx=zeros(size(target_field(3,:)));
target_dbzby=zeros(size(target_field(3,:)));
target_dbzbz=zeros(size(target_field(3,:)));
disp('Gradient Calcuation from Symbolic Target failed');
end
end


function in = intriangulation(vertices,faces,testp,heavytest)
% intriangulation: Test points in 3d wether inside or outside a (closed) triangulation
% usage: in = intriangulation(vertices,faces,testp,heavytest)
%
% arguments: (input)
%  vertices   - points in 3d as matrix with three columns
%
%  faces      - description of triangles as matrix with three columns.
%               Each row contains three indices into the matrix of vertices
%               which gives the three cornerpoints of the triangle.
%
%  testp      - points in 3d as matrix with three columns
%
%  heavytest  - int n >= 0. Perform n additional randomized rotation tests.
%
% IMPORTANT: the set of vertices and faces has to form a watertight surface!
%
% arguments: (output)
%  in - a vector of length size(testp,1), containing 0 and 1.
%       in(nr) =  0: testp(nr,:) is outside the triangulation
%       in(nr) =  1: testp(nr,:) is inside the triangulation
%       in(nr) = -1: unable to decide for testp(nr,:) 
%
% Thanks to Adam A for providing the FEX submission voxelise. The
% algorithms of voxelise form the algorithmic kernel of intriangulation.
%
% Thanks to Sven to discussions about speed and avoiding problems in
% special cases.
%
% Example usage:
%
%      n = 10;
%      vertices = rand(n, 3)-0.5; % Generate random points
%      tetra = delaunayn(vertices); % Generate delaunay triangulization
%      faces = freeBoundary(TriRep(tetra,vertices)); % use free boundary as triangulation
%      n = 1000;
%      testp = 2*rand(n,3)-1; % Generate random testpoints
%      in = intriangulation(vertices,faces,testp);
%      % Plot results
%      h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
%      set(h,'FaceColor','black','FaceAlpha',1/3,'EdgeColor','none');
%      hold on;
%      plot3(testp(:,1),testp(:,2),testp(:,3),'b.');
%      plot3(testp(in==1,1),testp(in==1,2),testp(in==1,3),'ro');
%
% See also: intetrahedron, tsearchn, inpolygon
%
% Author: Johannes Korsawe, heavily based on voxelise from Adam A.
% E-mail: johannes.korsawe@volkswagen.de
% Release: 1.3
% Release date: 25/09/2013
% check number of inputs
if nargin<3,
    fprintf('??? Error using ==> intriangulation\nThree input matrices are needed.\n');in=[];return;
end
if nargin==3,
    heavytest = 0;
end
% check size of inputs
if size(vertices,2)~=3 || size(faces,2)~=3 || size(testp,2)~=3,
    fprintf('??? Error using ==> intriagulation\nAll input matrices must have three columns.\n');in=[];return;
end
ipmax = max(faces(:));zerofound = ~isempty(find(faces(:)==0, 1));
if ipmax>size(vertices,1) || zerofound,
    fprintf('??? Error using ==> intriangulation\nThe triangulation data is defect. use trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3)) for test of deficiency.\n');return;
end
% loop for heavytest
inreturn = zeros(size(testp,1),1);VER = vertices;TESTP = testp;
for n = 1:heavytest+1,
    % Randomize
    if n>1,
        v=rand(1,3);D=rotmatrix(v/norm(v),rand*180/pi);vertices=VER*D;testp = TESTP*D;
    else,
        vertices=VER;
    end
    
    % Preprocessing data
    meshXYZ = zeros(size(faces,1),3,3);
    for loop_ind = 1:3,
        meshXYZ(:,:,loop_ind) = vertices(faces(:,loop_ind),:);
    end
    % Basic idea (ingenious from FeX-submission voxelise):
    % If point is inside, it will cross the triangulation an uneven number of times in each direction (x, -x, y, -y, z, -z).
    
    % The function VOXELISEinternal is about 98% identical to its version inside voxelise.m.
    % This includes the elaborate comments. Thanks to Adam A!
    
    % z-direction:
    % intialization of results and correction list
    [in,cl] = VOXELISEinternal(testp(:,1),testp(:,2),testp(:,3),meshXYZ);
    
    % x-direction:
    % has only to be done for those points, that were not determinable in the first step --> cl
    [in2,cl2] = VOXELISEinternal(testp(cl,2),testp(cl,3),testp(cl,1),meshXYZ(:,[2,3,1],:));
    % Use results of x-direction that determined "inside"
    in(cl(in2==1)) = 1;
    % remaining indices with unclear result
    cl = cl(cl2);
    
    % y-direction:
    % has only to be done for those points, that were not determinable in the first and second step --> cl
    [in3,cl3] = VOXELISEinternal(testp(cl,3),testp(cl,1),testp(cl,2),meshXYZ(:,[3,1,2],:));
    
    % Use results of y-direction that determined "inside"
    in(cl(in3==1)) = 1;
    % remaining indices with unclear result
    cl = cl(cl3);
    
    % mark those indices, where all three tests have failed
    in(cl) = -1;
    
    if n==1,
        inreturn = in;  % Starting guess
    else,
        % if ALWAYS inside, use as inside!
%        I = find(inreturn ~= in);
%        inreturn(I(in(I)==0)) = 0;
        
        % if AT LEAST ONCE inside, use as inside!
        I = find(inreturn ~= in);
        inreturn(I(in(I)==1)) = 1;
        
    end
end
in = inreturn;

%==========================================================================
function [OUTPUT,correctionLIST] = VOXELISEinternal(testx,testy,testz,meshXYZ)
% Prepare logical array to hold the logical data:
OUTPUT = false(size(testx,1),1);
%Identify the min and max x,y coordinates of the mesh:
meshZmin = min(min(meshXYZ(:,3,:)));meshZmax = max(max(meshXYZ(:,3,:)));
%Identify the min and max x,y,z coordinates of each facet:
meshXYZmin = min(meshXYZ,[],3);meshXYZmax = max(meshXYZ,[],3);
%======================================================
% TURN OFF DIVIDE-BY-ZERO WARNINGS
%======================================================
%This prevents the Y1predicted, Y2predicted, Y3predicted and YRpredicted
%calculations creating divide-by-zero warnings.  Suppressing these warnings
%doesn't affect the code, because only the sign of the result is important.
%That is, 'Inf' and '-Inf' results are ok.
%The warning will be returned to its original state at the end of the code.
warningrestorestate = warning('query', 'MATLAB:divideByZero');
%warning off MATLAB:divideByZero
%======================================================
% START COMPUTATION
%======================================================
correctionLIST = [];   %Prepare to record all rays that fail the voxelisation.  This array is built on-the-fly, but since
%it ought to be relatively small should not incur too much of a speed penalty.
% Loop through each testpoint.
% The testpoint-array will be tested by passing rays in the z-direction through
% each x,y coordinate of the testpoints, and finding the locations where the rays cross the mesh.
facetCROSSLIST = zeros(1,1e3);  % uses countindex: nf
nm = size(meshXYZmin,1);
for loop = 1:length(OUTPUT),
    
    nf = 0;
%    % - 1a - Find which mesh facets could possibly be crossed by the ray:
%    possibleCROSSLISTy = find( meshXYZmin(:,2)<=testy(loop) & meshXYZmax(:,2)>=testy(loop) );
%    % - 1b - Find which mesh facets could possibly be crossed by the ray:
%    possibleCROSSLIST = possibleCROSSLISTy( meshXYZmin(possibleCROSSLISTy,1)<=testx(loop) & meshXYZmax(possibleCROSSLISTy,1)>=testx(loop) );
    % Do - 1a - and - 1b - faster
    possibleCROSSLISTy = find((testy(loop)-meshXYZmin(:,2)).*(meshXYZmax(:,2)-testy(loop))>0);
    possibleCROSSLISTx = (testx(loop)-meshXYZmin(possibleCROSSLISTy,1)).*(meshXYZmax(possibleCROSSLISTy,1)-testx(loop))>0;
    possibleCROSSLIST = possibleCROSSLISTy(possibleCROSSLISTx);
    
    if isempty(possibleCROSSLIST)==0  %Only continue the analysis if some nearby facets were actually identified
        
        % - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:
        
        % GENERAL METHOD:
        % 1. Take each edge of the facet in turn.
        % 2. Find the position of the opposing vertex to that edge.
        % 3. Find the position of the ray relative to that edge.
        % 4. Check if ray is on the same side of the edge as the opposing vertex.
        % 5. If this is true for all three edges, then the ray definitely passes through the facet.
        %
        % NOTES:
        % 1. If the ray crosses exactly on an edge, this is counted as crossing the facet.
        % 2. If a ray crosses exactly on a vertex, this is also taken into account.
        
        for loopCHECKFACET = possibleCROSSLIST'
            
            %Check if ray crosses the facet.  This method is much (>>10 times) faster than using the built-in function 'inpolygon'.
            %Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.  If so, let testVn=1
            
            Y1predicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,1))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
            YRpredicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-testx(loop))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
            
            if (Y1predicted > meshXYZ(loopCHECKFACET,2,1) && YRpredicted > testy(loop)) || (Y1predicted < meshXYZ(loopCHECKFACET,2,1) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-testx(loop)) == 0
%                testV1 = 1;   %The ray is on the same side of the 2-3 edge as the 1st vertex.
            else
%                testV1 = 0;   %The ray is on the opposite side of the 2-3 edge to the 1st vertex.
                % As the check is for ALL three checks to be true, we can continue here, if only one check fails
                continue;
            end %if
            
            Y2predicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,2))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
            YRpredicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-testx(loop))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
            if (Y2predicted > meshXYZ(loopCHECKFACET,2,2) && YRpredicted > testy(loop)) || (Y2predicted < meshXYZ(loopCHECKFACET,2,2) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-testx(loop)) == 0
%                testV2 = 1;   %The ray is on the same side of the 3-1 edge as the 2nd vertex.
            else
%                testV2 = 0;   %The ray is on the opposite side of the 3-1 edge to the 2nd vertex.
                % As the check is for ALL three checks to be true, we can continue here, if only one check fails
                continue;
            end %if            
            
            Y3predicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,3))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
            YRpredicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-testx(loop))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
            if (Y3predicted > meshXYZ(loopCHECKFACET,2,3) && YRpredicted > testy(loop)) || (Y3predicted < meshXYZ(loopCHECKFACET,2,3) && YRpredicted < testy(loop)) || (meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-testx(loop)) == 0
%                testV3 = 1;   %The ray is on the same side of the 1-2 edge as the 3rd vertex.
            else
%                testV3 = 0;   %The ray is on the opposite side of the 1-2 edge to the 3rd vertex.
                % As the check is for ALL three checks to be true, we can continue here, if only one check fails
                continue;
            end %if
    
            nf=nf+1;facetCROSSLIST(nf)=loopCHECKFACET;
            
        end %for
        % Use only values ~=0
        facetCROSSLIST = facetCROSSLIST(1:nf);
        
        % - 3 - Find the z coordinate of the locations where the ray crosses each facet:
        gridCOzCROSS = zeros(1,nf);
        for loopFINDZ = facetCROSSLIST
            
            % METHOD:
            % 1. Define the equation describing the plane of the facet.  For a
            % more detailed outline of the maths, see:
            % http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
            %    Ax + By + Cz + D = 0
            %    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
            %           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
            %           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
            %           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
            % 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.
            
            planecoA = meshXYZ(loopFINDZ,2,1)*(meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,3,3)) + meshXYZ(loopFINDZ,2,2)*(meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,3,1)) + meshXYZ(loopFINDZ,2,3)*(meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,3,2));
            planecoB = meshXYZ(loopFINDZ,3,1)*(meshXYZ(loopFINDZ,1,2)-meshXYZ(loopFINDZ,1,3)) + meshXYZ(loopFINDZ,3,2)*(meshXYZ(loopFINDZ,1,3)-meshXYZ(loopFINDZ,1,1)) + meshXYZ(loopFINDZ,3,3)*(meshXYZ(loopFINDZ,1,1)-meshXYZ(loopFINDZ,1,2));
            planecoC = meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)-meshXYZ(loopFINDZ,2,3)) + meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)-meshXYZ(loopFINDZ,2,1)) + meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)-meshXYZ(loopFINDZ,2,2));
            planecoD = - meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,2)) - meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,3)) - meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,1));
            
            if abs(planecoC) < 1e-14
                planecoC=0;
            end
            
            gridCOzCROSS(facetCROSSLIST==loopFINDZ) = (- planecoD - planecoA*testx(loop) - planecoB*testy(loop)) / planecoC;
            
        end %for
        if isempty(gridCOzCROSS),continue;end
        
        %Remove values of gridCOzCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
        gridCOzCROSS = gridCOzCROSS( gridCOzCROSS>=meshZmin-1e-12 & gridCOzCROSS<=meshZmax+1e-12 );
        %Round gridCOzCROSS to remove any rounding errors, and take only the unique values:
        gridCOzCROSS = round(gridCOzCROSS*1e10)/1e10;
        
        % Replacement of the call to unique (gridCOzCROSS = unique(gridCOzCROSS);) by the following line:
        tmp = sort(gridCOzCROSS);I=[0,tmp(2:end)-tmp(1:end-1)]~=0;gridCOzCROSS = [tmp(1),tmp(I)];
        
        % - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:
        
        if rem(numel(gridCOzCROSS),2)==0  % Only rays which cross an even number of facets are voxelised
            
            for loopASSIGN = 1:(numel(gridCOzCROSS)/2)
                voxelsINSIDE = (testz(loop)>gridCOzCROSS(2*loopASSIGN-1) & testz(loop)<gridCOzCROSS(2*loopASSIGN));
                OUTPUT(loop) = voxelsINSIDE;
                if voxelsINSIDE,break;end
            end %for
            
        elseif numel(gridCOzCROSS)~=0    % Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
            correctionLIST = [ correctionLIST; loop ];
        end %if
        
    end %if
    
end %for
%======================================================
% RESTORE DIVIDE-BY-ZERO WARNINGS TO THE ORIGINAL STATE
%======================================================
warning(warningrestorestate)
% J.Korsawe: A correction is not possible as the testpoints need not to be
%            ordered in any way.
%            voxelise contains a correction algorithm which is appended here
%            without changes in syntax.
return
%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = size(correctionLIST,1);
if countCORRECTIONLIST>0
    
  %If necessary, add a one-pixel border around the x and y edges of the
  %array.  This prevents an error if the code tries to interpolate a ray at
  %the edge of the x,y grid.
  if min(correctionLIST(:,1))==1 || max(correctionLIST(:,1))==numel(gridCOx) || min(correctionLIST(:,2))==1 || max(correctionLIST(:,2))==numel(gridCOy)
    gridOUTPUT     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),gridOUTPUT,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
    correctionLIST = correctionLIST + 1;
  end
  
  for loopC = 1:countCORRECTIONLIST
    voxelsforcorrection = squeeze( sum( [ gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)-1,:) ,...
                                          gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2),:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)+1,:) ,...
                                          gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)-1,:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)+1,:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)-1,:) ,...
                                          gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2),:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)+1,:) ,...
                                         ] ) );
    voxelsforcorrection = (voxelsforcorrection>=4);
    gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),voxelsforcorrection) = 1;
  end %for
  %Remove the one-pixel border surrounding the array, if this was added
  %previously.
  if size(gridOUTPUT,1)>numel(gridCOx) || size(gridOUTPUT,2)>numel(gridCOy)
    gridOUTPUT = gridOUTPUT(2:end-1,2:end-1,:);
  end
  
end %if
%disp([' Ray tracing result: ',num2str(countCORRECTIONLIST),' rays (',num2str(countCORRECTIONLIST/(voxcountX*voxcountY)*100,'%5.1f'),'% of all rays) exactly crossed a facet edge and had to be computed by interpolation.'])
end %function
%==========================================================================
function D = rotmatrix(v,deg)
% calculate the rotation matrix about v by deg degrees
deg=deg/180*pi;
if deg~=0,
    v=v/norm(v);
    v1=v(1);v2=v(2);v3=v(3);ca=cos(deg);sa=sin(deg);
    D=[ca+v1*v1*(1-ca),v1*v2*(1-ca)-v3*sa,v1*v3*(1-ca)+v2*sa;
        v2*v1*(1-ca)+v3*sa,ca+v2*v2*(1-ca),v2*v3*(1-ca)-v1*sa;
        v3*v1*(1-ca)-v2*sa,v3*v2*(1-ca)+v1*sa,ca+v3*v3*(1-ca)];
else,
    D=eye(3,3);
end
end

end
