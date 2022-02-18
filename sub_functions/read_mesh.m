
function [coil_mesh,target_mesh,shielded_mesh]=read_mesh(input)


%Read the input mesh
if strcmp(input.sf_source_file,'none')

%Read the coil mesh surface
coil_mesh = stlread_local(strcat(input.geometry_source_path,'\',input.coil_mesh_file));
coil_mesh=create_unique_noded_mesh(coil_mesh);

%Change to format convenction used here
coil_mesh.vertices=coil_mesh.vertices'; 
coil_mesh.faces=coil_mesh.faces'; 


else
loaded_file=load(strcat(cd,'\',input.sf_source_file));
coil_mesh=loaded_file.coil_mesh;

end


%Read the target mesh surface
if ~strcmp(input.target_mesh_file,'none')
target_mesh = stlread_local(strcat(input.geometry_source_path,'\',input.target_mesh_file));
target_mesh=create_unique_noded_mesh(target_mesh);
else
target_mesh=[];
end

%Read the shielded mesh surface
if ~strcmp(input.secondary_target_mesh_file,'none')
shielded_mesh = stlread_local(strcat(input.geometry_source_path,'\',input.secondary_target_mesh_file));
shielded_mesh=create_unique_noded_mesh(shielded_mesh);
else
shielded_mesh=[];
end



% %Permute the nodes of the current carrying mesh as well as the target
% %points
% %Distribue the index of the nodes randomly
% permutet_indices=randperm(size(coil_mesh.vertices,2));
% coil_mesh.faces=[permutet_indices(coil_mesh.faces(:,1))' permutet_indices(coil_mesh.faces(:,2))' permutet_indices(coil_mesh.faces(:,3))'];
% non_permutet_vertices=coil_mesh.vertices;
% for vert_ind=1:size(coil_mesh.vertices,2)
%     coil_mesh.vertices(:,vert_ind)=non_permutet_vertices(:,permutet_indices==vert_ind);
% end





function unique_noded_mesh=create_unique_noded_mesh(non_unique_mesh)
faces=non_unique_mesh.faces;
verts=non_unique_mesh.vertices;
%Remove double nodes
[unique_verts,~,unique_inds] = unique(verts, 'rows','stable');
unique_assignments= [[1:numel(unique_inds)]' unique_inds];
for tri_ind=1:size(faces,1)
faces(tri_ind,1)=unique_assignments(faces(tri_ind,1),2);
faces(tri_ind,2)=unique_assignments(faces(tri_ind,2),2);
faces(tri_ind,3)=unique_assignments(faces(tri_ind,3),2);
end
unique_noded_mesh.faces=faces;
unique_noded_mesh.vertices=unique_verts;
end



function varargout = stlread_local(file)
% STLREAD imports geometry from an STL file into MATLAB.
%    FV = STLREAD(FILENAME) imports triangular faces from the ASCII or binary
%    STL file idicated by FILENAME, and returns the patch struct FV, with fields
%    'faces' and 'vertices'.
%
%    [F,V] = STLREAD(FILENAME) returns the faces F and vertices V separately.
%
%    [F,V,N] = STLREAD(FILENAME) also returns the face normal vectors.
%
%    The faces and vertices are arranged in the format used by the PATCH plot
%    object.
% Copyright 2011 The MathWorks, Inc.
    if ~exist(file,'file')
        error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
               ', be sure to specify the full path to the file.'], file);
    end
    
    fid = fopen(file,'r');    
    if ~isempty(ferror(fid))
        error(lasterror); %#ok
    end
    
    M = fread(fid,inf,'uint8=>uint8');
    fclose(fid);
    
    [f,v,n] = stlbinary(M);
    %if( isbinary(M) ) % This may not be a reliable test
    %    [f,v,n] = stlbinary(M);
    %else
    %    [f,v,n] = stlascii(M);
    %end
    
    varargout = cell(1,nargout);
    switch nargout        
        case 2
            varargout{1} = f;
            varargout{2} = v;
        case 3
            varargout{1} = f;
            varargout{2} = v;
            varargout{3} = n;
        otherwise
            varargout{1} = struct('faces',f,'vertices',v);
    end

function [F,V,N] = stlbinary(M)
    F = [];
    V = [];
    N = [];
    
    if length(M) < 84
        error('MATLAB:stlread:incorrectFormat', ...
              'Incomplete header information in binary STL file.');
    end
    
    % Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
    % that follow.
    numFaces = typecast(M(81:84),'uint32');
    %numFaces = double(numFaces);
    if numFaces == 0
        warning('MATLAB:stlread:nodata','No data in STL file.');
        return
    end
    
    T = M(85:end);
    F = NaN(numFaces,3);
    V = NaN(3*numFaces,3);
    N = NaN(numFaces,3);
    
    numRead = 0;
    while numRead < numFaces
        % Each facet is 50 bytes
        %  - Three single precision values specifying the face normal vector
        %  - Three single precision values specifying the first vertex (XYZ)
        %  - Three single precision values specifying the second vertex (XYZ)
        %  - Three single precision values specifying the third vertex (XYZ)
        %  - Two unused bytes
        i1    = 50 * numRead + 1;
        i2    = i1 + 50 - 1;
        facet = T(i1:i2)';
        
        n  = typecast(facet(1:12),'single');
        v1 = typecast(facet(13:24),'single');
        v2 = typecast(facet(25:36),'single');
        v3 = typecast(facet(37:48),'single');
        
        n = double(n);
        v = double([v1; v2; v3]);
        
        % Figure out where to fit these new vertices, and the face, in the
        % larger F and V collections.        
        fInd  = numRead + 1;        
        vInd1 = 3 * (fInd - 1) + 1;
        vInd2 = vInd1 + 3 - 1;
        
        V(vInd1:vInd2,:) = v;
        F(fInd,:)        = vInd1:vInd2;
        N(fInd,:)        = n;
        
        numRead = numRead + 1;
    end
    
end
function [F,V,N] = stlascii(M)
    warning('MATLAB:stlread:ascii','ASCII STL files currently not supported.');
    F = [];
    V = [];
    N = [];
end
% TODO: Change the testing criteria! Some binary STL files still begin with
% 'solid'.
function tf = isbinary(A)
% ISBINARY uses the first line of an STL file to identify its format.
    if isempty(A) || length(A) < 5
        error('MATLAB:stlread:incorrectFormat', ...
              'File does not appear to be an ASCII or binary STL file.');
    end    
    if strcmpi('solid',char(A(1:5)'))
        tf = false; % ASCII
    else
        tf = true;  % Binary
    end
end
end


end
