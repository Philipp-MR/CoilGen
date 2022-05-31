function  mesh = mesh_parameterization_iterative(mesh_in)
%beginn with the parameterization
%Schmidt  rms@dgp.toronto.edu" based on desbrun et al (2002), "Intrinsic Parameterizations of {Surface} Meshes"


%find all the usefull mesh properties
mesh = struct('v', [], 'n', [], 'u', [], 'f', [], 'e', [], 'bounds', []);
mesh.v=mesh_in.vertices';
mesh.u=[];
mesh.f=mesh_in.faces';
mesh.bounds=[ min(mesh.v) ; max(mesh.v) ; 0.5*(min(mesh.v) + max(mesh.v)) ];
mesh.version=1;
mesh.vidx=[1:size(mesh.v,1)]';
mesh.fidx=[1:size(mesh.f,1)]';
mesh.fn =mesh_in.fn;
for iii = 1:size(mesh.f,1)
fvv = mesh.v( mesh.f(iii,:), : );
ee1 = fvv(2,:) - fvv(1,:);
ee2 = fvv(3,:) - fvv(1,:);
M=cross(ee1, ee2) ;
mag = sum(M .* M, 2);
z = find(mag < eps);
mag(z) = 1;
Mlen = sqrt(mag);
mesh.fn(iii,:) = M ./ repmat(Mlen, 1, size(M,2));
end
e = zeros(size(mesh.f,1)*3*2,3);
for iii = 1:size(mesh.f,1)
for jjj = 0:2
t1 = [ mesh.f(iii,jjj+1), mesh.f(iii,mod(jjj+1,3)+1), 1 ];
t2 = [ mesh.f(iii,mod(jjj+1,3)+1), mesh.f(iii,jjj+1), 1 ];
e(((iii-1)*6)+(jjj*2)+1,:) = t1;
e(((iii-1)*6)+(jjj*2)+2,:) = t2;
end
end
mesh.e = sparse( e(:,1), e(:,2), e(:,3), size(mesh.v,1), size(mesh.v,1) );
[iiii,jjjj] = find(mesh.e==1);
mesh.isboundaryv = zeros(size(mesh.v,1),1);
mesh.isboundaryv(iiii) = 1;
be = sortrows( sort([iiii,jjjj],2) );
be = unique(be, 'rows');
mesh.isboundaryf = ( mesh.isboundaryv(mesh.f(:,1)) + mesh.isboundaryv(mesh.f(:,2)) + mesh.isboundaryv(mesh.f(:,3)) );
mesh.isboundaryf = mesh.isboundaryf > 0;
loops = [];
loopk = 1;
bloop = [];
while numel(be) > 0
bloop = [];
a1 = be(1,1); a2 = be(1,2);
be(1,:) = [];
bloop = cat(1,bloop,a2);
while size(be,1) > 0
nextrow = find( be(:,1) == a2 & be(:,2) ~= a1 );
if nextrow
  b2 = be(nextrow,1); b3 = be(nextrow,2);
else
  nextrow = find( be(:,2) == a2 & be(:,1) ~= a1 ); 
  b3 = be(nextrow,1); b2 = be(nextrow,2);
end
if isempty(nextrow)
   loops{loopk} = bloop;
   loopk = loopk+1;
   break;
else
   be(nextrow,:) = [];
   bloop = cat(1,bloop,b3);
   a1 = b2; a2 = b3;
end
end
end
if ~isempty(bloop) 
loops{loopk} = bloop;
loopk = loopk+1;
end
for kkkk = 1:numel(loops)
loop_sort = loops{kkkk};
prev_idx = [3,1,2];
loop1 = loop_sort(1);   loop2 = loop_sort(2);
[ffi,fj] = find(mesh.f == loop1);
for iii = 1:numel(ffi)
jp = prev_idx( fj(iii) );
if mesh.f(ffi(iii), jp) == loop2
nL = size(loop_sort,1);
loop_sort = loop_sort(nL:-1:1);
end
end   
loops{kkkk} = loop_sort;
end
if ~ isempty(loops)
for kkkk = 1:numel(loops)
loopsize(kkkk) = numel(loops{kkkk});
end
[~,idx] = sort(loopsize,'descend');
for kkkk = 1:numel(idx)
mesh.loops{kkkk} = loops{idx(kkkk)};
end
else
mesh.loops = [];
end
mesh.te = mesh.e;
mesh.te( mesh.e~=0 ) = 0;
for ti = 1:size(mesh.f,1)
for kkkk = 0:2
vv1 = mesh.f(ti,kkkk+1);  vv2 = mesh.f(ti, mod((kkkk+1),3)+1 );
if mesh.te( vv1, vv2 ) ~= 0
ttmp = vv1; vv1 = vv2; vv2 = ttmp;
end
mesh.te( vv1, vv2 ) = ti;
end
end
mesh.valence = zeros(size(mesh.v,1),1);
for vi = 1:size(mesh.v,1)
[~,jj] = find( mesh.e(vi,:) );    
mesh.valence(vi) = size(jj,2);
end
mesh.unique_vert_inds=mesh_in.unique_vert_inds;
mesh.n = vertexNormal(triangulation(mesh_in.faces',mesh_in.vertices'));
mesh.fn = faceNormal(triangulation(mesh_in.faces',mesh_in.vertices'));
iboundary = mesh.vidx(mesh.isboundaryv~=0);
dists = vmag2(vadd(mesh.v(iboundary,:), -mesh.v( iboundary(1),:)));
[~,maxi] = max(dists);
ifixed = [ iboundary(1),iboundary(maxi)];
fixedUV = [ iboundary(1), 0,0;  iboundary(maxi),1,0];
N = numel(mesh.vidx);
W = cotanWeights(mesh);
W = -W;
W(1:N+1:N*N) = -sum(W,2);
Ld = [W,sparse(N,N); sparse(N,N), W];
rhs = zeros(2*N,1);
A = sparse([],[],[],2*N,2*N,4*numel(iboundary));
for li = 1:numel(mesh.loops)
loop = mesh.loops{li};
for ii = 1:numel(loop)
jx = loop(ii);
jy = jx + N;
kx = loop( mod(ii,numel(loop)) + 1 );
ky = kx + N;
A(jx,ky) = A(jx,ky) + 1;
A(kx,jy) = A(kx,jy) - 1;   
end
end
A = (A + A');
Lc = Ld - A;
LcCons = Lc;
ifixed = [ifixed,ifixed+N];
LcCons(ifixed,:) = 0;
LcCons(sub2ind(size(Lc),ifixed, ifixed)) = 1;
rhs(fixedUV(:,1)) = fixedUV(:,2);
rhs(fixedUV(:,1)+N) = fixedUV(:,3);
rhsadd = zeros(size(rhs));
for k = 1:numel(ifixed)
ci = ifixed(k);
col = LcCons(:,ci);
col(ci) = 0;
rhsadd = rhsadd + rhs(ci)*col;
end
LcCons(:,ifixed) = 0;
LcCons(sub2ind(size(Lc),ifixed, ifixed)) = 1;
rhs = rhs-rhsadd;
mesh.uv = LcCons \ rhs;
mesh.uv = [mesh.uv(1:N), mesh.uv(N+1:2*N)];

mesh.vertices=mesh.v';
mesh.faces=mesh.f';
mesh.n=mesh.n';
mesh.uv=mesh.uv';
mesh.boundary=mesh.loops;



%%% Sub functions of mesh_parameterization_func
function [ Mout ] = vmag2( M )
Mout = sum( M .* M, 2 );
end
function [ Mout ] = vadd( M, v )
Mout = M + repmat(v, size(M,1), 1);
end
function [ cotAB ] = vcot( A, B )
tmp = dot(A,B);
cotAB = tmp / sqrt( dot(A,A)*dot(B,B) - tmp*tmp );
end
function [ vVertices ] = oneringv( mesh, nVertex )
vVertices = find(mesh.e(nVertex,:)~=0)';
end
function [ A ] = faceArea( mesh, faces )
if ~ exist('faces', 'var') || numel(faces) == 0
    faces = mesh.fidx;
end
n = numel(faces);
A = zeros(n,1);
for i = 1:n
    f = mesh.f(i,:);
    fv = mesh.v(f,:);
    A(i) = triarea(fv(1,:), fv(2,:), fv(3,:));
end
end
function [ A ] = triarea( p1, p2, p3 )
    u = p2 - p1;
    v = p3 - p1;
    A = 0.5 * sqrt( dot(u,u)*dot(v,v) - dot(u,v)^2 );
end
function [ W ] = cotanWeights( mesh, vertices, authalic, areaWeighted )
if ~ exist('vertices', 'var') || numel(vertices) == 0
    vertices = mesh.vidx;
end
if ~exist('authalic', 'var')
    authalic = 0;
end
if ~exist('areaWeighted', 'var')
    areaWeighted = 0;
end
n = numel(vertices);
W = mesh.e(vertices,:);
W(W~=0)=-1;
if areaWeighted
faceAreas = faceArea(mesh);
else
faceAreas = ones(numel(mesh.fidx),1);
end
for i = 1:n
qi = vertices(i);
ov = oneringv(mesh,qi);
for j = 1:numel(ov)
qj = ov(j);
faces = [mesh.te(qi,qj), mesh.te(qj,qi)];
faces = faces(faces~=0);
verts = zeros(numel(faces),1);
vertfaces = zeros(numel(faces),1);
for kk = 1:numel(faces)
f = mesh.f(faces(kk),:);
verts(kk) = f(f~=qi & f~=qj);
vertfaces(kk) = faces(kk);
end
sumAB = 0;
if authalic
for kk = 1:numel(verts)
qo = verts(kk);
v1 = mesh.v(qi,:)-mesh.v(qj,:);
v2 = mesh.v(qo,:)-mesh.v(qj,:);
sumAB = sumAB + vcot(v1,v2);
end               
sumAB = sumAB / vmag2( mesh.v(qi,:) - mesh.v(qj,:) );
else
for kk = 1:numel(verts)
qo = verts(kk);
v1 = mesh.v(qi,:)-mesh.v(qo,:);
v2 = mesh.v(qj,:)-mesh.v(qo,:);
sumAB = sumAB + vcot(v1,v2) / faceAreas(vertfaces(kk));
end
end
W(qi,qj) = sumAB;  
end
end
end
function [ vTris, vVerts ] = onering( mesh, nVert, mode )
vVerts = find(mesh.e(nVert,:)~=0)';
vTris = cat(1, full(mesh.te(nVert,vVerts))', full(mesh.te(vVerts,nVert)));
vTris = unique(vTris);
vTris(find(vTris==0)) = [];
if ~exist('mode','var')
return;
end
if strcmp(mode,'ccw') && ~isempty(vVerts)
if mesh.isboundaryv(nVert)
isb = mesh.isboundaryv(vVerts);
swapi = find(isb~=0);
tmp = vVerts(1); vVerts(1) = vVerts(swapi(1));  vVerts(swapi(1)) = tmp;
end
curv = vVerts(1);
vSorted = [curv];  
rest = vVerts(2:end);
tnbrs = [mesh.te(nVert,curv), mesh.te(curv,nVert)];
vnbrs = [0,0];
for j = 1:2
if tnbrs(j) ~= 0 
    vnbrs(j) = tripick2(mesh.f(tnbrs(j),:),nVert,curv);
end
end
prev = curv;
usev = 1;  if vnbrs(usev) == 0; usev = 2; end;
curv = vnbrs(usev);
tSorted = [tnbrs(usev)];
while ~isempty(rest)
vSorted = [vSorted,curv];
rest = rest(rest~=curv);
tnbrs = [mesh.te(nVert,curv), mesh.te(curv,nVert)];
if tnbrs(1) == 0 | tnbrs(2) == 0
    break;
end  
vnbrs = [0,0];
for j = 1:2
    if tnbrs(j) ~= 0 
        vnbrs(j) = tripick2(mesh.f(tnbrs(j),:),nVert,curv);
    end
end        
if vnbrs(1) == prev
    prev = curv;
    curv = vnbrs(2);
    tSorted = [tSorted, tnbrs(2)];
elseif vnbrs(2) == prev
    prev = curv;
    curv = vnbrs(1);
    tSorted = [tSorted, tnbrs(1)];
end
end
vTris = tSorted;
vVerts = vSorted;
end
end
function [ k ] = tripick2( face, i, j )
fi = find(face~=i & face~=j);
k = face(fi);
end
function [ Mout ] = vmag( M )
Mout = sqrt( sum( M .* M, 2 ) );
end
function [ nc ] = ncross( v1, v2 )
    nc = normalize( cross( v1, v2 ) );
end
end