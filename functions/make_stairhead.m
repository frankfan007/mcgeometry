function stackhead=make_stairhead(head)

% uses image erosion to create a multi-tissue layered head volume 

% input:
%   head: volume, dimension (nx,ny,nz)

% output:
%   stackhead: volume that was "eroded", dimension (nx,ny,nz)

% author: Stefan Carp, <stefan.carp@mgh.harvard.edu>
% this function is part of the mcgeometry toolbox,
%(https://github.com/wumelissa/mc_geometry)

%%

savehead=head;
head=savehead;
head(head>0)=1;

stackhead=head;

for I=1:20
    [xx,yy,zz] = ndgrid(-1:1);
    nhood = sqrt(xx.^2 + yy.^2 + zz.^2) <= 1.0;
    head = imerode(head,nhood);
    stackhead=stackhead+head; 
end

nVoxX=size(head,1);
nVoxY=size(head,2);
nVoxZ=size(head,3);

stackhead = reshape(stackhead,[1 nVoxX*nVoxY*nVoxZ]);
