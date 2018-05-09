function [flow1,flow1i]=consistentInverse(flow1,flow1i,grid1,it)
%consistent inverse according to MICCAI 2013 paper

if nargin<3
    grid1=1;
end
if nargin<4
    it=10;
end

flow1=flow1./grid1;
flow1i=flow1i./grid1;

for i=1:it

u2=flow1(:,:,1); v2=flow1(:,:,2);
u2i=flow1i(:,:,1); v2i=flow1i(:,:,2);

flow1(:,:,1)=0.5*(u2-imWarp(u2,v2,single(u2i)));
flow1(:,:,2)=0.5*(v2-imWarp(u2,v2,single(v2i)));

flow1i(:,:,1)=0.5*(u2i-imWarp(u2i,v2i,single(u2)));
flow1i(:,:,2)=0.5*(v2i-imWarp(u2i,v2i,single(v2)));

end

flow1=flow1.*grid1;
flow1i=flow1i.*grid1;