function [ B ] = imWarp( flowHor, flowVer, Bin, method )
if nargin<4
    method='linear';
end
%This function warps B towards A
[m,n,o]=size(Bin);

[x y] = meshgrid(1:size(Bin,2),1:size(Bin,1));
% xu is the grid of size nxm
xu=x+flowHor;
xu(xu<1)=1;
xu(xu>n)=n;
xu(isnan(xu))=0;

yv=y+flowVer;
yv(yv<1)=1;
yv(yv>m)=m;
yv(isnan(yv))=0;

%-- interp2: Bin(:,:,i) the value at the grid points
%-- xu, yu: the points where the interpolation values are needed
%-- output is the interpolated values at (xu,yv)
for i=1:size(Bin,3)
    B(:,:,i) = interp2(Bin(:,:,i), xu,yv, method);   
end
B(isnan(B)) = Bin(isnan(B));
end
