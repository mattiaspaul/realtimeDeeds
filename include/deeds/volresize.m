function vol2=volresize(vol1,newsize,method)
%-- VOLRESIZE: creates the meshgrid in vol1 range
%-- vol1: the value range
%-- newsize: the number of values
%
%--- if the size of the vol1 is same as the newsize: vol2=vol1
if sum(newsize==size(vol1))==3
    vol2=vol1;
else
    
a=single(linspace(1,size(vol1,1),newsize(1)));
b=single(linspace(1,size(vol1,2),newsize(2)));
c=single(linspace(1,size(vol1,3),newsize(3)));

%--- zi: 1s and 2s for xy-dim
[xi,yi,zi]=meshgrid(b,a,c);

%--- if the input arguments are less than 3, then do trilinear
%interpolation (no method selection)
if nargin<3
    vol2=trilinearSingle(single(vol1),single(xi),single(yi),single(zi));

else
    
    vol2=interp3(single(vol1),xi,yi,zi,method);
end

%--- change NaNs to 0
vol2(isnan(vol2))=0;
vol2=single(vol2);

end