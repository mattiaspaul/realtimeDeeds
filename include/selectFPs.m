function [ fpoints, fpointsImage ] = selectFPs( image, paramsHarris)
%% SELECTFPS: extracts feature points from the image
%% INPUT:
%   image: (w x h)
%   paramsHarris: [sigma, thresh, radius for NMS, display]
%% OUTPUT:
%   fpoints: coordinate of selected feature points (#points x 2)
%   fpointsImage: binary image with feature points set to 1 (w x h)

%-- initialize
fpoints = [];
fpointsImage=zeros(size(image));
len =0;  

%-- feature point selection
% harris detector
[cim, r, c] = harris(image, paramsHarris(1), paramsHarris(2), paramsHarris(3), paramsHarris(4));
fpoints=[fpoints; r c];                    
len = len+size(r,1);

cell=mat2cell(fpoints,size(fpoints,1),[1 1]);
idx=sub2ind(size(fpointsImage),cell{:});
fpointsImage(idx)=1;

fprintf('total %.0f feature points are found... \n',len);
end

