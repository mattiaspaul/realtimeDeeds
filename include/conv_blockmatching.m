function [distSSC] = conv_blockmatching(grefSum2,testImageSSC,fpoints,gRefBlock,gfilt,blockSize,searchSize,sscDim)
%% CONV_BLOCKMATCHING: computes Block-Matching on GPU using convolution (A²-2AB+B²)
%% INPUT: 
%   grefSum2: A² (searchSize x searchSize x Nfpoints)
%   im: target Image
%   grefBlock: reference image patches (blockSize x blockSize x sscDim*Nfpoints) 
%   gfilt: filter for B² computation
%   blockSize: size of image block
%   searchSize: size of search region
%   sscDim: number of SSC descriptor (for 2D = 6)
%% OUTPUT:
%   distSSC: (searchSize x searchSize x Nfpoints)

%% extract target image patch
Nfpoints = size(grefSum2,3);
div = (blockSize.*blockSize.*sscDim);
pad = (blockSize-1)/2 + (searchSize-1)/2;

tarSearch = zeros(searchSize+blockSize-1,searchSize+blockSize-1,sscDim,Nfpoints,'single');    
testImageSSC = padarray(testImageSSC,[pad, pad, 0]);

for j=1:Nfpoints
   tarSearch(:,:,:,j)=squeeze(testImageSSC(fpoints(j,1)+pad-pad:fpoints(j,1)+pad+pad,...
       fpoints(j,2)+pad-pad:fpoints(j,2)+pad+pad,:));   
end

% compute AB and B²
gTarSearch = gpuArray(reshape(tarSearch,searchSize+blockSize-1,searchSize+blockSize-1,[]));
gy_filt = vl_nnconv(gTarSearch,reshape(gRefBlock,blockSize,blockSize,sscDim,Nfpoints),[]);
gtargetSum2= vl_nnconv(gTarSearch.^2,gfilt,[]);
gtargetSum2= reshape(gtargetSum2,searchSize,searchSize,Nfpoints);

% compute SSD for all feature points
y_filt = reshape(gy_filt,searchSize,searchSize,Nfpoints)./div;
gdistSSC = gtargetSum2./div-2.*y_filt+grefSum2;   
distSSC = gather(gdistSSC);

end