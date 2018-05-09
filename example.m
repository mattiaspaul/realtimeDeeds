%% example code
% tested with MATLAB R2017a

clear all; clc; close all;

% add paths 
addpath('./extern');
addpath('./extern/NIfTI_20140122');
addpath('./extern/harris');
addpath('./include');

dataset_names={{'volunteerA','9'},{'volunteerA','26'},{'volunteerB','8'},...
    {'volunteerB','26'},{'volunteerC','8'},{'volunteerC','26'},{'volunteerD','7'},...
    {'volunteerD','27'}};

%% set dataset paths 
patString = dataset_names{1}{1};
sliceString = dataset_names{1}{2};

staticImageFilename = ['./Data/',patString,'/bhs/sl_ss',sliceString,'_bh.nii.gz'];
dynamicImagePath = ['./Data/',patString,'/dyn/images/'];
fieldsPath = ['./Data/',patString,'/dyn/deeds_bh/'];
lmksPath = ['./Data/',patString,'/dyn/lmks/'];

%% load static image
staticImage = load_nii(staticImageFilename);
spacing = staticImage.hdr.dime.pixdim(2);
mask = logical(ones(size(staticImage.img))-(regiongrowing(staticImage.img,1,1,3000)+...
    regiongrowing(staticImage.img,1,size(staticImage.img,2),3000)));
staticImageLM = fliplr(csvread([lmksPath,'sl_ss',sliceString,'_bh_lmks.txt']))+1;

%% find feature points
% harrisParams = [sigma, thresh, radius for NMS, display];
harrisParams = [2, 1000000, 8, 1];

% fpoints: (numOfFeaturePoints x 2)
% fpointsImage: binary image with size (imh x imw)
[fpoints, fpointsImage] = selectFPs(staticImage.img.*mask,harrisParams);
fpointsImageMask = fpointsImage(mask);
fpointsFieldMask = zeros(numel(fpointsImageMask)*2,1);
fpointsFieldMask(1:2:end) = fpointsImageMask;
fpointsFieldMask(2:2:end) = fpointsImageMask;
fpointsFieldMask = logical(fpointsFieldMask);

%% load images, fields and build model
numOfTrainingFields = 20; numOfsamples = 40; 
numOfTestdata = 20; variability = 0.98;

testImages = zeros([size(staticImage.img),numOfTestdata]);

for frame=1:numOfTestdata
    imageFilename = [dynamicImagePath,'sl_ss',sliceString,'_tt',num2str(frame+numOfTrainingFields),'.nii.gz'];
    image = load_nii(imageFilename);
    testImages(:,:,frame) = image.img;
end

fields = loadDeedsFieldsMasked(1:numOfsamples,sliceString,fieldsPath,mask,1);
[muStandard, UStandard, DStandard] = buildMotionModel( fields(:,1:numOfTrainingFields), variability);

%% load test images and perform block matching
BMBlockSize = 5;    BMSearchWindow = 15;
Nfpoints = numel(fpoints(:,1));    
blockSize = 2*BMBlockSize+1;   searchSize = BMSearchWindow*2+1;      
pad = BMSearchWindow+BMBlockSize;
[xs,ys] = meshgrid(-BMSearchWindow:BMSearchWindow,-BMSearchWindow:BMSearchWindow);

% compute SSC descriptor for reference image
SSCsigma = 1.2;  SSCdelta = 3;
refImageSSC = SSC_descriptor2D(staticImage.img,SSCsigma,SSCdelta);
refImageSSC = padarray(refImageSSC,[pad,pad,0]);
sscDim = size(refImageSSC,3);  

%---- prepare for the Block-Matching using GPU
BMIndex = [];
refBlock = zeros(blockSize,blockSize,sscDim,Nfpoints,'single');
filt = ones(blockSize,blockSize,sscDim,Nfpoints,'single');
gfilt = gpuArray(filt);

for j=1:Nfpoints
   refBlock(:,:,:,j) = refImageSSC(fpoints(j,1)+pad-BMBlockSize:fpoints(j,1)+pad+BMBlockSize,...
                       fpoints(j,2)+pad-BMBlockSize:fpoints(j,2)+pad+BMBlockSize,:);
end

gRefBlock = gpuArray(reshape(refBlock,blockSize,blockSize,[]));
gRefSum2 = vl_nnconv(gRefBlock.^2,gfilt,[]);
grefSum2 = repmat(gRefSum2,[searchSize,searchSize,1])./(blockSize*blockSize*sscDim);
grefSum2 = reshape(grefSum2,searchSize,searchSize,Nfpoints);
fprintf('variables for GPU block-matching are prepared...\n');    

% sparse model
muSparse=muStandard(fpointsFieldMask);
USparse=UStandard(fpointsFieldMask,:);

%% compute block-matching, coupled convex optimization with and without temporal term
TREidx = [];    TREBM = [];     TREBMCoupled = [];  TREBMCoupledTemporal = [];
tstart = tic;
for frame=1:numOfTestdata
    filenameLMs = [lmksPath,'sl_ss',sliceString,'_tt',num2str(frame+numOfTrainingFields),'_lmks.txt'];    
    oldField = fields(:,1);
        
    %% Block-Matching: results in distance Map (distSSC)
%     disp(['BM Frame: ',num2str(frame)]);

    testImageSSC = SSC_descriptor2D(squeeze(testImages(:,:,frame)),SSCsigma,SSCdelta);  
    distSSC = conv_blockmatching(grefSum2,testImageSSC,fpoints,...
              gRefBlock,gfilt,blockSize,searchSize,sscDim);
    distSSC = reshape(permute(distSSC,[3,1,2]),Nfpoints,[]);
    [~,bestind] = min(distSSC,[],2);
    
    % compute sparse displacement field
    xDisp = xs(bestind(:)); yDisp = ys(bestind(:));
    
    %% Block-Matching field reconstruction
    bmField = zeros(size(fields(:,1),1),1);
    bmDisp = zeros(size(xDisp,1)*2,1);    
    bmDisp(1:2:end) = xDisp;     bmDisp(2:2:end) = yDisp;
    bmField(fpointsFieldMask) = bmDisp;
    
    % reconstruction of dense field from block-matching result
    approxWeights = USparse\(bmDisp-muSparse);
    stddevLimits = 3*sqrt(DStandard(1:size(UStandard,2)));
    stddevLimitsIndicator = abs(approxWeights)>stddevLimits;
    approxWeights(stddevLimitsIndicator) = sign(approxWeights(stddevLimitsIndicator)).*...
                                           stddevLimits(stddevLimitsIndicator);
    approxFieldBM = UStandard*approxWeights+muStandard;
    
    %% coupled convex optimization
    lambdas=logspace(-1.5,0,4);
    maxComp = [1,1,1,1,1,1,1,1,1,1,1,1,1,1].*size(UStandard,2);
    coupledDisp = zeros(size(xDisp,1)*2,1);
    
    [coupledDisp(1:2:end),coupledDisp(2:2:end),bestindCS] = coupledSearch(distSSC,...
                           USparse,lambdas,BMSearchWindow,muSparse,DStandard,maxComp);
    
    % field reconstruction from coupled convex optimization result
    approxWeights = regularizedRegression(USparse*diag(sqrt(DStandard(1:size(UStandard,2)))),...
                                            coupledDisp-muSparse,1);        
    stddevLimits = 3.*ones(size(approxWeights));
    stddevLimitsIndicator = abs(approxWeights)>stddevLimits;
    approxWeights(stddevLimitsIndicator) = sign(approxWeights(stddevLimitsIndicator)).*...
                                           stddevLimits(stddevLimitsIndicator);      
    approxFieldBMCoupled = (UStandard*diag(sqrt(DStandard(1:size(UStandard,2))))*...
                            approxWeights)+muStandard;
                        
    %% coupled convex optimization (temporal)

    beta = [0.5,0.25,0,0,0,0];
    coupledDispTemporal = zeros(size(xDisp,1)*2,1);
    [coupledDispTemporal(1:2:end),coupledDispTemporal(2:2:end),bestindCS] = ...
        coupledSearchTemporal(distSSC,USparse,lambdas,BMSearchWindow,muSparse,...
        DStandard,maxComp,beta,oldField(fpointsFieldMask));
    
    % field reconstruction from coupled convex optimization (temporal) result
    approxWeights = regularizedRegression(USparse*diag(sqrt(DStandard(1:size(UStandard,2)))),...
                                          coupledDispTemporal-muSparse,1);        
    stddevLimits = 3.*ones(size(approxWeights));
    stddevLimitsIndicator = abs(approxWeights)>stddevLimits;
    approxWeights(stddevLimitsIndicator) = sign(approxWeights(stddevLimitsIndicator)).*...
                                           stddevLimits(stddevLimitsIndicator);            
    approxFieldBMCoupledTemporal = (UStandard*diag(sqrt(DStandard(1:size(UStandard,2))))*...
                                    approxWeights)+muStandard;

    
    %% compute TRE
    if exist(filenameLMs, 'file') == 2
        currentLMs = fliplr(csvread(filenameLMs))+1;
        TREidx = [TREidx frame+numOfTestdata];

        % load deeds dense field
        fullField = loadDeedsField(frame,sliceString,fieldsPath,1);
        xFullField = fullField(:,:,1);
        yFullField = fullField(:,:,2); 

        % TRE Block-Matching
        xFullField(mask) = approxFieldBM(1:2:end);
        yFullField(mask) = approxFieldBM(2:2:end);
        fullField(:,:,1) = xFullField;
        fullField(:,:,2) = yFullField;
        griddedInterpolantX = griddedInterpolant(fullField(:,:,1));
        griddedInterpolantY = griddedInterpolant(fullField(:,:,2));
        displ(:,1) = griddedInterpolantX(staticImageLM(:,2),staticImageLM(:,1));
        displ(:,2) = griddedInterpolantY(staticImageLM(:,2),staticImageLM(:,1));
        warpedLMs = displ  +staticImageLM;
        tre = mean(sqrt(sum((warpedLMs-currentLMs).^2,2))*spacing);
        TREBM = [TREBM tre];     

        % TRE coupled convex optimization
        xFullField(mask) = approxFieldBMCoupled(1:2:end);
        yFullField(mask) = approxFieldBMCoupled(2:2:end);
        fullField(:,:,1) = xFullField;
        fullField(:,:,2) = yFullField;        
        griddedInterpolantX = griddedInterpolant(fullField(:,:,1));
        griddedInterpolantY = griddedInterpolant(fullField(:,:,2));
        displ(:,1) = griddedInterpolantX(staticImageLM(:,2),staticImageLM(:,1));
        displ(:,2) = griddedInterpolantY(staticImageLM(:,2),staticImageLM(:,1));
        warpedLMs = displ + staticImageLM;
        tre = mean(sqrt(sum((warpedLMs-currentLMs).^2,2))*spacing);
        TREBMCoupled = [TREBMCoupled tre];    

        % TRE coupled convex optimization (temporal)
        xFullField(mask) = approxFieldBMCoupledTemporal(1:2:end);
        yFullField(mask) = approxFieldBMCoupledTemporal(2:2:end);
        fullField(:,:,1) = xFullField;
        fullField(:,:,2) = yFullField;        
        griddedInterpolantX = griddedInterpolant(fullField(:,:,1));
        griddedInterpolantY = griddedInterpolant(fullField(:,:,2));
        displ(:,1) = griddedInterpolantX(staticImageLM(:,2),staticImageLM(:,1));
        displ(:,2) = griddedInterpolantY(staticImageLM(:,2),staticImageLM(:,1));
        warpedLMs = displ + staticImageLM;
        tre = mean(sqrt(sum((warpedLMs-currentLMs).^2,2))*spacing);
        TREBMCoupledTemporal = [TREBMCoupledTemporal tre];  
    end

end
tend = toc(tstart);
fprintf('\nmean TRE: \n');
fprintf('block-matching \t\t = %.2f mm\n',mean(TREBM));
fprintf('coupled convex \t\t = %.2f mm\n',mean(TREBMCoupled));
fprintf('coupled convex temporal  = %.2f mm\n',mean(TREBMCoupledTemporal));
fprintf('computation time \t = %.2f s (%.2f s/frame)\n',tend,tend/numOfTestdata);
