%% deeds symmetric 2d image registration
clear all;

%before running for first time compile mex-files
mex trilinearSingle.cpp
mex regMarg2dincr.cpp
%load example_images.mat
%200 images as a matrix
%load slice17.mat

%% initalise flow fields to zero
flow0=zeros(2,2,2,'single');
flow0i=zeros(2,2,2,'single');

% specify registration parameters (using four scale levels)
grid=[7,6,5,4]; %grid-point spacing
maxdisp=[16,15,8,5]; %maximum displacement in each dimension/direction
quant=[4,3,2,1]; %quantisation step in pixel
alpha=10; %regularisation parameter (greater values=smoother field)

% this step has to be done for 199 times.
for n=1:199
    %im1 is the reference image
    im1(:,:) = X(1,:,:);
    im2(:,:) = X(n+1,:,:);
    
    for i=1:4
    grid1=grid(i); quant1=quant(i);  maxdisp1=maxdisp(i);
    flow1=deedsReg1(im1,im2,flow0,grid1,quant1,alpha,maxdisp1); %forward transform
    flow1i=deedsReg1(im2,im1,flow0i,grid1,quant1,alpha,maxdisp1); %backward transform
    [flow1,flow1i]=consistentInverse(flow1,flow1i,grid1,10); %consistent Inverse

    flow0=flow1; flow0i=flow1i;

    end

    flow2=volresize(flow1,[size(im1),2]);
    ux=flow2(:,:,1); vx=flow2(:,:,2);
    magn=sqrt(ux.^2+vx.^2);
    im2w=imWarp(ux,vx,im2);
    
    str = sprintf('/share/data_tiffy3/ha/MotionEstimation/Results/2DRegistration/sl%3d/displacement_sl35_%03d.mat',pt,n);
    save(str, 'flow0','flow0i','flow1','flow1i','flow2','ux','vx','im1','im2','im2w','magn');
end

%this is just for the display of the result
% flow_rgb=flowToColor(cat(3,ux,vx));
% figure; subplot(1,3,1); imshow(flow_rgb); title('colour-coding of deformation field','FontSize',16);
% subplot(1,3,2); imagesc(magn); title('magnitude of deformation field','FontSize',16);
% axis image; axis off;
% subplot(1,3,3); imagesc(im1-im2w,[-200,200]); title('difference image after registration','FontSize',16);
%  axis image; axis off;

