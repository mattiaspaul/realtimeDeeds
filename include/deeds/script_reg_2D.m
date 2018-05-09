%% deeds symmetric 2d image registration
%clear all;
function [UX_2D, WX_2D] = script_reg_2D (len,ref_nii,mov_nii)
%before running for first time compile mex-files
mex ../deedsRegistration2D/trilinearSingle.cpp
mex ../deedsRegistration2D/regMarg2dincr.cpp
%load example_images.mat

% displacement field
UX_2D = zeros(len-1,224,224);
WX_2D = UX_2D;

%% initalise flow fields to zero
flow0=zeros(2,2,2,'single');
flow0i=zeros(2,2,2,'single');

% specify registration parameters (using four scale levels)
grid=[7,6,5,4]; %grid-point spacing
maxdisp=[32,21,12,5]; %maximum displacement in each dimension/direction
quant=[4,3,2,1]; %quantisation step in pixel
alpha=10; %regularisation parameter (greater values=smoother field)

% this step has to be done for length(timestep) times.
for l=1:len
    %im1 is the reference image
    im1(:,:) = single(ref_nii.img(l,:,:));
    im2(:,:) = single(mov_nii.img(l,:,:));
    
    for i=1:4
    grid1=grid(i); quant1=quant(i);  maxdisp1=maxdisp(i);
    flow1=deedsReg1(im1,im2,flow0,grid1,quant1,alpha,maxdisp1); %forward transform
    flow1i=deedsReg1(im2,im1,flow0i,grid1,quant1,alpha,maxdisp1); %backward transform
    [flow1,flow1i]=consistentInverse(flow1,flow1i,grid1,10); %consistent Inverse

    flow0=flow1; flow0i=flow1i;

    end

    flow2=volresize(flow1,[size(im1),2]);
    ux=flow2(:,:,1); vx=flow2(:,:,2);
    magn=sqrt(ux.^2+vx.^2);     %---magnitude
    im2w=imWarp(ux,vx,im2);
    
    UX_2D(l,:,:) = vx;
    WX_2D(l,:,:) = ux;
    %M(:,:,p) = magn;
    %IM(:,:,p+1) = im2;
end

% str = sprintf('C:/Users/InYoung/Documents/1. Uni-Luebeck/0. WS1516/Praktikum/data/evaluation/sl%02d_reg/result_sl%02d_%03d.mat',pt,pt,snum);
% save(str, 'UX','VX');
fprintf('slice registration is finished..\n');
%save(str,'flow1','flow1i','flow2','ux','vx','im1','im2','magn');
end
%this is just for the display of the result
% flow_rgb=flowToColor(cat(3,ux,vx));
% figure; subplot(1,3,1); imshow(flow_rgb); title('colour-coding of deformation field','FontSize',16);
% subplot(1,3,2); imagesc(magn); title('magnitude of deformation field','FontSize',16);
% axis image; axis off;
% subplot(1,3,3); imagesc(im1-im2w,[-200,200]); title('difference image after registration','FontSize',16);
%  axis image; axis off;

