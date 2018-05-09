%% deeds symmetric 2d image registration
%clear all;

%before running for first time compile mex-files
mex ../deedsRegistration2D/trilinearSingle.cpp
mex ../deedsRegistration2D/regMarg2dincr.cpp
%load example_images.mat
%load sl04_slice16.mat
%load slice17.mat
% load(loadfile);

dataMat = zeros(224*224*2,len);
% save image
%IM = zeros(224,224,len);
%IM(:,:,1) = im1;

%% initalise flow fields to zero


% specify registration parameters (using four scale levels)
grid=[7,6,5,4]; %grid-point spacing
maxdisp=[16,15,8,5]; %maximum displacement in each dimension/direction
quant=[4,3,2,1]; %quantisation step in pixel
alpha=10; %regularisation parameter (greater values=smoother field)

% this step has to be done for length(timestep) times.
im1(:,:) = X(1,:,:);

for k=1:len
    %im1 is the reference image
    im2(:,:) = X(k+1,:,:);
    flow0=zeros(2,2,2,'single');
    flow0i=zeros(2,2,2,'single');
    for i=1:4
    grid1=grid(i); quant1=quant(i);  maxdisp1=maxdisp(i);
    flow1=deedsReg1ssc(im1,im2,flow0,grid1,quant1,alpha,maxdisp1); %forward transform
    flow1i=deedsReg1ssc(im2,im1,flow0i,grid1,quant1,alpha,maxdisp1); %backward transform
    [flow1,flow1i]=consistentInverse(flow1,flow1i,grid1,10); %consistent Inverse

    flow0=flow1; flow0i=flow1i;

    end

    flow2=volresize(flow1,[size(im1),2]);
    ux=flow2(:,:,1); vx=flow2(:,:,2);
    magn=sqrt(ux.^2+vx.^2);     %---magnitude
    im2w=imWarp(ux,vx,im2);
    
    dataMat(:,k) = [vx(:);ux(:)]; % x,y 
end
% 
% str = sprintf('C:/Users/InYoung/Documents/1. Uni-Luebeck/0. WS1516/Praktikum/data/evaluation/sl%02d_reg/result_sl%02d_%03d.mat',pt,pt,snum);
 save('dataMat.mat', 'dataMat');
 clc;
fprintf('slice registration is finished..\n');
%save(str,'flow1','flow1i','flow2','ux','vx','im1','im2','magn');
    
%this is just for the display of the result
% flow_rgb=flowToColor(cat(3,ux,vx));
% figure; subplot(1,3,1); imshow(flow_rgb); title('colour-coding of deformation field','FontSize',16);
% subplot(1,3,2); imagesc(magn); title('magnitude of deformation field','FontSize',16);
% axis image; axis off;
% subplot(1,3,3); imagesc(im1-im2w,[-200,200]); title('difference image after registration','FontSize',16);
%  axis image; axis off;

