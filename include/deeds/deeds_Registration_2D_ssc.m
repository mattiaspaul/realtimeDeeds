%% deeds symmetric 2d image registration
close all

for f=2:200;
    im1=frames(:,:,1); im2=frames(:,:,f);
%%
% initalise flow fields to zero
flow0=zeros(2,2,2,'single');
flow0i=zeros(2,2,2,'single');

% specify registration parameters (using four scale levels)
grid=[6,5,4,3]; %grid-point spacing
quant=[3,2,1,1]; %quantisation step in pixel
maxdisp=[7,6,5,4].*quant;%[21,14,7,3.5]; %maximum displacement in each dimension/direction

alpha=1; %regularisation parameter (greater values=smoother field)
figure; subplot(1,3,1); imagesc(im1-im2,[-100,100]); axis image; axis off;
title('difference image before registration','FontSize',16);

tic;
for i=1:4
grid1=grid(i); quant1=quant(i);  maxdisp1=maxdisp(i);

[flow1_,flow1]=deedsReg1ssc(im1,im2,flow0,grid1,quant1,alpha,maxdisp1); %forward transform
[flow1i_,flow1i]=deedsReg1ssc(im2,im1,flow0i,grid1,quant1,alpha,maxdisp1); %backward transform
[flow1,flow1i]=consistentInverse(flow1,flow1i,grid1,10); %consistent Inverse

flow0=flow1; flow0i=flow1i;
toc;
end

flow2=volresize(flow1,[size(im1),2]);
ux=flow2(:,:,1); vx=flow2(:,:,2);
magn=sqrt(ux.^2+vx.^2);
im2w=imWarp(ux,vx,single(im2));
% grey2jet(normalize
% subplot(1,3,2); imagesc((magn).*256); title('magnitude of deformation field','FontSize',16);
% axis image; axis off;
% subplot(1,3,3); imagesc(single(im1)-im2w,[-100,100]); title('difference image after registration','FontSize',16);
% colormap gray; axis image; axis off; drawnow;

end

%[aae,tre]=angularError(gtx,gty,ux,vx,21)
%show2map(gtx,gty,ux,vx);
