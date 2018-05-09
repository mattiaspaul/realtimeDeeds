% show results:
clear all;
load result/result_im028.mat;

flow_rgb=flowToColor(cat(3,ux,vx));
figure; subplot(1,3,1); imshow(flow_rgb); title('colour-coding of deformation field','FontSize',16);
subplot(1,3,2); imagesc(magn); title('magnitude of deformation field','FontSize',16);
axis image; axis off;
subplot(1,3,3); imagesc(im1-im2w,[-200,200]); title('difference image after registration','FontSize',16);
 axis image; axis off;
