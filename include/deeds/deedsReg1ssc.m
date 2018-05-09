function [u1,us]=deedsReg1ssc(im1,im2,u0,step1,quant,alpha,maxdisp)

im1=single(im1); im2=single(im2);
alpha1=step1/(quant*alpha);

%displacement space definition
[xs,ys]=meshgrid(-maxdisp:quant:maxdisp,-maxdisp:quant:maxdisp);
sizegrid=floor(size(im1)./step1);
u0=single(volresize(u0,[sizegrid,2]));

%upsampling of flow to image resolution for warping of moving image
uv=volresize(u0,[size(im1),2]);
im2w=imWarp(uv(:,:,1),uv(:,:,2),im2);

ssc_im1=SSC_descriptor2D(im1,1.2,3);
ssc_im2w=SSC_descriptor2D(im2w,1.2,3);

%calculate data-term: sum of absolute differences, downsample to grid-level
dataD=zeros([sizegrid,numel(xs)],'single');
for i=1:numel(xs)
    dataD(:,:,i)=imresizeG(alpha1*sum(abs(ssc_im1-imshift2(ssc_im2w,xs(i),ys(i))),3),step1);
end

%infer regularisation using MST-based belief propagation
marg=regMarg2dincr(dataD,im1,step1,u0(:,:,1),u0(:,:,2),quant);

%select minimum and add to previous field
[~,indx]=min(marg,[],3);
u1=xs(indx); u1(:,:,2)=ys(indx);
u1=u1+u0;

if nargout>1
xcoord={1:size(u1,1),1:size(u1,2)};
us=csaps(xcoord,double(u1(:,:,1)),[0.1,0.1],xcoord);
us(:,:,2)=csaps(xcoord,double(u1(:,:,2)),[0.1,0.1],xcoord);
end
