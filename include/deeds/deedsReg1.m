function u1=deedsReg1(im1,im2,u0,step1,quant,alpha,maxdisp)
%--DEEDSREG1 
%-- im1: reference image
%-- im2: target image
%-- u0: displacement matrix (2:2:2) -> (32:32:2)... 
%-- step1: grid [7,6,5,4]
%-- quant: quant [4,3,2,1]
%-- alpha: 10 regularisation parameter
%-- maxdisp: [32,21,12,5] maxmum displacement in each direction/dimension

%-- changes the precision to single. 
%-- alpha1 = 7/40, 6/30, 5/20, 4/10
im1=single(im1); im2=single(im2);
alpha1=step1/(quant*alpha);

%displacement space definition
[xs,ys]=meshgrid(-maxdisp:quant:maxdisp,-maxdisp:quant:maxdisp);

%--- volresize changes the grid size of u0 to fit sizegrid
sizegrid=floor(size(im1)./step1);
u0=single(volresize(u0,[sizegrid,2]));

%upsampling of flow to image resolution for warping of moving image
uv=volresize(u0,[size(im1),2]);
im2w=imWarp(uv(:,:,1),uv(:,:,2),im2);

%calculate data-term: sum of absolute differences, downsample to grid-level
dataD=zeros([sizegrid,numel(xs)],'single');
for i=1:numel(xs)
    dataD(:,:,i)=imresizeG(alpha1*abs(im1-imshift2(im2w,xs(i),ys(i))),step1);
end

%infer regularisation using MST-based belief propagation
marg=regMarg2dincr(dataD,im1,step1,u0(:,:,1),u0(:,:,2),quant);

%select minimum and add to previous field
[~,indx]=min(marg,[],3);
u1=xs(indx); u1(:,:,2)=ys(indx);
u1=u1+u0;
