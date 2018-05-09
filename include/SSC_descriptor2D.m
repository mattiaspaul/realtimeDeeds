function ssc=SSC_descriptor2D(I,sigma,delta)
% Calculation of SSC (self-similarity context)
%
% If you use this implementation please cite:
% M.P. Heinrich et al.: "Towards Realtime Multimodal Fusion for 
% Image-Guided Interventions Using Self-similarities"
% MICCAI (2013) LNCS Springer
%
% M.P. Heinrich et al.: "MIND: Modality Independent Neighbourhood
% Descriptor for Multi-Modal Deformable Registration"
% Medical Image Analysis (2012)
%
% Contact: heinrich(at)imi(dot)uni-luebeck(dot)de
%
% I: input volume (2D)
% sigma: Gaussian weighting for patches 
% delta: Distance between patch centres
%
% ssc: output descriptor (3D)

I=single(I);
if nargin<2
    sigma=0.8;
end
if nargin<3
    delta=2;
end

% Filter for efficient patch SSD calculation
filt=fspecial('gaussian',[ceil(sigma*3/2)*2+1,ceil(sigma*3/2)*2+1],sigma);

%displacements of patches
dx=[+0,+1,+0,-1,+1,+0].*delta; % 1, 2, 1, 0, 2, 1
dy=[-1,+0,+1,+0,+0,+1].*delta; % 0, 1, 2, 1, 1, 2
sx=[-1,+0,+1,+0,-1,+0].*delta; % 0, 1, 2, 1, 0, 1
sy=[+0,-1,+0,+1,+0,-1].*delta; % 1, 0, 1, 2, 1, 0

% Self-similarity Distances
ssc=zeros([size(I),numel(dx)],'single');

% Calculating Gaussian weighted patch SSD using convolution
for i=1:numel(dx)
    ssc(:,:,i)=imfilter((imshift(I,sx(i),sy(i))-imshift(I,dx(i),dy(i))).^2,filt);
end

% Remove minimal distance to scale descriptor to max=1
ssc=ssc-repmat(min(ssc,[],3),1,1,numel(dx));

% Variance measure (standard absolute deviation)
V=mean(ssc,3);
val1=[0.001*(mean(V(:))),1000*mean(V(:))];
V=(min(max(V,min(val1)),max(val1)));

% descriptor calculation according
ssc=exp(-ssc./repmat(V,1,1,numel(dx)));


function im1shift=imshift(im1,x,y)

[m,n,o]=size(im1);

im1shift=im1;
x1s=max(1,x+1);
x2s=min(n,n+x);

y1s=max(1,y+1);
y2s=min(m,m+y);

x1=max(1,-x+1);
x2=min(n,n-x);

y1=max(1,-y+1);
y2=min(m,m-y);

im1shift(y1:y2,x1:x2,:)=im1(y1s:y2s,x1s:x2s,:);

