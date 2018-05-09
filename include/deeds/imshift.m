function im1shift=imshift(im1,x,y,pad)
if nargin<4
    pad=0;
end
[m,n,o]=size(im1);


im1shift=im1;%ones(size(im1)).*pad;

x1s=max(1,x+1);
x2s=min(n,n+x);

y1s=max(1,y+1);
y2s=min(m,m+y);

x1=max(1,-x+1);
x2=min(n,n-x);

y1=max(1,-y+1);
y2=min(m,m-y);

% length1=x2-x1+1;
% length2=x2s-x1s+1;
% length3=y2-y1+1;
% length4=y2s-y1s+1;


im1shift(y1:y2,x1:x2,:)=im1(y1s:y2s,x1s:x2s,:);