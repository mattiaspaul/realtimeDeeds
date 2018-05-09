function [imshift2b]=imshift2(im1,x,y,pad)

% tic;
% im1shift=imWarp(ones(size(im1)).*x,ones(size(im1)).*y,im1);
% toc;


imshift2b=zeros(size(im1));
x1=x-floor(x);
x2=1-x1; %ceil(x)-x;
y1=y-floor(y);
y2=1-y1; %ceil(y)-y;

imshift2b=imshift2b+x2*y2.*imshift(im1,floor(x),floor(y));
imshift2b=imshift2b+x1*y2.*imshift(im1,ceil(x),floor(y));
imshift2b=imshift2b+x1*y1.*imshift(im1,ceil(x),ceil(y));
imshift2b=imshift2b+x2*y1.*imshift(im1,floor(x),ceil(y));


