function im2=imresizeG(im1,factor)
[m,n,p]=size(im1);
m2=floor(m/factor);
n2=floor(n/factor);
im2=zeros(m2,n2,p,'single');
for i=1:factor
    for j=1:factor
        im2=im2+im1(i:factor:m2*factor+i-1,j:factor:n2*factor+j-1,:);
    end
end
im2=im2./(factor^2);