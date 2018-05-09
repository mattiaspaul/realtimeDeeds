function [ b ] = regularizedRegression( A, x, lambda )
%% REGULARIZEDREGRESSION: 
%% INPUT:
%   A: model (#elements x #eigenvectors)
%   x: displacement field   (#elements)
%   lambda: regularization parameter
%% OUTPUT:
%   b: new weight for the model (#eigenvectors)

[U,S,V]=svd(A,'econ');
W=diag(S)./(diag(S).^2+lambda);
b=V*diag(W)*U'*x;
end

