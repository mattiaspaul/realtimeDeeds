function [ mu, U, D ] = buildMotionModel( inputFields, variability, mu )
%% BUILDMOTIONMODEL generates motion model
%% INPUT:
%  inputFields: dense motion fields (#Elements x #TrainingSamples)
%  variability: model variability (scalar between 0 and 1)
%  (optional) mu: mean motion field 
%% OUTPUT:
%  mu: mean motion field
%  U: eigenvectors, i.e. motion model (#Elements x #Eigenvectors)
%  D: squared eigenvalues (#Eigenvectors)

if nargin < 3
  mu=mean(inputFields,2);
end
N=size(inputFields,2);

inputFields=bsxfun(@minus,inputFields,mu);
[U,D]=eig(((1/sqrt(N-1))*inputFields)'*((1/sqrt(N-1))*inputFields),'vector');

% -- sort eigenvectors/values
[D,idx] = sort(D,'descend');      
U = U(:,idx);

% -- compute the number of eigenvalues needed
numOfEV = N-1;
if variability < 1
s = 0;
numOfEV = 0;
    while s<variability
        numOfEV = numOfEV+1;
        s = sum(D(1:numOfEV))/sum(D);
    end
end
%-- compute Eigenvectors for cov(dataMat')
U = normc(inputFields*U(:,1:numOfEV));
D = D(1:numOfEV);

fprintf('%.0f Eigenvectors are needed to retain %.2f %% variability..\n',numOfEV,variability*100);

end

