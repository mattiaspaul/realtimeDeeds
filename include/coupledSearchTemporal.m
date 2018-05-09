function [xDisp,yDisp,bestind]=coupledSearchTemporal(distSSC,U,lambdas,search,mu,D,maxComp,beta,oldField)
%% COUPLEDSEARCHTEMPORAL: performs coupled convex optimization with temporal term
%% INPUT:
%   distSSC: matrix of dissimilarities (#keypoint x #displacements)
%   PCA model of learning phase (including mean and serparate for x-,y-motion)
%   lambdas: weighting of regularisation terms (e.g. [0.03,0.1,0.3,1])
%   maxComp: number of PCA components to be use (e.g. [3,5,7,9])
%   search: half-width of search window
%% OUTPUT:
%   xDisp, yDisp: regularised (sparse) displacements
%   bestind: indicies with smallest SSD cost

%search-grid and initial displacement
[xs,ys] = meshgrid(-search:search,-search:search);
[~,bestind] = min(distSSC,[],2);
yDisp = ys(bestind(:)); xDisp = xs(bestind(:));
bestindOld = bestind;

field = zeros(numel(xDisp)*2,1);
xDisp_old = xDisp; yDisp_old = yDisp;
intermedFields = field;
%iterate over coupled optimisation
for i=1:length(lambdas)
        field(1:2:end) = xDisp;
        field(2:2:end) = yDisp;
        weights = regularizedRegression(U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i)))),field-mu,.1);        
        approxField = round((U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i))))*weights)+mu);
        xDisp_reg = approxField(1:2:end);
        yDisp_reg = approxField(2:2:end);
        
        %coupling term (squared Euclidean distance to model displacement)
        r = ((repmat(xs(:)',length(xDisp),1)-repmat(xDisp_reg(:),1,numel(xs))).^2+...
            (repmat(ys(:)',length(yDisp),1)-repmat(yDisp_reg(:),1,numel(xs))).^2)./(search.^2);
        
        %coupling term (squared Euclidean distance to model old field)
        r2 = ((repmat(xs(:)',length(xDisp),1)-repmat(oldField(1:2:end),1,numel(xs))).^2+...
            (repmat(ys(:)',length(yDisp),1)-repmat(oldField(2:2:end),1,numel(xs))).^2)./(search.^2);
        
        [~,bestind] = min(distSSC+lambdas(i).*r+beta(i).*r2,[],2);        
        
        yDisp = ys(bestind(:)); xDisp = xs(bestind(:));

%         fprintf('changed: %f %f\n',sum(sqrt((xDisp_reg-xDisp).^2+...
%             (yDisp_reg-yDisp).^2)),sum(sqrt((xDisp_old-xDisp).^2+...
%             (yDisp_old-yDisp).^2))./size(distSSC,1));
        bestindOld=bestind;
        
        field(1:2:end) = xDisp;
        field(2:2:end) = yDisp;
        weights = regularizedRegression(U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i)))),field-mu,.1);        
        approxField = round((U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i))))*weights)+mu);
        xDisp_reg = approxField(1:2:end);
        yDisp_reg = approxField(2:2:end); 
                
        indices = sub2ind(size(distSSC),(1:size(distSSC,1))',bestind);       
        xDisp_old = xDisp; yDisp_old = yDisp;
end
