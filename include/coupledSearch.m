function [xDisp,yDisp,bestind]=coupledSearch(distSSC,U,lambdas,search,mu,D,maxComp)
%% COUPLEDSEARCH: performs coupled convex optimization
%% INPUT: 
% distSSC: matrix of dissimilarities (#keypoint x #displacements)
% U: PCA model of learning phase (including mean and serparate for x-,y-motion)
% lambdas: weighting of regularisation terms (e.g. [0.03,0.1,0.3,1])
% search: half-width of search window
% mu: mean of training data
% D: squared eigenvalues
% maxComp: number of PCA components to be use (e.g. [3,5,7,9])
%% OUTPUT:
% xDisp, yDisp: regularised (sparse) displacements
% bestind: indicies of minimum values in distSSC

%% search-grid and initial displacement
[xs,ys]=meshgrid(-search:search,-search:search);
[~,bestind]=min(distSSC,[],2);
yDisp=ys(bestind(:)); xDisp=xs(bestind(:));
bestindOld=bestind;

field=zeros(numel(xDisp)*2,1);
xDisp_old=xDisp; yDisp_old=yDisp;
intermedFields=field;
%% iterate over coupled optimisation
for i=1:length(lambdas)
        field(1:2:end)=xDisp;
        field(2:2:end)=yDisp;
        
        % compute the weights based on the given sparse field ans model
        weights=regularizedRegression(U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i)))),field-mu,1);        
        stddevLimits=100.*ones(size(weights));
        stddevLimitsIndicator=abs(weights)>stddevLimits;
        weights(stddevLimitsIndicator)=sign(weights(stddevLimitsIndicator)).*stddevLimits(stddevLimitsIndicator);
        
        % reconstruct the sparse displacement field using the computed
        % weight vector and the model.
        approxField=round((U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i))))*weights)+mu);
        xDisp_reg=approxField(1:2:end);
        yDisp_reg=approxField(2:2:end);
        
        %% coupling term (squared Euclidean distance to model displacement)
        r=((repmat(xs(:)',length(xDisp),1)-repmat(xDisp_reg(:),1,numel(xs))).^2 +...
        (repmat(ys(:)',length(yDisp),1)-repmat(yDisp_reg(:),1,numel(xs))).^2)./(search.^2);
        [~,bestind2]=min(r,[],2);
        [~,bestind]=min(distSSC+lambdas(i).*r,[],2);        
        yDisp=ys(bestind(:)); xDisp=xs(bestind(:));
        intermedFields(1:2:end,i)=xDisp;
        intermedFields(2:2:end,i)=yDisp;

%         fprintf('changed: %f %f\n',sum(sqrt((xDisp_reg-xDisp).^2+(yDisp_reg-yDisp).^2)),...
%             sum(sqrt((xDisp_old-xDisp).^2+(yDisp_old-yDisp).^2))./size(distSSC,1));
        bestindOld=bestind;
        
        field(1:2:end)=xDisp;
        field(2:2:end)=yDisp;
        weights=regularizedRegression(U(:,1:maxComp(i))*diag(sqrt(D(1:maxComp(i)))),field-mu,.1);        

        xDisp_old=xDisp; yDisp_old=yDisp;
end


