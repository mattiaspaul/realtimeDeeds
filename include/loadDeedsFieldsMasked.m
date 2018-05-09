function [ fields ] = loadDeedsFieldsMasked( fieldIdx,sliceString,fieldsPath,mask,scale)
%% LOADDEEDSFIELDSMASKED: returns the masked deeds registration fields 
%% INPUT:
%   fieldIdx: array of field indicies (#indicies)
%   sliceString: data slice
%   fieldsPath: path to the fields files
%   scale: image scale
%% OUTPUT:
% fields: masked dense registration fields (#image points*2 x #fields)

fields=[];
for i=1:numel(fieldIdx)
    fieldNum=fieldIdx(i);        
    fieldFilename=[fieldsPath,'sl_ss',sliceString,'_tt',num2str(fieldNum),'_displ.mat'];
    load(fieldFilename);
    
    ux=imresize(ux,scale).*scale;
    vx=imresize(vx,scale).*scale;
    
    ux=ux(mask);
    vx=vx(mask);        
    fieldVec=[ux(:); vx(:)];
    fieldVec(1:2:end)=ux(:);
    fieldVec(2:2:end)=vx(:);
    fields=[fields fieldVec];
end

end

