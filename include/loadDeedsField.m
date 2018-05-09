function [ field ] = loadDeedsField( fieldIdx,sliceString,fieldsPath,scale)

field=[];

fieldFilename=[fieldsPath,'sl_ss',sliceString,'_tt',num2str(fieldIdx),'_displ.mat'];
load(fieldFilename);

field(:,:,1)=imresize(ux,scale).*scale;
field(:,:,2)=imresize(vx,scale).*scale;
end
