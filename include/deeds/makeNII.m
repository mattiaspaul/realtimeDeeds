
for p=1:199
    str = sprintf('/share/data_tiffy3/ha/MotionEstimation/Results/result_sl01/result_im%03d.mat',p);
    load(str);
    X_result(p,:,:) = im2w;
%    clear im2w;
end

% extra ausführen in NifTI folder
%  nii = make_nii(X_result);
%  view_nii(nii);
