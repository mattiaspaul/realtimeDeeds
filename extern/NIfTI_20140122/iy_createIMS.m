% falls die Files noch in DICOM ist, dicm1nii.m ausführen

anz = 200;
IMS = zeros(anz,224,224);
for p=1:anz
    str = sprintf('/share/data_tiffy3/ha/MotionEstimation/Results/nii_sl04/reconstruction_best4_matched_timestep_%04d.nii',p);
    load_nii(str);
    im(:,:) = ans.img(35,:,:);
    IMS(p,:,:) = fliplr(fliplr(im)');
end

sprintf('done.')
nii = make_nii(IMS);
view_nii(nii);

