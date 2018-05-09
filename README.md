# realtimeDeeds
**Code**<p>
Example code for TBME Publication *"Model-based Sparse-to-dense Image Registration for Realtime Respiratory Motion Estimation in Image-guided Interventions"*.<p>


**Datasets**<p>
Datasets used in the experments: <p>
- *2D MRI dataset [1]* used in the code can be found under: https://zenodo.org/record/55345#.WvLINd8xBhF <p>
- *3D MRI dataset [2]* and *3D US dataset [3]* can be found under:<p>
  [2] http://www.vision.ee.ethz.ch/en/datasets/ (4D MRI lung data) <p>
  [3] http://clust.ethz.ch/data.html <p>

*[1]: CF Baumgartner, C Kolbitsch, JR McClelland, D Rueckert, AP King, Autoadaptive motion modelling for MR-based respiratory motion estimation, Medical Image Analysis (2016), http://dx.doi.org/10.1016/j.media.2016.06.005* <p>
*[2]: Boye, D. et al. - Population based modeling of respiratory lung motion and prediction from partial information - Proc. SPIE 8669, Medical Imaging 2013: Image Processing, 86690U (March 13, 2013); doi:10.1117/12.2007076.*<p>
<p>

<b>Landmarks</b><p>
Manually selected landmarks for 2D/3D MRI and 3D US dataset are available as .txt (for 2D MRI and 3D US) and .m (for 3D MRI) files. 
- 2D MRI: each txt-file contains landmark coordinates for one frame.
- 3D MRI: by executing MATLAB script file, 3 matrices (refLM, LM01px/LM04px, frames) are generated. 
- 3D US:  each txt-file contains landmark coordinates for all frames with frame numbers in the first column.
