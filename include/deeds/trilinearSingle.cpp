// Example of a trilinear interpolation as a MEX file
#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
mwSize  ndims, ndims1;
const mwSize *dims, *dims1;
float *in1,*in2, *in3,*in4, *out;
int rows1, cols1, lyrs1, rows2, cols2, lyrs2;
// Check number of arguments
if (nrhs != 4) {
    mexErrMsgTxt("Trilinear requires four input arguments.");
} else if (nlhs > 1) {
    mexErrMsgTxt("One output argument required");
}
// check the type of arguments
if (mxIsSingle(prhs[0]) && mxIsSingle(prhs[1]) && mxIsSingle(prhs[2]) && mxIsSingle(prhs[3]))
{
    ndims1 = mxGetNumberOfDimensions(prhs[0]);
    if (ndims1<3) {
        mexErrMsgTxt("Input must be 3 dimensional");
        return;
    }
    dims1 = mxGetDimensions(prhs[0]);
    rows1 = dims1[0];
    cols1 = dims1[1];
    lyrs1 = dims1[2];

    //in1 = mxGetPr(prhs[0]); // get a pointers to the first input argument
    in1= (float *)mxGetData(prhs[0]);
    in2= (float *)mxGetData(prhs[1]);
    in3= (float *)mxGetData(prhs[2]);
    in4= (float *)mxGetData(prhs[3]);
    
    // create the output
    ndims = mxGetNumberOfDimensions(prhs[1]);
    dims = mxGetDimensions(prhs[1]);
    plhs[0] = mxCreateNumericArray(ndims,dims,mxSINGLE_CLASS,mxREAL);
    // get a pointer to the beginning
    out = (float *)mxGetData(plhs[0]); //mxGetPr(plhs[0]);
    
    rows2 = dims[0];
    cols2 = dims[1];
    if (ndims<3) {
        lyrs2 = 1;
    }
    else {
        lyrs2 = dims[2];
    }
    // compute the trilinear interpolation ...
    int  position, npos, prod=rows1*cols1, prod2=rows2*cols2, izx, izy, izz;
    float row, col, lay, xd, yd, zd, i1, i2, j1, j2, w1, w2;
    for ( int k=0; k<lyrs2; k++) {      //2
        for (int i=0; i<rows2; i++) {   //32
            for (int j=0; j<cols2; j++) {       //32
            // elements of a matrix in the array are stored collumn by column so
            position = j*rows2 + i + k*prod2;
            
            col=in2[position]-1;
            row=in3[position]-1;
            lay=in4[position]-1;
            
            if ( row>=0 & row<=(rows1-1) & col>=0 & col<=(cols1-1) & lay>=0 & lay<=(lyrs1-1) ) {                
                xd = col-int(col); izx=(xd != 0);
                yd = row-int(row); izy=(yd != 0);
                zd = lay-int(lay); izz=(zd != 0);
                
                npos=int(row) + int(col)*rows1 + int(lay)*prod;

                i1 = in1[npos                        ]*(1-yd) + in1[npos + 1*izy                        ]*yd;
                i2 = in1[npos + rows1*izx            ]*(1-yd) + in1[npos + 1*izy + rows1*izx            ]*yd;
                j1 = in1[npos + prod*izz             ]*(1-yd) + in1[npos + 1*izy + prod*izz             ]*yd;
                j2 = in1[npos + rows1*izx  + prod*izz]*(1-yd) + in1[npos + 1*izy + rows1*izx  + prod*izz]*yd;
                
                w1=i1*(1-xd) + i2*xd;
                w2=j1*(1-xd) + j2*xd;

                *(out+position) = w1*(1-zd) + w2*zd;
            }
            else {
                *(out+position) = 0;
            }
            }; //for j
        }; //for i
    };//for k
}
else { mexErrMsgTxt("Arguments must be of type single.");
};
return;
} ;