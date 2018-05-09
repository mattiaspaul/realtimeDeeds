#include "mex.h"
#include <math.h>
#include <iostream>
//#include <sys/time.h>
#include <vector>
#include <algorithm>


using namespace std;
#define printf mexPrintf

float beta=1;

#include "primsMST2d.h"

void dt1sq(float *val,int* ind,int len,float offset,int k,int* v,float* z,float* f,int* ind1){
	float INF=1e10;
    
	int j=0;
	z[0]=-INF;
	z[1]=INF;
	v[0]=0;
	for(int q=1;q<len;q++){
		float s=((val[q*k]+(q+offset)*(q+offset))-(val[v[j]*k]+(v[j]+offset)*(v[j]+offset)))/(2.0*(q-v[j]));
        while(s<=z[j]){
			j--;
            s=((val[q*k]+(q+offset)*(q+offset))-(val[v[j]*k]+(v[j]+offset)*(v[j]+offset)))/(2.0*(q-v[j]));
		}
		j++;
		v[j]=q;
		z[j]=s;
		z[j+1]=INF;
	}
	j=0;
	for(int q=0;q<len;q++){
		f[q]=val[q*k]; 
		ind1[q]=ind[q*k];
	}
	for(int q=0;q<len;q++){
		while(z[j+1]<q){
			j++;
		}
		ind[q*k]=ind1[v[j]];
        val[q*k]=(q-v[j]-offset)*(q-v[j]-offset)+f[v[j]];
	}
    
}

void dt2x(float* r,int* indr,int rl,float dx,float dy){
	//rl is length of one side
	for(int i=0;i<rl*rl;i++){
		indr[i]=i;
	}
	int* v=new int[rl]; //slightly faster if not intitialised in each loop
	float* z=new float[rl+1];
	float* f=new float[rl];
	int* i1=new int[rl];
	
    for(int i=0;i<rl;i++){
        dt1sq(r+i,indr+i,rl,-dx,rl,v,z,f,i1);//);
    }
	
    for(int j=0;j<rl;j++){
        dt1sq(r+j*rl,indr+j*rl,rl,-dy,1,v,z,f,i1);//);
    }
	
    delete []i1;
	delete []f;
    
	delete []v;
	delete []z;
    
	
}

void dt2xNaive(float* r,int* indr,int rl,float dx,float dy){
    float* data2=new float[rl*rl];
    for(int i=0;i<rl*rl;i++){
        data2[i]=r[i];
    }
    for(int x1=0;x1<rl;x1++){
        for(int y1=0;y1<rl;y1++){
            float minval=1e20; int minind=-1;
            for(int x2=0;x2<rl;x2++){
                for(int y2=0;y2<rl;y2++){
                    float data=data2[y2+x2*rl];
                    float regularisation=pow(x1-x2+dx,2)+pow(y1-y2+dy,2);
                    float cost=data+regularisation;
                    if(cost<minval){
                        minval=cost;
                        minind=y2+x2*rl;
                    }
                }
            }
            r[y1+x1*rl]=minval;
            indr[y1+x1*rl]=minind;
        }
    }
    delete data2;
}


void MSTmarginals(float* marginals,int* ordered,int* parents,float* u0,float* v0,int sz,int len2,float quant){
 
	//timeval time1,time2;
    

    int len1=sqrt(len2);
    
	float *cost1=new float[len2];
    int *indr=new int[len2];

	float *vals=new float[len2];
    
	//gettimeofday(&time1, NULL);
	
    
	for(int i=0;i<len2;i++){
		cost1[i]=0;
	}
	
	float* message=new float[sz*len2];
	for(int i=0;i<sz*len2;i++){
		message[i]=0.0;
	}
    
    	//calculate mst-cost
	for(int i=(sz-1);i>0;i--){ //do for each control point
		
		int ochild=ordered[i];
		int oparent=parents[ordered[i]];
		
		for(int l=0;l<len2;l++){
			cost1[l]=marginals[ochild+l*sz];
		}
        //for incremental regularisation
        float dx1=(u0[oparent]-u0[ochild])/(float)quant;
		float dy1=(u0[oparent]-u0[ochild])/(float)quant;
		//fast distance transform
		dt2x(cost1,indr,len1,dx1,dy1);
		
		for(int l=0;l<len2;l++){
			message[ochild+l*sz]=cost1[l];
			marginals[oparent+l*sz]+=cost1[l];
		}
	}
	
	for(int i=1;i<sz;i++){
		int ochild=ordered[i];
		int oparent=parents[ordered[i]];
		for(int l=0;l<len2;l++){
			cost1[l]=marginals[oparent+l*sz]-message[ochild+l*sz]+message[oparent+l*sz];
		}
        //for incremental regularisation (now parent->child)
        float dx1=(u0[ochild]-u0[oparent])/(float)quant;
		float dy1=(u0[ochild]-u0[oparent])/(float)quant;
		//fast distance transform
		dt2x(cost1,indr,len1,dx1,dy1);
        
		for(int l=0;l<len2;l++){
			message[ochild+l*sz]=cost1[l];
		}
	}
	
    
	//printf("backward pass finished\n");
	
	for(int i=0;i<sz*len2;i++){
		marginals[i]+=message[i];
	}
	
	delete message;
	
	delete cost1;
	delete vals;
	
}
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
//timeval time1,time2;
	float* data=(float*)mxGetData(prhs[0]);
	float* im1=(float*)mxGetData(prhs[1]);
    int step1=(int)mxGetScalar(prhs[2]);
    
    float* u0=(float*)mxGetData(prhs[3]);
	float* v0=(float*)mxGetData(prhs[4]);
    float quant=(float)mxGetScalar(prhs[5]);

    //-- dimension of the data matrix and im1
	const mwSize* dims1=mxGetDimensions(prhs[0]);
	int m=dims1[0]; int n=dims1[1]; int len2=dims1[2];
    const mwSize* dims2=mxGetDimensions(prhs[1]);
	int m2=dims2[0]; int n2=dims2[1];
    
    int sz=m*n;
	//printf("sz: %d, len2: %d, step1: %d\n",sz,len2,step1);
	//gettimeofday(&time1, NULL);

	int* parents=new int[sz];
    int* ordered=new int[sz];
    primsGraph(im1,ordered,parents,m2,n2,step1); //m2xn2 (size of image)

 //   printf("ordered max: %d, min %d\nparents max: %d, min %d\n",*max_element(ordered,ordered+sz)
 //          ,*min_element(ordered,ordered+sz),*max_element(parents,parents+sz),*min_element(parents,parents+sz));
	int* dimsout=new int[3];
	dimsout[0]=m; dimsout[1]=n; dimsout[2]=len2;
	
    plhs[0]=mxCreateNumericArray(3,dimsout,mxSINGLE_CLASS,mxREAL);
	float* marginals=(float*)mxGetData(plhs[0]);
	
	for(int i=0;i<sz*len2;i++){
        marginals[i]=data[i];
    }
	
	
    MSTmarginals(marginals,ordered,parents,u0,v0,sz,len2,quant);

	//gettimeofday(&time2, NULL);
	//float timeA=time2.tv_sec+time2.tv_usec/1e6-(time1.tv_sec+time1.tv_usec/1e6);
	//printf("Computation time deeds-MST reg: %f secs.\n",timeA);
	
	
	
    return;
}