// Tutorial: 
// https://www.caam.rice.edu/~optimization/L1/optseminar/Getreuer-slides-cmex.pdf
// Usage:
// 1. Compile 
//      >> mex computeEdgeCost.cpp
// 2. Call the function 
//      >> computeEdgeCost(R1, R2, xxF, xyF, xzF, yyF, yzF, zzF)

#include <iostream>
#include <math.h>
#include "mex.h"
    
#define R1_in prhs[0]
#define R2_in prhs[1]
#define xxF_in prhs[2]
#define xyF_in prhs[3]
#define xzF_in prhs[4]
#define yyF_in prhs[5]
#define yzF_in prhs[6]
#define zzF_in prhs[7]
#define cost_out plhs[0]

void DoComputation(double *cost, double *R1,  double *R2, double *xxF,
        double *xyF, double *xzF, double *yyF, double *yzF, double *zzF)
{
    double r1[3], r2[3], r3[3];
    
    r1[0] = R1[0] * R2[0] + R1[3] * R2[3] + R1[6] * R2[6];
    r1[1] = R1[0] * R2[1] + R1[3] * R2[4] + R1[6] * R2[7];
    r1[2] = R1[0] * R2[2] + R1[3] * R2[5] + R1[6] * R2[8];
    
    r2[0] = R1[1] * R2[0] + R1[4] * R2[3] + R1[7] * R2[6];
    r2[1] = R1[1] * R2[1] + R1[4] * R2[4] + R1[7] * R2[7];
    r2[2] = R1[1] * R2[2] + R1[4] * R2[5] + R1[7] * R2[8];
    
    r3[0] = R1[2] * R2[0] + R1[5] * R2[3] + R1[8] * R2[6];
    r3[1] = R1[2] * R2[1] + R1[5] * R2[4] + R1[8] * R2[7];
    r3[2] = R1[2] * R2[2] + R1[5] * R2[5] + R1[8] * R2[8];

    double m11 = 0;
    double m12 = 0;
    double m13 = 0;
    double m21 = 0;
    double m22 = 0;
    double m23 = 0;
    double m31 = 0;
    double m32 = 0;
    double m33 = 0;
    
    // m11 = r3*yyF*r3'-2*r3*yzF*r2'+r2*zzF*r2';
    m11 += r3[0]*(r3[0]*yyF[0]+r3[1]*yyF[3]+r3[2]*yyF[6]);
    m11 += r3[1]*(r3[0]*yyF[1]+r3[1]*yyF[4]+r3[2]*yyF[7]);
    m11 += r3[2]*(r3[0]*yyF[2]+r3[1]*yyF[5]+r3[2]*yyF[8]);
    m11 -= 2*r2[0]*(r3[0]*yzF[0]+r3[1]*yzF[3]+r3[2]*yzF[6]);
    m11 -= 2*r2[1]*(r3[0]*yzF[1]+r3[1]*yzF[4]+r3[2]*yzF[7]);
    m11 -= 2*r2[2]*(r3[0]*yzF[2]+r3[1]*yzF[5]+r3[2]*yzF[8]);
    m11 += r2[0]*(r2[0]*zzF[0]+r2[1]*zzF[3]+r2[2]*zzF[6]);
    m11 += r2[1]*(r2[0]*zzF[1]+r2[1]*zzF[4]+r2[2]*zzF[7]);
    m11 += r2[2]*(r2[0]*zzF[2]+r2[1]*zzF[5]+r2[2]*zzF[8]);
    
    // m22 = r1*zzF*r1'-2*r1*xzF*r3'+r3*xxF*r3';
    m22 += r1[0]*(r1[0]*zzF[0]+r1[1]*zzF[3]+r1[2]*zzF[6]);
    m22 += r1[1]*(r1[0]*zzF[1]+r1[1]*zzF[4]+r1[2]*zzF[7]);
    m22 += r1[2]*(r1[0]*zzF[2]+r1[1]*zzF[5]+r1[2]*zzF[8]);
    m22 -= 2*r3[0]*(r1[0]*xzF[0]+r1[1]*xzF[3]+r1[2]*xzF[6]);
    m22 -= 2*r3[1]*(r1[0]*xzF[1]+r1[1]*xzF[4]+r1[2]*xzF[7]);
    m22 -= 2*r3[2]*(r1[0]*xzF[2]+r1[1]*xzF[5]+r1[2]*xzF[8]);
    m22 += r3[0]*(r3[0]*xxF[0]+r3[1]*xxF[3]+r3[2]*xxF[6]);
    m22 += r3[1]*(r3[0]*xxF[1]+r3[1]*xxF[4]+r3[2]*xxF[7]);
    m22 += r3[2]*(r3[0]*xxF[2]+r3[1]*xxF[5]+r3[2]*xxF[8]);
    
    // m33 = r2*xxF*r2'-2*r1*xyF*r2'+r1*yyF*r1';
    m33 += r2[0]*(r2[0]*xxF[0]+r2[1]*xxF[3]+r2[2]*xxF[6]);
    m33 += r2[1]*(r2[0]*xxF[1]+r2[1]*xxF[4]+r2[2]*xxF[7]);
    m33 += r2[2]*(r2[0]*xxF[2]+r2[1]*xxF[5]+r2[2]*xxF[8]);
    m33 -= 2*r2[0]*(r1[0]*xyF[0]+r1[1]*xyF[3]+r1[2]*xyF[6]);
    m33 -= 2*r2[1]*(r1[0]*xyF[1]+r1[1]*xyF[4]+r1[2]*xyF[7]);
    m33 -= 2*r2[2]*(r1[0]*xyF[2]+r1[1]*xyF[5]+r1[2]*xyF[8]);
    m33 += r1[0]*(r1[0]*yyF[0]+r1[1]*yyF[3]+r1[2]*yyF[6]);
    m33 += r1[1]*(r1[0]*yyF[1]+r1[1]*yyF[4]+r1[2]*yyF[7]);
    m33 += r1[2]*(r1[0]*yyF[2]+r1[1]*yyF[5]+r1[2]*yyF[8]);
    
    // m12 = r1*yzF*r3'-r1*zzF*r2'-r3*xyF*r3'+r3*xzF*r2';
    m12 += r3[0]*(r1[0]*yzF[0]+r1[1]*yzF[3]+r1[2]*yzF[6]);
    m12 += r3[1]*(r1[0]*yzF[1]+r1[1]*yzF[4]+r1[2]*yzF[7]);
    m12 += r3[2]*(r1[0]*yzF[2]+r1[1]*yzF[5]+r1[2]*yzF[8]);
    m12 -= r2[0]*(r1[0]*zzF[0]+r1[1]*zzF[3]+r1[2]*zzF[6]);
    m12 -= r2[1]*(r1[0]*zzF[1]+r1[1]*zzF[4]+r1[2]*zzF[7]);
    m12 -= r2[2]*(r1[0]*zzF[2]+r1[1]*zzF[5]+r1[2]*zzF[8]);
    m12 -= r3[0]*(r3[0]*xyF[0]+r3[1]*xyF[3]+r3[2]*xyF[6]);
    m12 -= r3[1]*(r3[0]*xyF[1]+r3[1]*xyF[4]+r3[2]*xyF[7]);
    m12 -= r3[2]*(r3[0]*xyF[2]+r3[1]*xyF[5]+r3[2]*xyF[8]);
    m12 += r2[0]*(r3[0]*xzF[0]+r3[1]*xzF[3]+r3[2]*xzF[6]);
    m12 += r2[1]*(r3[0]*xzF[1]+r3[1]*xzF[4]+r3[2]*xzF[7]);
    m12 += r2[2]*(r3[0]*xzF[2]+r3[1]*xzF[5]+r3[2]*xzF[8]);
    
    // m13 = r2*xyF*r3'-r2*xzF*r2'-r1*yyF*r3'+r1*yzF*r2';
    m13 += r3[0]*(r2[0]*xyF[0]+r2[1]*xyF[3]+r2[2]*xyF[6]);
    m13 += r3[1]*(r2[0]*xyF[1]+r2[1]*xyF[4]+r2[2]*xyF[7]);
    m13 += r3[2]*(r2[0]*xyF[2]+r2[1]*xyF[5]+r2[2]*xyF[8]);
    m13 -= r2[0]*(r2[0]*xzF[0]+r2[1]*xzF[3]+r2[2]*xzF[6]);
    m13 -= r2[1]*(r2[0]*xzF[1]+r2[1]*xzF[4]+r2[2]*xzF[7]);
    m13 -= r2[2]*(r2[0]*xzF[2]+r2[1]*xzF[5]+r2[2]*xzF[8]);
    m13 -= r3[0]*(r1[0]*yyF[0]+r1[1]*yyF[3]+r1[2]*yyF[6]);
    m13 -= r3[1]*(r1[0]*yyF[1]+r1[1]*yyF[4]+r1[2]*yyF[7]);
    m13 -= r3[2]*(r1[0]*yyF[2]+r1[1]*yyF[5]+r1[2]*yyF[8]);
    m13 += r2[0]*(r1[0]*yzF[0]+r1[1]*yzF[3]+r1[2]*yzF[6]);
    m13 += r2[1]*(r1[0]*yzF[1]+r1[1]*yzF[4]+r1[2]*yzF[7]);
    m13 += r2[2]*(r1[0]*yzF[2]+r1[1]*yzF[5]+r1[2]*yzF[8]);

    // m23 = r1*xzF*r2'-r1*yzF*r1'-r3*xxF*r2'+r3*xyF*r1';
    m23 += r2[0]*(r1[0]*xzF[0]+r1[1]*xzF[3]+r1[2]*xzF[6]);
    m23 += r2[1]*(r1[0]*xzF[1]+r1[1]*xzF[4]+r1[2]*xzF[7]);
    m23 += r2[2]*(r1[0]*xzF[2]+r1[1]*xzF[5]+r1[2]*xzF[8]);
    m23 -= r1[0]*(r1[0]*yzF[0]+r1[1]*yzF[3]+r1[2]*yzF[6]);
    m23 -= r1[1]*(r1[0]*yzF[1]+r1[1]*yzF[4]+r1[2]*yzF[7]);
    m23 -= r1[2]*(r1[0]*yzF[2]+r1[1]*yzF[5]+r1[2]*yzF[8]);
    m23 -= r2[0]*(r3[0]*xxF[0]+r3[1]*xxF[3]+r3[2]*xxF[6]);
    m23 -= r2[1]*(r3[0]*xxF[1]+r3[1]*xxF[4]+r3[2]*xxF[7]);
    m23 -= r2[2]*(r3[0]*xxF[2]+r3[1]*xxF[5]+r3[2]*xxF[8]);
    m23 += r1[0]*(r3[0]*xyF[0]+r3[1]*xyF[3]+r3[2]*xyF[6]);
    m23 += r1[1]*(r3[0]*xyF[1]+r3[1]*xyF[4]+r3[2]*xyF[7]);
    m23 += r1[2]*(r3[0]*xyF[2]+r3[1]*xyF[5]+r3[2]*xyF[8]);

 
    double b1 = -m11-m22-m33;
    double b2 = -m13*m13-m23*m23-m12*m12+m11*m22+m11*m33+m22*m33;
    double b3 = m22*m13*m13+m11*m23*m23+m33*m12*m12-m11*m22*m33-2*m12*m23*m13;
    double s = 0.5*(2*b1*b1*b1-9*b1*b2+27*b3);
    double t = sqrt(pow(b1*b1-3*b2, 3));
    double k = pow(t, 1.0/3.0)*cos(acos(s/t)/3);
    double lambda = (-b1-2*k)/3;
    if (lambda < 0)
        lambda = -lambda;
    
    cost[0] = sqrt(lambda);
    
   

}

/* The gateway function. */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

    double *cost, *R1, *R2, *xxF, *xyF, *xzF, *yyF, *yzF, *zzF;
    
    /*** Check Inputs **/
    if(nrhs != 8)
        mexErrMsgTxt("Must have 8 input arguments.");
    if(nlhs != 1)
        mexErrMsgTxt("Must have 1 output arguments.");
    
    if(mxIsComplex(R1_in) || mxGetNumberOfDimensions(R1_in) != 2 
            || mxIsSparse(R1_in) || !mxIsDouble(R1_in))
        mexErrMsgTxt("Sorry! R1 must be a real 2D full double matrix.");
    if(mxIsComplex(R2_in) || mxGetNumberOfDimensions(R2_in) != 2 
            || mxIsSparse(R2_in) || !mxIsDouble(R2_in))
        mexErrMsgTxt("Sorry! R2 must be a real 2D full double matrix.");
    if(mxIsComplex(xxF_in) || mxGetNumberOfDimensions(xxF_in) != 2 
            || mxIsSparse(xxF_in) || !mxIsDouble(xxF_in))
        mexErrMsgTxt("Sorry! xxF must be a real 2D full double matrix.");
    if(mxIsComplex(xyF_in) || mxGetNumberOfDimensions(xyF_in) != 2 
            || mxIsSparse(xyF_in) || !mxIsDouble(xyF_in))
        mexErrMsgTxt("Sorry! xyF must be a real 2D full double matrix.");
    if(mxIsComplex(xzF_in) || mxGetNumberOfDimensions(xzF_in) != 2 
            || mxIsSparse(xzF_in) || !mxIsDouble(xzF_in))
        mexErrMsgTxt("Sorry! xzF must be a real 2D full double matrix.");
    if(mxIsComplex(yyF_in) || mxGetNumberOfDimensions(yyF_in) != 2 
            || mxIsSparse(yyF_in) || !mxIsDouble(yyF_in))
        mexErrMsgTxt("Sorry! yyF must be a real 2D full double matrix.");
    if(mxIsComplex(yzF_in) || mxGetNumberOfDimensions(yzF_in) != 2 
            || mxIsSparse(yzF_in) || !mxIsDouble(yzF_in))
        mexErrMsgTxt("Sorry! yzF must be a real 2D full double matrix.");
    if(mxIsComplex(zzF_in) || mxGetNumberOfDimensions(zzF_in) != 2 
            || mxIsSparse(zzF_in) || !mxIsDouble(zzF_in))
        mexErrMsgTxt("Sorry! zzF must be a real 2D full double matrix.");
    
    
//     
//   p = mxGetScalar(P_IN); /* Get the value of p */
    
    
    /*** Read input ***/
    R1 = mxGetPr(R1_in); 
    R2 = mxGetPr(R2_in); 
    xxF = mxGetPr(xxF_in); 
    xyF = mxGetPr(xyF_in); 
    xzF = mxGetPr(xzF_in); 
    yyF = mxGetPr(yyF_in); 
    yzF = mxGetPr(yzF_in); 
    zzF = mxGetPr(zzF_in); 
    
    /*** Create the output matrix ***/
    cost_out = mxCreateDoubleMatrix(1,1,mxREAL);
    cost = mxGetPr(cost_out);
    
    /* Do the computation */
    DoComputation(cost, R1, R2, xxF, xyF, xzF, yyF, yzF, zzF);
    
    
}
