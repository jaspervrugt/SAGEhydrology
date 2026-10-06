#include "mex.h"
#include <cmath>
#include <limits>

static double logistic(double x) {
    if (x >= 0.0) { const double e=std::exp(-x); return 1.0/(1.0+e); }
    const double e=std::exp(x); return e/(1.0+e);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
    if (nrhs != 4 || nlhs > 1)
        mexErrMsgIdAndTxt("SAGE:dualKosugiObjective:Args","Expected eta, q, u, delta.");
    if (!mxIsDouble(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 5 ||
        !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) ||
        mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
        mexErrMsgIdAndTxt("SAGE:dualKosugiObjective:Shape","Invalid inputs.");
    const double *e=mxGetDoubles(prhs[0]), *q=mxGetDoubles(prhs[1]);
    const double *u=mxGetDoubles(prhs[2]);
    const double delta=mxGetScalar(prhs[3]);
    const size_t n=mxGetNumberOfElements(prhs[1]);
    double value=std::numeric_limits<double>::max()/1e100;
    bool valid=std::isfinite(delta) && delta>0.0;
    for (int j=0;j<5;++j) valid=valid && std::isfinite(e[j]) && std::abs(e[j])<40.0;
    if (valid && n>0) {
        const double a1=std::exp(e[0]), b1=std::exp(e[1]);
        const double a2=a1+std::exp(e[2]), b2=std::exp(e[3]);
        const double w=logistic(e[4]), root2=std::sqrt(2.0);
        double sum=0.0;
        for (size_t i=0;i<n;++i) {
            if (!(q[i]>0.0) || !std::isfinite(q[i]) || !std::isfinite(u[i])) { valid=false; break; }
            const double uh=w*0.5*std::erfc(std::log(q[i]/a1)/(root2*b1))+
                (1.0-w)*0.5*std::erfc(std::log(q[i]/a2)/(root2*b2));
            const double r=(u[i]-uh)/delta;
            sum += delta*delta*(std::sqrt(1.0+r*r)-1.0);
        }
        if (valid) value=sum/static_cast<double>(n);
    }
    plhs[0]=mxCreateDoubleScalar(value);
}
