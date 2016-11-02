#ifndef ISOCUT5_H
#define ISOCUT5_H

struct isocut5_opts {
    bool already_sorted=false;
};

void isocut5(double *dipscore_out,double *cutpoint_out,long N,float *samples,isocut5_opts opts);

/*
 * MCWRAP [ dipscore[1,1], cutpoint[1,1] ] = isocut5_mex(samples[1,N])
 * SET_INPUT N = size(samples,2)
 * SOURCES isocut5.cpp jisotonic.cpp
 * HEADERS isocut5.h jisotonic.h
 */
void isocut5_mex(double *dipscore,double *cutpoint,int N,double *samples);

#endif // ISOCUT5_H
