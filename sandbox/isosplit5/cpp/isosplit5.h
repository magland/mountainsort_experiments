#ifndef ISOSPLIT5_H
#define ISOSPLIT5_H

struct isosplit5_opts {
    float isocut_threshold=1.5;
    int min_cluster_size=10;
};

void isosplit5(int *labels_out,int M, long N,float *X,isosplit5_opts opts);

#endif // ISOSPLIT5_H
