#include "isosplit5.h"
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "isocut5.h"

namespace ns_isosplit5 {
    struct kmeans_opts {
        int num_iterations=5;
    };

    int compute_max(long N,int *labels);
    void kmeans_multistep(int *labels,int M,long N,float *X,int K1,int K2,int K3,kmeans_opts opts);
    void kmeans_maxsize(int *labels,int M,long N,float *X,int maxsize,kmeans_opts opts);
    void compare_clusters(double *dip_score,std::vector<int> *new_labels1,std::vector<int> *new_labels2,int M,long N1,long N2,float *X1,float *X2,float *centroid1,float *centroid2);
}

class isosplit5_data {
public:
    isosplit5_data(int M_in,long N_in,float* X_in) {
        M=M_in;
        N=N_in;
        labels=(int*)malloc(sizeof(int)*N);
        X=X_in;
    }

    virtual ~isosplit5_data() {
        if (labels) free(labels);
        if (centroids) free(centroids);
        if (active_labels_vec) free(active_labels_vec);
    }

    void initialize_labels();

    void setKAndAllocate(int K_in) {
        K=K_in;
        centroids=(double*)malloc(sizeof(double)*M*K);
        for (long i=0; i<M*K; i++)
            centroids[i]=0;
        active_labels_vec=(int*)malloc(sizeof(int)*K);
        for (int k=1; k<=K; k++)
            active_labels_vec[k-1]=1;
    }

    void compute_all_centroids() {
        for (int k=1; k<=K; k++) {
            recompute_centroid(k);
        }
    }
    void recompute_centroid(int k) {
        double count=0;
        double *ptr=&centroids[M*(k-1)];
        for (int m=0; m<M; m++) {
            ptr[m]=0;
        }
        for (long i=0; i<N; i++) {
            if (labels[i]==k) {
                for (int m=0; m<M; m++) {
                    ptr[m]+=X[m+M*i];
                }
                count++;
            }
        }
        if (count) {
            for (int m=0; m<M; m++) {
                ptr[m]/=count;
            }
        }
    }
    void get_pairs_to_compare(std::vector<int> *k1s,std::vector<int> *k2s);
    long compare_pairs(const std::vector<int> &k1s,const std::vector<int> &k2s,float isocut_threshold); //return number of changes
    std::vector<int> get_active_labels();

    float *X;
    int M=0;
    long N=0;
    int K=0;
    int *labels=0;
    double *centroids=0;
    int* active_labels_vec=0;
};

void isosplit5_mex(double *labels_out, int M, int N, double *X)
{
    float *Xf=(float*)malloc(sizeof(float)*M*N);
    int *labelsi=(int*)malloc(sizeof(int)*N);
    for (long i=0; i<M*N; i++)
        Xf[i]=X[i];
    isosplit5_opts opts;
    isosplit5(labelsi,M,N,Xf,opts);
    for (long i=0; i<N; i++)
        labels_out[i]=labelsi[i];
    free(Xf);
    free(labelsi);
}

void isosplit5(int *labels_out,int M, long N,float *X,isosplit5_opts opts) {
    for (long i=0; i<N; i++) {
        labels_out[i]=1;
    }

    isosplit5_data DD(M,N,X);

    DD.initialize_labels();
    DD.compute_all_centroids();

    int max_iterations=500;
    int max_iterations_without_merges=5;

    int iteration_number=1;
    int num_iterations_without_merges=0;
    while (true) {
        iteration_number++;
        if (iteration_number>max_iterations) {
            printf ("isosplit5: Exceeded maximum number of iterations. Breaking.");
            break;
        }

        printf("Number of active labels: %ld\n",DD.get_active_labels().size());

        std::vector<int> k1s,k2s;
        DD.get_pairs_to_compare(&k1s,&k2s);

        printf ("compare %ld pairs\n",k1s.size());
        std::vector<int> old_active_labels=DD.get_active_labels();
        long num_changes=DD.compare_pairs(k1s,k2s,opts.isocut_threshold);
        std::vector<int> new_active_labels=DD.get_active_labels();
        printf ("  %ld changes\n",num_changes);

        if (new_active_labels.size()==old_active_labels.size())
            num_iterations_without_merges++;
        else
            num_iterations_without_merges=0;

        if (num_iterations_without_merges>=max_iterations_without_merges)
            break;
    }

    for (int pass=1; pass<=2; pass++) {
        std::vector<int> active_labels=DD.get_active_labels();
        for (int i1=0; i1<(int)active_labels.size(); i1++) {
            for (int i2=i1+1; i2<(int)active_labels.size(); i2++) {
                int k1=active_labels[i1];
                int k2=active_labels[i2];
                if ((DD.active_labels_vec[k1-1])&&(DD.active_labels_vec[k2-1])) {
                    printf("Number of active labels: %ld\n",DD.get_active_labels().size());
                    printf ("compare %d/%d (pass %d)\n",k1,k2,pass);
                    std::vector<int> k1s,k2s;
                    k1s.push_back(k1);
                    k2s.push_back(k2);
                    long num_changes=DD.compare_pairs(k1s,k2s,opts.isocut_threshold);
                    printf ("  %ld changes\n",num_changes);
                }
            }
        }
    }

    std::vector<int> active_labels=DD.get_active_labels();
    std::vector<int> labels_map(ns_isosplit5::compute_max(N,DD.labels)+1);
    for (int i=0; i<(int)active_labels.size(); i++) {
        labels_map[active_labels[i]]=i+1;
    }
    for (long i=0; i<N; i++) {
        labels_out[i]=labels_map[DD.labels[i]];
    }
}

namespace ns_isosplit5 {
    int compute_max(long N,int *labels) {
        if (N==0) return 0;
        int ret=labels[0];
        for (long i=0; i<N; i++) {
            if (labels[i]>ret)
                ret=labels[i];
        }
        return ret;
    }

    void kmeans_initialize(double *centroids,int M,long N,int K,float *X) {
        std::vector<int> used(N);
        for (long i=0; i<N; i++)
            used[i]=0;
        for (int k=0; k<K; k++)
            used[k]=1;
        for (int k=0; k<K; k++) {
            long ii=rand()%N;
            int tmp=used[k];
            used[k]=used[ii];
            used[ii]=tmp;
        }
        std::vector<long> inds;
        for (long i=0; i<N; i++) {
            if (used[i])
                inds.push_back(i);
        }
        for (int k=0; k<(int)inds.size(); k++) {
            for (int m=0; m<M; m++) {
                centroids[m+M*k]=X[m+M*inds[k]];
            }
        }
    }
    double compute_dist(int M,float *X,double *Y) {
        double sumsqr=0;
        for (int m=0; m<M; m++) {
            double val=X[m]-Y[m];
            sumsqr+=val*val;
        }
        return sqrt(sumsqr);
    }

    int kmeans_assign2(int M,int K,float *X0,double *centroids) {
        int ret=0;
        double best_dist=0;
        for (int k=1; k<=K; k++) {
            double dist=compute_dist(M,X0,&centroids[M*(k-1)]);
            if ((ret==0)||(dist<best_dist)) {
                best_dist=dist;
                ret=k;
            }
        }
        return ret;
    }
    void kmeans_assign(int *labels,int M,long N,int K,float *X,double *centroids) {
        for (long i=0; i<N; i++) {
            labels[i]=kmeans_assign2(M,K,&X[M*i],centroids);
        }
    }
    void kmeans_centroids(double *centroids,int M,long N,int K,float *X,int *labels) {
        std::vector<long> counts(K);
        for (int k=1; k<=K; k++) {
            counts[k-1]=0;
            for (int m=0; m<M; m++) {
                centroids[m+(k-1)*M]=0;
            }
        }
        for (long i=0; i<N; i++) {
            int k=labels[i];
            for (int m=0; m<M; m++) {
                centroids[m+(k-1)*M]+=X[m+i*M];
            }
            counts[k-1]++;
        }
        for (int k=1; k<=K; k++) {
            if (counts[k-1]) {
                for (int m=0; m<M; m++) {
                    centroids[m+(k-1)*M]/=counts[k-1];
                }
            }
        }
    }

    void kmeans(int *labels,int M,long N,float *X,int K,kmeans_opts opts) {
        if (K>N) K=N;
        double *centroids=(double*)malloc(sizeof(double)*M*K);
        kmeans_initialize(centroids,M,N,K,X);
        for (int it=1; it<=opts.num_iterations; it++) {
            kmeans_assign(labels,M,N,K,X,centroids);
            kmeans_centroids(centroids,M,N,K,X,labels);
        }
        kmeans_assign(labels,M,N,K,X,centroids);
        free(centroids);
    }

    void extract_subarray(float *X_sub,int M,float *X,const std::vector<long> &inds) {
        for (long i=0; i<(long)inds.size(); i++) {
            for (int m=0; m<M; m++) {
                X_sub[m+i*M]=X[m+inds[i]*M];
            }
        }
    }

    void kmeans_maxsize(int *labels,int M,long N,float *X,int maxsize,kmeans_opts opts) {
        if (N<=maxsize) {
            for (long i=0; i<N; i++)
                labels[i]=1;
            return;
        }
        int K=ceil(N*1.0/maxsize);
        int *labels1=(int*)malloc(sizeof(int)*N);
        kmeans(labels1,M,N,X,K,opts);
        int L1=compute_max(N,labels1);
        int current_max_k=0;
        for (int k=1; k<=L1; k++) {
            std::vector<long> inds_k;
            for (long i=0; i<N; i++) {
                if (labels1[i]==k)
                    inds_k.push_back(i);
            }
            if (inds_k.size()>0) {
                float *X2=(float*)malloc(sizeof(float)*M*inds_k.size());
                int *labels2=(int*)malloc(sizeof(int)*inds_k.size());
                extract_subarray(X2,M,X,inds_k);
                kmeans_maxsize(labels2,M,inds_k.size(),X2,maxsize,opts);
                for (long j=0; j<(long)inds_k.size(); j++) {
                    labels[inds_k[j]]=current_max_k+labels2[j];
                }
                current_max_k+=compute_max(inds_k.size(),labels2);
                free(X2);
                free(labels2);
            }
        }
        free(labels1);
    }

    void kmeans_multistep(int *labels,int M,long N,float *X,int K1,int K2,int K3,kmeans_opts opts) {
        if (K2>1) {
            int *labels1=(int*)malloc(sizeof(int)*N);
            for (long i=0; i<N; i++)
                labels1[i]=0;
            kmeans_multistep(labels1,M,N,X,K2,K3,0,opts);
            int L1=compute_max(N,labels1);
            int current_max_k=0;
            for (int k=1; k<=L1; k++) {
                std::vector<long> inds_k;
                for (long i=0; i<N; i++) {
                    if (labels1[i]==k)
                        inds_k.push_back(i);
                }
                if (inds_k.size()>0) {
                    float *X2=(float*)malloc(sizeof(float)*M*inds_k.size());
                    int *labels2=(int*)malloc(sizeof(int)*inds_k.size());
                    extract_subarray(X2,M,X,inds_k);
                    kmeans_multistep(labels2,M,inds_k.size(),X2,K1,0,0,opts);
                    for (long j=0; j<(long)inds_k.size(); j++) {
                        labels[inds_k[j]]=current_max_k+labels2[j];
                    }
                    current_max_k+=compute_max(inds_k.size(),labels2);
                    free(X2);
                    free(labels2);
                }
            }
            free(labels1);
        }
        else {
            kmeans(labels,M,N,X,K1,opts);
        }
    }

    double dot_product(long N,float *X,float *Y) {
        double ret=0;
        for (long i=0; i<N; i++) {
            ret+=X[i]*Y[i];
        }
        return ret;
    }

    void normalize_vector(int N,float *V) {
        double norm=sqrt(dot_product(N,V,V));
        if (!norm) return;
        for (long i=0; i<N; i++)
            V[i]/=norm;
    }

    void compare_clusters(double *dip_score,std::vector<int> *new_labels1,std::vector<int> *new_labels2,int M,long N1,long N2,float *X1,float *X2,double *centroid1,double *centroid2) {
        float *V=(float*)malloc(sizeof(float)*M);
        float *projection=(float*)malloc(sizeof(float)*(N1+N2));
        for (int m=0; m<M; m++) {
            V[m]=centroid2[m]-centroid1[m];
        }
        normalize_vector(M,V);
        for (long i=0; i<N1; i++) {
            projection[i]=dot_product(M,V,&X1[M*i]);
        }
        for (long i=0; i<N2; i++) {
            projection[N1+i]=dot_product(M,V,&X2[M*i]);
        }
        isocut5_opts icopts;
        icopts.already_sorted=false;
        double cutpoint;
        isocut5(dip_score,&cutpoint,N1+N2,projection,icopts);
        new_labels1->resize(N1);
        new_labels2->resize(N2);
        for (long i=0; i<N1; i++) {
            if (projection[i]<cutpoint) (*new_labels1)[i]=1;
            else (*new_labels1)[i]=2;
        }
        for (long i=0; i<N2; i++) {
            if (projection[N1+i]<cutpoint) (*new_labels2)[i]=1;
            else (*new_labels2)[i]=2;
        }
        free(projection);
        free(V);
    }
}

void get_pairs_to_compare3(std::vector<int> *i1s, std::vector<int> *i2s,int M,int N,double *centroids) {
    float distances[N][N];
    int used[N];
    for (int i=0; i<N; i++)
        used[i]=0;
    for (int i=0; i<N; i++) {
        for (int j=i; j<N; j++) {
            double sumsqr=0;
            for (int m=0; m<M; m++) {
                float diff0=centroids[m+M*i]-centroids[m+M*j];
                sumsqr+=diff0*diff0;
            }
            distances[i][j]=distances[j][i]=sqrt(sumsqr);
        }
    }

    bool something_changed=true;
    while (something_changed) {
        something_changed=false;
        int closest[N];
        for (int i=0; i<N; i++)
            closest[i]=-1;
        for (int i=0; i<N; i++) {
            double best_distance=-1;
            if (!used[i]) {
                for (int j=0; j<N; j++) {
                    if ((!used[j])&&(j!=i)) {
                        double dist=distances[i][j];
                        if ((best_distance<0)||(dist<best_distance)) {
                            best_distance=dist;
                            closest[i]=j;
                        }
                    }
                }
            }
        }
        for (int i=0; i<N; i++) {
            if (!used[i]) {
                if ((closest[i]>=0)&&(closest[i]>i)&&(closest[closest[i]]==i)) { //mutual nearest neighbor among those that have not yet been used
                    i1s->push_back(i);
                    i2s->push_back(closest[i]);
                    used[i]=1;
                    used[closest[i]]=1;
                    something_changed=true;
                }
            }
        }
    }
}

void get_pairs_to_compare2(std::vector<int> *i1s, std::vector<int> *i2s,int M,int N,double *centroids) {
    float *centroidsf=(float*)malloc(sizeof(float)*M*N);
    for (long i=0; i<M*N; i++)
        centroidsf[i]=centroids[i];
    int *groups=(int*)malloc(sizeof(int)*N);
    int maxsize=1000;
    ns_isosplit5::kmeans_opts oo;
    ns_isosplit5::kmeans_maxsize(groups,M,N,centroidsf,maxsize,oo);
    int num_groups=ns_isosplit5::compute_max(N,groups);
    printf("num_groups=%d, num_centroids=%d\n",num_groups,N);

    for (int group_number=1; group_number<=num_groups; group_number++) {
        std::vector<long> inds_group;
        for (long i=0; i<N; i++)
            if (groups[i]==group_number)
                inds_group.push_back(i);
        long N0=inds_group.size();
        if (N0>0) {
            float *centroids0f=(float*)malloc(sizeof(float)*M*N0);
            double *centroids0=(double*)malloc(sizeof(double)*M*N0);
            ns_isosplit5::extract_subarray(centroids0f,M,centroidsf,inds_group);
            for (long i=0; i<M*N0; i++)
                centroids0[i]=centroids0f[i];
            std::vector<int> i1s0,i2s0;
            get_pairs_to_compare3(&i1s0,&i2s0,M,N0,centroids0);
            for (long jj=0; jj<(long)i1s0.size(); jj++) {
                i1s->push_back(inds_group[i1s0[jj]]);
                i2s->push_back(inds_group[i2s0[jj]]);
            }
            free(centroids0f);
            free(centroids0);
        }
    }

    free(groups);
    free(centroidsf);
}

void isosplit5_data::initialize_labels()
{
    ns_isosplit5::kmeans_opts oo;
    ns_isosplit5::kmeans_multistep(labels,M,N,X,10,10,10,oo);
    //ns_isosplit5::kmeans_maxsize(labels,M,N,X,10000,oo);
    setKAndAllocate(ns_isosplit5::compute_max(N,labels));
}

void isosplit5_data::get_pairs_to_compare(std::vector<int> *k1s, std::vector<int> *k2s)
{
    k1s->clear();
    k2s->clear();
    std::vector<int> active_labels=get_active_labels();
    if (active_labels.size()<=1) return;
    int A=active_labels.size();
    double *active_centroids=(double*)malloc(sizeof(double)*M*A);
    for (int i=0; i<A; i++) {
        int k=active_labels[i];
        for (int m=0; m<M; m++) {
            active_centroids[m+M*i]=centroids[m+M*(k-1)];
        }
    }
    std::vector<int> i1s,i2s;
    get_pairs_to_compare2(&i1s,&i2s,M,A,active_centroids);
    for (int i=0; i<(int)i1s.size(); i++) {
        k1s->push_back(active_labels[i1s[i]]);
        k2s->push_back(active_labels[i2s[i]]);
    }
    free(active_centroids);
}

long isosplit5_data::compare_pairs(const std::vector<int> &k1s, const std::vector<int> &k2s,float isocut_threshold)
{
    long num_changes=0;
    for (int i=0; i<(int)k1s.size(); i++) {
        int k1=k1s[i];
        int k2=k2s[i];
        std::vector<long> inds_1,inds_2;
        for (long i=0; i<this->N; i++) {
            if (this->labels[i]==k1)
                inds_1.push_back(i);
            if (this->labels[i]==k2)
                inds_2.push_back(i);
        }
        float *X1=(float*)malloc(sizeof(float)*M*inds_1.size());
        float *X2=(float*)malloc(sizeof(float)*M*inds_2.size());
        ns_isosplit5::extract_subarray(X1,M,X,inds_1);
        ns_isosplit5::extract_subarray(X2,M,X,inds_2);
        double dip_score=0;
        std::vector<int> new_labels1,new_labels2;
        ns_isosplit5::compare_clusters(&dip_score,&new_labels1,&new_labels2,M,inds_1.size(),inds_2.size(),X1,X2,&centroids[M*(k1-1)],&centroids[M*(k2-1)]);
        if (k1s.size()==1) printf("dip score = %g\n",dip_score);
        if (dip_score<isocut_threshold) {
            //merge
            for (long a=0; a<(long)inds_2.size(); a++) {
                this->labels[inds_2[a]]=k1;
            }
            this->active_labels_vec[k2-1]=0;
            this->recompute_centroid(k1);
            num_changes+=inds_2.size();
        }
        else {
            //redistribute
            for (long a=0; a<(long)inds_1.size(); a++) {
                if (new_labels1[a]==2) {
                    this->labels[inds_1[a]]=k2;
                    num_changes++;
                }
            }
            for (long a=0; a<(long)inds_2.size(); a++) {
                if (new_labels2[a]==1) {
                    this->labels[inds_2[a]]=k1;
                    num_changes++;
                }
            }
            this->recompute_centroid(k1);
            this->recompute_centroid(k2);
        }
        free(X1);
        free(X2);
    }
    return num_changes;
}

std::vector<int> isosplit5_data::get_active_labels()
{
    std::vector<int> ret;
    for (int i=0; i<K; i++) {
        if (active_labels_vec[i])
            ret.push_back(i+1);
    }
    return ret;
}


