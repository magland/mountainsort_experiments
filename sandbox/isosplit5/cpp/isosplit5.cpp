#include "isosplit5.h"
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <math.h>
#include "isocut5.h"

typedef std::vector<std::vector<int>> intarray2d;
void alloc(intarray2d &X,int N1,int N2) {
    X.resize(N1);
    for (int i=0; i<N1; i++) {
        X[i].resize(N2);
    }
}

namespace ns_isosplit5 {
    struct kmeans_opts {
        int num_iterations=0;
    };

    int compute_max(long N,int *labels);
    long compute_max(long N,long *inds);
    void kmeans_multistep(int *labels,int M,long N,float *X,int K1,int K2,int K3,kmeans_opts opts);
    void kmeans_maxsize(int *labels,int M,long N,float *X,int maxsize,kmeans_opts opts);
    void compare_clusters(double *dip_score,std::vector<int> *new_labels1,std::vector<int> *new_labels2,int M,long N1,long N2,float *X1,float *X2,float *centroid1,float *centroid2);
    void compute_centroids(float *centroids,int M,long N,int Kmax,float *X,int *labels);
    void get_pairs_to_compare(std::vector<int> *inds1,std::vector<int> *inds2,int M,int K,float *active_centroids,const intarray2d &active_comparisons_made);
    void compare_pairs(std::vector<int> *clusters_changed,long *total_num_label_changes,int M,long N,float *X,int *labels,const std::vector<int> &inds1,const std::vector<int> &inds2,const isosplit5_opts &opts); //the labels are updated
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

struct parcelate2_opts {
    bool final_reassign=false; //not yet implemented
};

struct p2_parcel {
    std::vector<long> indices;
    std::vector<float> centroid;
    double radius;
};

std::vector<float> p2_compute_centroid(int M,float *X,const std::vector<long> &indices) {
    std::vector<double> ret(M);
    double count=0;
    for (int m=0; m<M; m++) {
        ret[m]=0;
    }
    for (long i=0; i<(long)indices.size(); i++) {
        for (int m=0; m<M; m++) {
            ret[m]+=X[m+M*indices[i]];
        }
        count++;
    }
    if (count) {
        for (int m=0; m<M; m++) {
            ret[m]/=count;
        }
    }
    std::vector<float> retf(M);
    for (int m=0; m<M; m++)
        retf[m]=ret[m];
    return retf;
}

double p2_compute_max_distance(const std::vector<float> &centroid,int M,float *X,const std::vector<long> &indices) {
    double max_dist=0;
    for (long i=0; i<(long)indices.size(); i++) {
        double dist=0;
        for (int m=0; m<M; m++) {
            double val=centroid[m]-X[m+M*indices[i]];
            dist+=val*val;
        }
        dist=sqrt(dist);
        if (dist>max_dist) max_dist=dist;
    }
    return max_dist;
}

std::vector<long> p2_randsample(long N,long K) {
    std::vector<long> inds;
    for (int a=0; a<K; a++)
        inds.push_back(a);
    return inds;
    /*
    if (K>N) K=N;std::vector<long> inds;
    std::vector<int> used(N);
    for (long i=0; i<N; i++)
        used[i]=0;
    for (long k=0; k<K; k++)
        used[k]=1;
    for (long k=0; k<K; k++) {
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
    return inds;
    */
}

void parcelate2(int *labels,int M,long N,float *X,int target_parcel_size,int target_num_parcels,const parcelate2_opts &p2opts) {
    printf("parcelate2...\n");
    std::vector<p2_parcel> parcels;

    for (long i=0; i<N; i++)
        labels[i]=1;

    p2_parcel P;
    P.indices.resize(N);
    for (long i=0; i<N; i++)
        P.indices[i]=i;
    P.centroid=p2_compute_centroid(M,X,P.indices);
    P.radius=p2_compute_max_distance(P.centroid,M,X,P.indices);
    parcels.push_back(P);

    printf("Radius of initial cluster = %g\n",P.radius);

    int split_factor=3; // split factor around 2.71 is in a sense ideal

    double target_radius;
    while ((long)parcels.size()<target_num_parcels) {
        bool candidate_found=false;
        for (int i=0; i<(int)parcels.size(); i++) {
            if ((parcels[i].radius>0)&&(parcels[i].indices.size()>(unsigned long)target_parcel_size))
                candidate_found=true;
        }
        if (!candidate_found) {
            printf("No more candidates found for splitting.\n");
            // nothing else will ever be split
            break;
        }

        target_radius=0;
        for (int i=0; i<(int)parcels.size(); i++) {
            if (parcels[i].indices.size()>target_parcel_size) {
                double tmp=parcels[i].radius*0.95;
                if (tmp>target_radius) target_radius=tmp;
            }
        }
        if (target_radius==0) {
            printf("Unexpected target radius of zero.\n");
            break;
        }

        int p_index=0;
        while (p_index<(int)parcels.size()) {
            std::vector<long> inds=parcels[p_index].indices;
            double rad=parcels[p_index].radius;
            long sz=parcels[p_index].indices.size();
            if ((sz>target_parcel_size)&&(rad>=target_radius)) {
                printf("sz = %ld\n",sz);
                std::vector<int> assignments(inds.size());
                std::vector<long> iii=p2_randsample(sz,split_factor);
                for (long i=0; i<(long)inds.size(); i++) {
                    int best_pt=-1;
                    double best_dist=0;
                    for (int j=0; j<(int)iii.size(); j++) {
                        double dist=0;
                        for (int m=0; m<M; m++) {
                            double val=X[m+M*inds[iii[j]]]-X[m+M*inds[i]];
                            dist+=val*val;
                        }
                        dist=sqrt(dist);
                        if ((best_pt<0)||(dist<best_dist)) {
                            best_dist=dist;
                            best_pt=j;
                        }
                    }
                    assignments[i]=best_pt;
                }
                parcels[p_index].indices.clear();
                for (long i=0; i<(long)inds.size(); i++) {
                    if (assignments[i]==0) {
                        parcels[p_index].indices.push_back(inds[i]);
                        labels[inds[i]]=p_index+1;
                    }
                }
                parcels[p_index].centroid=p2_compute_centroid(M,X,parcels[p_index].indices);
                parcels[p_index].radius=p2_compute_max_distance(parcels[p_index].centroid,M,X,parcels[p_index].indices);
                for (int jj=1; jj<(int)iii.size(); jj++) {
                    p2_parcel PP;
                    for (long i=0; i<(long)inds.size(); i++) {
                        if (assignments[i]==jj) {
                            PP.indices.push_back(inds[i]);
                            labels[inds[i]]=parcels.size()+1;
                        }
                    }
                    PP.centroid=p2_compute_centroid(M,X,PP.indices);
                    PP.radius=p2_compute_max_distance(PP.centroid,M,X,PP.indices);
                    if (PP.indices.size()>0)
                        parcels.push_back(PP);
                    else
                        printf("Unexpected problem. New parcel has no points -- original size = %ld.\n",sz);
                }
                if ((long)parcels[p_index].indices.size()==sz) {
                    printf("Warning: Size did not change after splitting parcel.\n");
                    p_index++;
                }
            }
            else {
                p_index++;
            }
        }
    }

     //final reassign not yet implemented
    if (p2opts.final_reassign) {
        //centroids=get_parcel_centroids(parcels);
        //labels=knnsearch(centroids',X','K',1)';
    }
}

void isosplit5(int *labels,int M, long N,float *X,isosplit5_opts opts) {

    // compute the initial clusters
    printf("Computing initial clusters (M = %d, N=%ld)...\n",M,N);
    int target_parcel_size=opts.min_cluster_size;
    int target_num_parcels=opts.K_init;
    // !! important not to do a final reassign because then the shapes will not be conducive to isosplit iterations -- hexagons are not good for isosplit!
    parcelate2_opts p2opts;
    p2opts.final_reassign=false;
    parcelate2(labels,M,N,X,target_parcel_size,target_num_parcels,p2opts);
    int Kmax=ns_isosplit5::compute_max(N,labels);

    printf("Kmax = %d\n",Kmax);

    printf("Computing centroids...\n");
    float *centroids=(float *)malloc(sizeof(float)*M*Kmax);
    ns_isosplit5::compute_centroids(centroids,M,N,Kmax,X,labels);

    // The active labels are those that are still being used -- for now, everything is active
    int active_labels_vec[Kmax];
    for (long i=0; i<Kmax; i++)
        active_labels_vec[i]=1;
    std::vector<int> active_labels;
    for (int i=0; i<Kmax; i++)
        active_labels.push_back(i+1);

    printf("Number of active labels = %ld\n",active_labels.size());

    // Repeat while something has been merged in the pass
    bool final_pass=false; // plus we do one final pass at the end
    intarray2d comparisons_made; // Keep a matrix of comparisons that have been made in this pass
    alloc(comparisons_made,Kmax,Kmax);
    for (int i1=0; i1<Kmax; i1++)
        for (int i2=0; i2<Kmax; i2++)
            comparisons_made[i1][i2]=0;
    while (true) { //passes
        bool something_merged=false; //Keep track of whether something has merged in this pass. If not, do a final pass.
        int clusters_changed_vec[Kmax]; //Keep track of the clusters that have changed in this pass so that we can update the comparisons_made matrix at the end
        for (int i=0; i<Kmax; i++) clusters_changed_vec[i]=0;
        int iteration_number=0;
        while (true) { //iterations
            iteration_number++;
            if (iteration_number>opts.max_iterations_per_pass) {
                printf("Warning: max iterations per pass exceeded.\n");
                break;
            }
            printf("Iteration %d\n",iteration_number);

            printf("num active labels = %ld\n",active_labels.size());

            if (active_labels.size()>0) {
                // Create an array of active centroids and comparisons made, for determining the pairs to compare
                float *active_centroids=(float *)malloc(sizeof(float)*M*active_labels.size());
                for (int i=0; i<(int)active_labels.size(); i++) {
                    for (int m=0; m<M; m++) {
                        active_centroids[m+M*i]=centroids[m+M*(active_labels[i]-1)];
                    }
                }
                intarray2d active_comparisons_made;
                alloc(active_comparisons_made,active_labels.size(),active_labels.size());
                for (int i1=0; i1<(int)active_labels.size(); i1++) {
                    for (int i2=0; i2<(int)active_labels.size(); i2++) {
                        active_comparisons_made[i1][i2]=comparisons_made[active_labels[i1]-1][active_labels[i2]-1];
                    }
                }

                // Find the pairs to compare on this iteration
                // These will be closest pairs of active clusters that have not yet
                // been compared in this pass
                std::vector<int> inds1,inds2;
                ns_isosplit5::get_pairs_to_compare(&inds1,&inds2,M,active_labels.size(),active_centroids,active_comparisons_made);
                std::vector<int> inds1b,inds2b; //remap the clusters to the original labeling
                for (int i=0; i<(int)inds1.size(); i++) {
                    inds1b.push_back(active_labels[inds1[i]-1]);
                    inds2b.push_back(active_labels[inds2[i]-1]);
                }

                // If we didn't find any, break from this iteration
                if (inds1b.size()==0) {
                    break;
                }

                printf("Comparing %ld cluster pairs...\n",inds1.size());
                // Actually compare the pairs -- in principle this operation could be parallelized
                std::vector<int> clusters_changed;
                long total_num_label_changes=0;
                ns_isosplit5::compare_pairs(&clusters_changed,&total_num_label_changes,M,N,X,labels,inds1b,inds2b,opts); //the labels are updated
                printf("Number of clusters changed is %ld (%ld labels).\n",clusters_changed.size(),total_num_label_changes);
                for (int i=0; i<(int)clusters_changed.size(); i++) clusters_changed_vec[clusters_changed[i]]=1;

                printf(".\n");

                // Update which comparisons have been made
                for (int j=0; j<(int)inds1b.size(); j++) {
                    comparisons_made[inds1b[j]-1][inds2b[j]-1]=1;
                    comparisons_made[inds2b[j]-1][inds1b[j]-1]=1;
                }

                // Recompute the centers -- note: maybe this should only apply to those that changed? That would speed things up
                ns_isosplit5::compute_centroids(centroids,M,N,Kmax,X,labels);

                // For diagnostics
                printf("total num label changes = %ld\n",total_num_label_changes);

                // Determine whether something has merged and update the active labels
                for (int i=0; i<Kmax; i++) active_labels_vec[i]=0;
                for (long i=0; i<N; i++)
                    active_labels_vec[labels[i]-1]=1;
                std::vector<int> new_active_labels;
                for (int i=0; i<Kmax; i++)
                    if (active_labels_vec[i])
                        new_active_labels.push_back(i+1);
                if (new_active_labels.size()<active_labels.size())
                    something_merged=true;
                active_labels=new_active_labels;

                free(active_centroids);
            }
            //break;
        }



        // zero out the comparisons made matrix only for those that have changed
        printf("Clusters changed -------------------------------- : ");
        for (int i=0; i<Kmax; i++) {
            if (clusters_changed_vec[i]) {
                printf("%d ",i+1);
                for (int j=0; j<Kmax; j++) {
                    comparisons_made[i][j]=0;
                    comparisons_made[j][i]=0;
                }
            }
        }
        printf("\n");

        if (something_merged) final_pass=false;
        if (final_pass) break; // This was the final pass and nothing has merged
        if (!something_merged) final_pass=true; // If we are done, do one last pass for final redistributes

        //break;
    }

    // We should remap the labels to occupy the first natural numbers
    int labels_map[Kmax];
    for (int i=0; i<Kmax; i++)
        labels_map[i]=0;
    for (int i=0; i<(int)active_labels.size(); i++) {
        labels_map[active_labels[i]-1]=i+1;
    }
    for (long i=0; i<N; i++) {
        labels[i]=labels_map[labels[i]-1];
    }

    // If the user wants to refine the clusters, then we repeat isosplit on each
    // of the new clusters, recursively. Unless we only found only one cluster.
    int K=ns_isosplit5::compute_max(N,labels);

    if ((opts.refine_clusters)&&(K>1)) {
        int *labels_split=(int*)malloc(sizeof(int)*N);
        isosplit5_opts opts2=opts;
        opts2.refine_clusters=true; // Maybe we should provide an option on whether to do recursive refinement
        int k_offset=0;
        for (int k=1; k<=K; k++) {
            std::vector<long> inds_k;
            for (long i=0; i<N; i++)
                if (labels[i]==k)
                    inds_k.push_back(i);
            if (inds_k.size()>0) {
                float *X_k=(float*)malloc(sizeof(float)*M*inds_k.size()); //Warning: this may cause memory problems -- especially for recursive case
                int *labels_k=(int *)malloc(sizeof(int)*inds_k.size());
                for (long i=0; i<(long)inds_k.size(); i++) {
                    for (int m=0; m<M; m++) {
                        X_k[m+M*i]=X[m+M*inds_k[i]];
                    }
                }
                isosplit5(labels_k,M,inds_k.size(),X_k,opts2);
                for (long i=0; i<(long)inds_k.size(); i++) {
                    labels_split[inds_k[i]]=k_offset+labels_k[i];
                }
                k_offset+=ns_isosplit5::compute_max(inds_k.size(),labels_k);
                free(labels_k);
                free(X_k);
            }
        }
        for (long i=0; i<N; i++)
            labels[i]=labels_split[i];
        free(labels_split);
    }

    free(centroids);
}

/*



*/

/*
void isosplit5_old(int *labels_out,int M, long N,float *X,isosplit5_opts opts) {
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
*/

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

    long compute_max(long N,long *inds) {
        if (N==0) return 0;
        long ret=inds[0];
        for (long i=0; i<N; i++) {
            if (inds[i]>ret)
                ret=inds[i];
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

    void compute_centroids(float *centroids,int M,long N,int Kmax,float *X,int *labels) {
        std::vector<double> C(M*Kmax);
        std::vector<double> counts(Kmax);
        for (long i=0; i<N; i++) {
            int i0=labels[i]-1;
            for (int m=0; m<M; m++) {
                C[m+M*i0]+=X[m+M*i];
            }
            counts[i0]++;
        }
        for (int k=0; k<Kmax; k++) {
            if (counts[k]) {
                for (int m=0; m<M; m++) {
                    C[m+M*k]/=counts[k];
                }
            }
        }
        for (int jj=0; jj<M*Kmax; jj++)
            centroids[jj]=C[jj];
    }

    void get_pairs_to_compare(std::vector<int> *inds1,std::vector<int> *inds2,int M,int K,float *active_centroids,const intarray2d &active_comparisons_made) {
        inds1->clear();
        inds2->clear();
        double dists[K][K];
        for (int k1=0; k1<K; k1++) {
            for (int k2=0; k2<K; k2++) {
                if ((active_comparisons_made[k1][k2])||(k1==k2))
                    dists[k1][k2]=-1;
                else {
                    double dist=0;
                    for (int m=0; m<M; m++) {
                        double val=active_centroids[m+M*k1]-active_centroids[m+M*k2];
                        dist+=val*val;
                    }
                    dist=sqrt(dist);
                    dists[k1][k2]=dist;
                }
            }
        }
        bool something_changed=true;
        while (something_changed) {
            something_changed=false;
            std::vector<int> best_inds(K);
            for (int k=0; k<K; k++) {
                int best_ind=-1;
                double best_distance=-1;
                for (int k2=0; k2<K; k2++) {
                    if (dists[k][k2]>=0) {
                        if ((best_distance<0)||(dists[k][k2]<best_distance)) {
                            best_distance=dists[k][k2];
                            best_ind=k2;
                        }
                    }
                }
                best_inds[k]=best_ind;
            }
            for (int j=0; j<K; j++) {
                if (best_inds[j]>j) {
                    if (best_inds[best_inds[j]]==j) { //mutual
                        if (dists[j][best_inds[j]]>=0) {
                            inds1->push_back(j+1);
                            inds2->push_back(best_inds[j]+1);
                            for (int aa=0; aa<K; aa++) {
                                dists[j][aa]=-1;
                                dists[aa][j]=-1;
                                dists[best_inds[j]][aa]=-1;
                                dists[aa][best_inds[j]]=-1;
                            }
                            something_changed=true;
                        }
                    }
                }
            }
        }
    }

    std::vector<float> compute_centroid(int M,long N,float *X) {
        std::vector<double> ret(M);
        double count=0;
        for (int m=0; m<M; m++) {
            ret[m]=0;
        }
        for (long i=0; i<N; i++) {
            for (int m=0; m<M; m++) {
                ret[m]+=X[m+M*i];
            }
            count++;
        }
        if (count) {
            for (int m=0; m<M; m++) {
                ret[m]/=count;
            }
        }
        std::vector<float> retf(M);
        for (int m=0; m<M; m++)
            retf[m]=ret[m];
        return retf;
    }

    bool merge_test(std::vector<int> *L12,int M,long N1,long N2,float *X1,float *X2,const isosplit5_opts &opts) {
        L12->resize(N1+N2);
        for (long i=0; i<N1+N2; i++)
            (*L12)[i]=1;
        if ((N1==0)||(N2==0)) {
            printf("Error in merge test: N1 or N2 is zero.\n");
            return true;
        }
        std::vector<float> centroid1=compute_centroid(M,N1,X1);
        std::vector<float> centroid2=compute_centroid(M,N2,X2);
        std::vector<float> V(M);
        double sumsqr=0;
        for (int m=0; m<M; m++) {
            V[m]=centroid2[m]-centroid1[m];
            sumsqr+=V[m]*V[m];
        }
        if (sumsqr) {
            for (int m=0; m<M; m++)
                V[m]/=sqrt(sumsqr);
        }

        std::vector<float> projection1(N1),projection2(N2),projection12(N1+N2);
        for (long i=0; i<N1; i++) {
            double tmp=0;
            for (int m=0; m<M; m++)
                tmp+=V[m]*X1[m+i*M];
            projection1[i]=tmp;
            projection12[i]=tmp;
        }
        for (long i=0; i<N2; i++) {
            double tmp=0;
            for (int m=0; m<M; m++)
                tmp+=V[m]*X2[m+i*M];
            projection2[i]=tmp;
            projection12[N1+i]=tmp;
        }

        bool do_merge;
        isocut5_opts oo;
        oo.already_sorted=false;
        double dipscore,cutpoint;
        isocut5(&dipscore,&cutpoint,N1+N2,projection12.data(),oo);
        if (dipscore<opts.isocut_threshold) {
            do_merge=true;
        }
        else {
            do_merge=false;
        }
        for (long i=0; i<N1+N2; i++) {
            if (projection12[i]<cutpoint)
                (*L12)[i]=1;
            else
                (*L12)[i]=2;
        }

        return do_merge;
    }

    void compare_pairs(std::vector<int> *clusters_changed,long *total_num_label_changes,int M,long N,float *X,int *labels,const std::vector<int> &k1s,const std::vector<int> &k2s,const isosplit5_opts &opts) {
        printf("compare_pairs (M=%d, N=%ld, num pairs = %ld...\n",M,N,k1s.size());
        int Kmax=ns_isosplit5::compute_max(N,labels);
        int clusters_changed_vec[Kmax];
        for (int i=0; i<Kmax; i++)
            clusters_changed_vec[i]=0;
        int *new_labels=(int *)malloc(sizeof(int)*N);
        *total_num_label_changes=0;
        for (long i=0; i<N; i++)
            new_labels[i]=labels[i];
        for (int i1=0; i1<(int)k1s.size(); i1++) {
            int k1=k1s[i1];
            int k2=k2s[i1];
            printf("k1/k2 = %d/%d\n",k1,k2);
            std::vector<long> inds1,inds2;
            for (long i=0; i<N; i++) {
                if (labels[i]==k1)
                    inds1.push_back(i);
                if (labels[i]==k2)
                    inds2.push_back(i);
            }
            if ((inds1.size()>0)&&(inds2.size()>0)) {
                std::vector<long> inds12;
                inds12.insert(inds12.end(),inds1.begin(),inds1.end());
                inds12.insert(inds12.end(),inds2.begin(),inds2.end());
                std::vector<int> L12_old(inds12.size());
                for (long i=0; i<(long)inds1.size(); i++)
                    L12_old[i]=1;
                for (long i=0; i<(long)inds2.size(); i++)
                    L12_old[inds1.size()+i]=2;
                std::vector<int> L12(inds12.size());

                bool do_merge;
                if (((long)inds1.size()<opts.min_cluster_size)||((long)inds2.size()<opts.min_cluster_size)) {
                    do_merge=true;
                }
                else {
                    float *X1=(float*)malloc(sizeof(float)*M*inds1.size());
                    float *X2=(float*)malloc(sizeof(float)*M*inds2.size());
                    extract_subarray(X1,M,X,inds1);
                    extract_subarray(X2,M,X,inds2);
                    do_merge=merge_test(&L12,M,inds1.size(),inds2.size(),X1,X2,opts);
                    free(X1);
                    free(X2);
                }
                if (do_merge) {
                    for (long i=0; i<(long)inds2.size(); i++) {
                        new_labels[inds2[i]]=k1;
                    }
                    *total_num_label_changes+=inds2.size();
                    clusters_changed_vec[k1-1]=1;
                    clusters_changed_vec[k2-1]=1;
                }
                else {
                    //redistribute
                    bool something_was_redistributed=false;
                    for (long i=0; i<(long)inds1.size(); i++) {
                        if (L12[i]==2) {
                            new_labels[inds1[i]]=k2;
                            (*total_num_label_changes)++;
                            something_was_redistributed=true;
                        }
                    }
                    for (long i=0; i<(long)inds2.size(); i++) {
                        if (L12[inds1.size()+i]==1) {
                            new_labels[inds2[i]]=k1;
                            (*total_num_label_changes)++;
                            something_was_redistributed=true;
                        }
                    }
                    if (something_was_redistributed) {
                        clusters_changed_vec[k1-1]=1;
                        clusters_changed_vec[k2-1]=1;
                    }
                }
            }
        }
        clusters_changed->clear();
        for (int k=0; k<Kmax; k++)
            if (clusters_changed_vec[k])
                clusters_changed->push_back(k+1);
        for (long i=0; i<N; i++)
            labels[i]=new_labels[i];
        free(new_labels);
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
    ns_isosplit5::kmeans_multistep(labels,M,N,X,1000,0,0,oo);
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



