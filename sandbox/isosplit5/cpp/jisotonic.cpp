#include "jisotonic.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

void jisotonic(int N, float* BB, float* MSE, float* AA, float* WW)
{
    if (N < 1)
        return;

    float* unweightedcount = (float*)malloc(sizeof(float) * N);
    float* count = (float*)malloc(sizeof(float) * N);
    float* sum = (float*)malloc(sizeof(float) * N);
    float* sumsqr = (float*)malloc(sizeof(float) * N);
    int last_index = -1;

    last_index++;
    unweightedcount[last_index] = 1;
    float w0 = 1;
    if (WW)
        w0 = WW[0];
    count[last_index] = w0;
    sum[last_index] = AA[0] * w0;
    sumsqr[last_index] = AA[0] * AA[0] * w0;
    MSE[0] = 0;

    for (int j = 1; j < N; j++) {
        last_index++;
        unweightedcount[last_index] = 1;
        w0 = 1;
        if (WW)
            w0 = WW[j];
        count[last_index] = w0;
        sum[last_index] = AA[j] * w0;
        sumsqr[last_index] = AA[j] * AA[j] * w0;
        MSE[j] = MSE[j - 1];
        while (true) {
            if (last_index <= 0)
                break;
            if (sum[last_index - 1] / count[last_index - 1] < sum[last_index] / count[last_index]) {
                break;
            }
            else {
                float prevMSE = sumsqr[last_index - 1] - sum[last_index - 1] * sum[last_index - 1] / count[last_index - 1];
                prevMSE += sumsqr[last_index] - sum[last_index] * sum[last_index] / count[last_index];
                unweightedcount[last_index - 1] += unweightedcount[last_index];
                count[last_index - 1] += count[last_index];
                sum[last_index - 1] += sum[last_index];
                sumsqr[last_index - 1] += sumsqr[last_index];
                float newMSE = sumsqr[last_index - 1] - sum[last_index - 1] * sum[last_index - 1] / count[last_index - 1];
                MSE[j] += newMSE - prevMSE;
                last_index--;
            }
        }
    }

    int ii = 0;
    for (int k = 0; k <= last_index; k++) {
        for (int cc = 0; cc < unweightedcount[k]; cc++) {
            BB[ii + cc] = sum[k] / count[k];
        }
        ii += unweightedcount[k];
    }

    free(unweightedcount);
    free(count);
    free(sum);
    free(sumsqr);
}

void jisotonic_updown(int N, float* out, float* in, float* weights)
{
    float* B1 = (float*)malloc(sizeof(float) * N);
    float* MSE1 = (float*)malloc(sizeof(float) * N);
    float* B2 = (float*)malloc(sizeof(float) * N);
    float* MSE2 = (float*)malloc(sizeof(float) * N);
    float* in_reversed = (float*)malloc(sizeof(float) * N);
    float* weights_reversed = 0;

    for (int j = 0; j < N; j++) {
        in_reversed[j] = in[N - 1 - j];
    }
    if (weights) {
        weights_reversed = (float*)malloc(sizeof(float) * N);
        for (int j = 0; j < N; j++) {
            weights_reversed[j] = weights[N - 1 - j];
        }
    }
    jisotonic(N, B1, MSE1, in, weights);
    jisotonic(N, B2, MSE2, in_reversed, weights_reversed);
    for (int j = 0; j < N; j++)
        MSE1[j] += MSE2[N - 1 - j];
    float bestval = MSE1[0];
    int best_ind = 0;
    for (int j = 0; j < N; j++) {
        if (MSE1[j] < bestval) {
            bestval = MSE1[j];
            best_ind = j;
        }
    }
    jisotonic(best_ind + 1, B1, MSE1, in, weights);
    jisotonic(N - best_ind, B2, MSE2, in_reversed, weights_reversed);
    for (int j = 0; j <= best_ind; j++)
        out[j] = B1[j];
    for (int j = 0; j < N - best_ind - 1; j++)
        out[N - 1 - j] = B2[j];

    free(B1);
    free(MSE1);
    free(B2);
    free(MSE2);
    free(in_reversed);
    if (weights_reversed)
        free(weights_reversed);
}

void jisotonic_downup(int N, float* out, float* in, float* weights)
{
    float* in_neg = (float*)malloc(sizeof(float) * N);

    for (int j = 0; j < N; j++)
        in_neg[j] = -in[j];
    jisotonic_updown(N, out, in_neg, weights);
    for (int j = 0; j < N; j++)
        out[j] = -out[j];

    free(in_neg);
}

void jisotonic_sort(int N, float* out, const float* in)
{
    std::copy(in, in + N, out);
    std::sort(out, out + N);

    //	QVector<float> in0(N);
    //	for (int j=0; j<N; j++) in0[j]=in[j];
    //	qSort(in0);
    //	for (int j=0; j<N; j++) out[j]=in0[j];
}
