/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef jisotonic_h
#define jisotonic_h

void jisotonic(int N, float* BB, float* MSE, float* AA, float* WW);
void jisotonic_updown(int N, float* out, float* in, float* weights);
void jisotonic_downup(int N, float* out, float* in, float* weights);
void jisotonic_sort(int N, float* out, const float* in);

#endif
