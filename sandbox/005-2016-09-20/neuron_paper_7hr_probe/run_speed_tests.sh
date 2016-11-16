#!/bin/bash

echo "10min $(date)"
kron-run ms2_10min np7 --_force_run > out_ms2_10min.txt
echo "30min $(date)"
kron-run ms2_30min np7 --_force_run > out_ms2_30min.txt
echo "60min $(date)"
kron-run ms2_60min np7 --_force_run > out_ms2_60min.txt
echo "120min $(date)"
kron-run ms2_120min np7 --_force_run > out_ms2_120min.txt
echo "240min $(date)"
kron-run ms2_240min np7 --_force_run > out_ms2_240min.txt
echo "420min $(date)"
kron-run ms2 np7 --_force_run > out_ms2.txt

echo "ch4 $(date)"
kron-run ms2_120min_ch4 np7 --_force_run > out_ms2_120min_ch4.txt
echo "ch8 $(date)"
kron-run ms2_120min_ch8 np7 --_force_run > out_ms2_120min_ch8.txt
echo "ch12 $(date)"
kron-run ms2_120min_ch12 np7 --_force_run > out_ms2_120min_ch12.txt

echo thr32 "$(date)$
kron-run ms2_120min_thr32 np7 --_force_run > out_ms2_120min_thr32.txt
echo thr16 "$(date)$
kron-run ms2_120min_thr16 np7 --_force_run > out_ms2_120min_thr16.txt
echo thr8 "$(date)$
kron-run ms2_120min_thr8 np7 --_force_run > out_ms2_120min_thr8.txt
echo thr4 "$(date)$
kron-run ms2_120min_thr4 np7 --_force_run > out_ms2_120min_thr4.txt
echo thr2 "$(date)$
kron-run ms2_120min_thr2 np7 --_force_run > out_ms2_120min_thr2.txt

