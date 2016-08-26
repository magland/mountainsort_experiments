firings1=${1:-manclust1_firings.mda.prv}
firings2=${2:-msfirings1.accepted.mda.prv}

mountaincompare --samplerate=30000 --raw=raw.mda.prv --pre=output/pre2.mda.prv --filt=output/pre1b.mda.prv --firings1=$firings1 --firings2=$firings2
