firings=${1:-msfirings1.accepted.mda.prv}

mountainview --samplerate=30000 --raw=raw.mda.prv --pre=output/pre2.mda.prv --filt=output/pre1b.mda.prv --firings=$firings
