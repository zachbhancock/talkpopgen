#This code is for simulating populations using msprime and tskit and is from the video "The Causes of Evolution: Variation". For this code to run, you will need to have downloaded msprime, tskit, numpy, and pyslim. Much of this code is from the tskit and msprime tutorial websites, where you can find more details to run more complicated simulations and the general API.

import numpy as np
import msprime, pyslim
import tskit
from IPython.display import SVG, display

#simulate samples
ts = msprime.sim_ancestry(
 samples=3,
 recombination_rate=1e-8,
 sequence_length=5_000,
 population_size=10_000,
 random_seed=123456)

# Visualise the simulated ancestral history.
SVG(ts.draw_svg())

#add mutations to the simulated sample set
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=54321)
SVG(ts.draw_svg())

#view each of the variants and compare them to the ancestral state
for variant in ts.variants():
    print(variant)
    
#estimate nucleotide diversity
d = ts.diversity()

#estimate nucleotide diversity for discrete windows
windows = np.linspace(0, ts.sequence_length, num=3)
d = ts.diversity(windows=windows)
print(windows, d, sep="\n")

#compute the site frequency spectrum
afs = ts.allele_frequency_spectrum(polarised=True, span_normalise=False)
print(afs)