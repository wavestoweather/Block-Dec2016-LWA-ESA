# Number of samples for the Monte-Carlo confidence estimation
NSAMPLES := 1000

# How many ensemble members in the sampled ensemble?
SAMPLE_SIZE := 100

# Threshold for the correlation distribution width below which datapoints are
# considered to be statistically significant
SIGNIFICANCE_LEVEL := 0.09

# How many forecast members per cluster?
NCLUSTER := 15

# Determine filename suffix for C-extensions to installed Python
PYEXT := $(shell python3-config --extension-suffix)

# Names of all figures to be plotted
FIGURES := \
		figures/forecast-combination.pdf \
		figures/event24-reanalysis+nh18fig5.pdf \
		figures/event24-plume.pdf \
		figures/event24-budget.pdf \
		figures/event24-evaluation.pdf \
		figures/idealized.pdf \
		figures/event24-maps+hovmoeller.pdf \
		figures/event24-cluster+scatter.pdf \
		figures/event24-maps-separate.pdf \
		figures/event24-cluster-f2.pdf \
		figures/event24-nh18fig4.pdf


.PHONY: all clean

all: $(FIGURES)

clean:
	rm -f $(FIGURES)
	rm -f scripts/common/extensions/*.o
	rm -f scripts/common/_*.c
	rm -f scripts/common/_*.o
	rm -f scripts/common/_*.so
	rm -f scripts/hn2016_falwa/*.so


# Simulations with the traffic jam model

SIMULATIONS := \
	data/idealized-0.40-60-2.0.nc \
	data/idealized-0.40-55-2.0.nc \
	data/idealized-0.40-65-2.0.nc \
	data/idealized-0.40-60-1.8.nc \
	data/idealized-0.40-60-2.2.nc

figures/idealized.pdf: $(SIMULATIONS) scripts/plot-idealized.py
	python3 -m scripts.plot-idealized $(SIMULATIONS) --output $@ --data-output "data/idealized.json"

data/idealized-%-60-2.0.nc: scripts/run-idealized.py
	python3 -m scripts.run-idealized $@ --alpha="$*" --urefcg="60" --forcing-peak="-2.0"

data/idealized-0.40-%-2.0.nc: scripts/run-idealized.py
	python3 -m scripts.run-idealized $@ --alpha="0.40" --urefcg="$*" --forcing-peak="-2.0"

data/idealized-0.40-60-%.nc: scripts/run-idealized.py
	python3 -m scripts.run-idealized $@ --alpha="0.40" --urefcg="60" --forcing-peak="-$*"


# Python scipts for data processing and plotting, C-extensions

scripts/calculate-lwa.py: scripts/common/regrid.py scripts/hn2016_falwa/oopinterface.py
	touch $@

scripts/plot.py: scripts/common/plotting.py scripts/common/esa.py scripts/common/metrics.py scripts/common/texify.py
	touch $@

scripts/plot-budget.py: scripts/plot.py scripts/common/plotting.py scripts/common/texify.py
	touch $@

scripts/plot-evaluation.py: scripts/common/metrics.py scripts/common/plotting.py
	touch $@

scripts/common/metrics.py: scripts/common/texify.py
	touch $@

scripts/common/esa.py: scripts/common/_esa$(PYEXT)
	touch $@

scripts/common/regrid.py: scripts/common/_regrid$(PYEXT)
	touch $@

scripts/common/_%$(PYEXT): scripts/common/extensions/build_%.py scripts/common/extensions/%.c
	python3 $<


# hn2016_falwa package

scripts/hn2016_falwa/oopinterface.py: \
		scripts/hn2016_falwa/constant.py \
		scripts/hn2016_falwa/utilities.py \
		scripts/hn2016_falwa/interpolate_fields$(PYEXT) \
		scripts/hn2016_falwa/compute_lwa_and_barotropic_fluxes$(PYEXT) \
		scripts/hn2016_falwa/compute_reference_states$(PYEXT)
	touch $@

scripts/hn2016_falwa/%$(PYEXT): scripts/hn2016_falwa/%.f90
	f2py -c -m $* $<
	mv $(notdir $@) scripts/hn2016_falwa



# Data Processing (Part 1):
# Calculating LWA and related quantities with the hn2016_falwa package.

# Write processed data to netcdf files:
event24_ana := data/ERA-2016-12-10-to-2016-12-24-lwa-qg.nc
event24_ens := data/ENS-2016-12-09T12Z-lwa-qg.nc data/ENS-2016-12-10T00Z-lwa-qg.nc data/ENS-2016-12-10T12Z-lwa-qg.nc data/ENS-2016-12-11T00Z-lwa-qg.nc
event24_eval := data/EVAL-2016-12-18T00Z-lwa-qg.nc

# Don't remove these intermediate files
.SECONDARY: $(event24_ana) $(event24_ens) $(event24_eval)
.SECONDARY: data/ERA-2016-12-10-to-2016-12-24-uvt.nc

# Reanalysis data
data/ERA-%-lwa-qg.nc: data/ERA-%-uvt.nc scripts/calculate-lwa.py
	python3 -m scripts.calculate-lwa --budget --half-resolution $<

# Ensemble data
data/ENS-%-lwa-qg.nc: scripts/calculate-lwa.py data/IFS-ENS/ENS-%/ENS-*-[uvt].nc
	python3 -m scripts.calculate-lwa --half-resolution $(filter-out $<,$^)

# Forecast evaluation data for 18 Dec 00 UTC
data/EVAL-2016-12-18T00Z-lwa-qg.nc: scripts/calculate-lwa.py data/IFS-ENS/EVAL/ENS-DEC18-EVAL-*.nc
	python3 -m scripts.calculate-lwa --half-resolution --eval data/IFS-ENS/EVAL/ENS-DEC18-EVAL-*.nc


# Data Processing (Part 2):
# Analysis and plotting of figures.

# Figure 1: Ensemble combination
figures/forecast-combination.pdf: figures/forecast-combination/forecast-combination.tex
	cd $(dir $<) && lualatex $(notdir $<)
	mv $(<:.tex=.pdf) $@


# Figure 2: Reanalysis overview
figures/event24-reanalysis+nh18fig5.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot reanalysis+nh18fig5 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--hovmoeller-extent "target" \
			--end 2016-12-24T00


# Figure 3: Ensemble target metric overview
# Figure 7: ESA maps and Hovmöller with full flux
# Figure 8: Cluster and scatter analysis of upstream region (needs data from idealized experiment)
figures/event24-%.pdf: scripts/plot.py figures/idealized.pdf $(event24_ana) $(event24_ens)
	python3 -m scripts.plot $* 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--ncluster $(NCLUSTER) \
			--delta-hours "0,-24,-48,-72,-96" \
			--source "f1+f2+f3" \
			--scatter 2016-12-16T00,65,-90,35,-15 \
			--scatter-idealized "data/idealized.json" \
			--hovmoeller-extent "target" \
			--end 2016-12-21T00


# Figure 4: Integrated budget in 2.5 days prior to onset
figures/event24-budget.pdf: scripts/plot-budget.py $(event24_ana)
	python3 -m scripts.plot-budget $(event24_ana) 2016-12-15T12 2016-12-18T00,65,10,35,-30 \
			--statistics-output "data/budget-statistics.json" \
			--output $@


# Figure 5: Target forecast evaluation
figures/event24-evaluation.pdf: scripts/plot-evaluation.py $(event24_eval) $(event24_ana)
	python3 -m scripts.plot-evaluation 65,10,35,-30 $(event24_eval) \
		--reanalysis $(event24_ana) \
		--highlight 2016-12-09T12 2016-12-10T00 2016-12-10T12 2016-12-11T00 \
		--output $@


# Figure 9: ESA maps with flux terms
figures/event24-maps-separate.pdf: figures/event24-maps-separate/event24-maps-separate.tex \
		figures/event24-maps-separate/f1f3.pdf figures/event24-maps-separate/f2.pdf
	cd $(dir $<) && lualatex $(notdir $<)
	mv $(<:.tex=.pdf) $@

figures/event24-maps-separate/f1f3.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot maps4 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--ncluster $(NCLUSTER) \
			--delta-hours "24, 0,-24, -48" \
			--source "f1+f3" \
			--panel-letters "abcd"

figures/event24-maps-separate/f2.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot maps4 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--ncluster $(NCLUSTER) \
			--delta-hours "24, 0,-24, -48" \
			--source "f2" \
			--panel-letters "efgh"


# Figure 10: Cluster analysis of F2 in block
figures/event24-cluster-f2.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot cluster 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--ncluster $(NCLUSTER) \
			--source "f2"


# Figure 11: Flux-LWA relationship in block (similar to NH18's Fig. 4)
figures/event24-nh18fig4.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot nh18fig4 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--ncluster $(NCLUSTER) \
			--scatter 2016-12-18T00,50,-5,40,-15
