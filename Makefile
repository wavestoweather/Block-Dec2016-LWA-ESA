# Number of samples for the Monte-Carlo confidence estimation
NSAMPLES := 1000

# How many ensemble members in the sampled ensemble?
SAMPLE_SIZE := 75

# Threshold for the correlation distribution width below which datapoints are
# considered to be statistically significant
SIGNIFICANCE_LEVEL := 0.1

# Determine filename suffix for C-extensions to installed Python
PYEXT := $(shell python3-config --extension-suffix)

# Names of all figures to be plotted
FIGURES := \
		figures/forecast-combination.pdf \
		figures/event24-reanalysis+nh18fig5.pdf \
		figures/event24-plume.pdf \
		figures/event24-maps+hovmoeller.pdf \
		figures/event24-cluster+scatter.pdf \
		figures/event24-maps-separate.pdf \
		figures/event24-cluster-f2.pdf \
		figures/event24-nh18fig4.pdf


.PHONY: all clean local

all: $(FIGURES)

clean:
	rm -f $(FIGURES)
	rm -f scripts/common/extensions/*.o
	rm -f scripts/common/_*.c
	rm -f scripts/common/_*.o
	rm -f scripts/common/_*.so
	rm -f scripts/hn2016_falwa/*.so


# Python scipts for data processing and plotting, C-extensions

scripts/calculate-lwa.py: scripts/common/regrid.py scripts/hn2016_falwa/oopinterface.py
	touch $@

scripts/plot.py: scripts/common/plotting.py scripts/common/esa.py scripts/common/metrics.py scripts/common/texify.py
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
event24_ens := data/ENS-2016-12-10T00Z-lwa-qg.nc data/ENS-2016-12-10T12Z-lwa-qg.nc data/ENS-2016-12-11T00Z-lwa-qg.nc
# Don't remove these intermediate files
.SECONDARY: $(event24_ana) $(event24_ens)

# Reanalysis data
data/ERA-%-lwa-qg.nc: data/ERA-%-uvt.nc scripts/calculate-lwa.py
	python3 -m scripts.calculate-lwa --half-resolution $<

# Ensemble data
data/ENS-%-lwa-qg.nc: scripts/calculate-lwa.py data/IFS-ENS/ENS-%/ENS-*-[uvt].nc
	python3 -m scripts.calculate-lwa --half-resolution $(filter-out $<,$^)



# Data Processing (Part 2):
# Analysis and plotting of figures.

# Figure 1: Ensemble combination
figures/forecast-combination.pdf: figures/forecast-combination/forecast-combination.tex
	cd $(dir $<) && lualatex $(notdir $<)
	mv $(<:.tex=.pdf) $@

# Figure 2: Reanalysis overview
# Figure 3: Ensemble target metric overview
# Figure 4: ESA maps and Hovmöller with full flux
# Figure 5: Cluster and scatter analysis of upstream region
figures/event24-%.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot $* 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--delta-hours "0,-24,-48,-72,-96" \
			--source "f1+f2+f3" \
			--scatter 2016-12-16T00,65,-90,35,-15 \
			--hovmoeller-extent "target" \
			--end 2016-12-24T00

# Figure 6: ESA maps with flux terms
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
			--delta-hours "24, 0,-24, -48" \
			--source "f2" \
			--panel-letters "efgh"

# Figure 7: Cluster analysis of F2 in block
figures/event24-cluster-f2.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot cluster 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--source "f2"

# Figure 8: Flux-LWA relationship in block (similar to NH18's Fig. 4)
figures/event24-nh18fig4.pdf: scripts/plot.py $(event24_ana) $(event24_ens)
	python3 -m scripts.plot nh18fig4 2016-12-18T00,65,10,35,-30 $(event24_ens) \
			--reanalysis $(event24_ana) \
			--output $@ \
			--nsamples $(NSAMPLES) \
			--sample-size $(SAMPLE_SIZE) \
			--significance-level $(SIGNIFICANCE_LEVEL) \
			--scatter 2016-12-18T00,50,-5,40,-15

