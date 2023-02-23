.PHONY: all thesis_images paper_images

THESISIMAGES = 	 \
			 images/thesis/ch2-diff-theta-histogram.png \
			 images/thesis/ch2-pearson-correlation.png \
			 images/thesis/ch4-extinction-frequencies.png \
			 images/thesis/ch4-extinction-frequencies-vertex.png \
			 images/thesis/ch4-extinction-frequencies-extrusion.png \
			 images/thesis/ch4-effective-g1-durations.png \
			 images/thesis/ch4-effective-g1-duration-sample-sizes.png \
			 images/thesis/ch5-exponential-g1-proportion.png \
			 images/thesis/ch5-exponential-g1-proportion-well-mixed.png \
			 images/thesis/ch5-exponential-g1-proportion-vertex.png \
			 images/thesis/ch5-exponential-g1-proportion-vertex-no-extrusions.png \
			 images/thesis/ch5-homotypic-exponential-proliferation-regimes.png \
			 images/thesis/ch5-homotypic-exponential-estimated-proliferation-regimes.png \
			 images/thesis/ch5-quadratic-bounds.png \
			 images/thesis/ch5-uniform-parameter-space.png \
			 images/thesis/ch5-uniform-g1-proportion.png \
			 images/thesis/ch5-uniform-g1-proportion-well-mixed.png \
			 images/thesis/ch5-uniform-g1-proportion-vertex.png \
			 images/thesis/ch5-homotypic-uniform-proliferation-regimes.png \
			 images/thesis/ch5-homotypic-uniform-estimated-proliferation-regimes.png \
			 images/thesis/ch6-heterotypic-survival-difference-exponential-well-mixed.png \
			 images/thesis/ch6-heterotypic-survival-difference-exponential-vertex.png \
			 images/thesis/ch6-homotypic-survival-difference-exponential-A-well-mixed.png \
			 images/thesis/ch6-homotypic-survival-difference-exponential-B-well-mixed.png \
			 images/thesis/ch6-homotypic-survival-difference-exponential-A-vertex.png \
			 images/thesis/ch6-homotypic-survival-difference-exponential-B-vertex.png \
			 images/thesis/ch6-classification-of-competitive-interactions.png \
			 images/thesis/ch6-heterotypic-asymptotic-behaviour.png \
			 images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-well-mixed.png \
			 images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-vertex-random.png \
			 images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-vertex-segregated.png \
			 images/thesis/ch6-heterotypic-competition-regimes.png \
			 images/thesis/ch6-heterotypic-competition-regimes-transformed.png \
			 images/thesis/ch7-constant-emission.png \

all: thesis_images paper_images

thesis_images: $(THESISIMAGES)

paper_images: $(PAPERIMAGES)

# create directories
$(shell mkdir -p images/thesis)
$(shell mkdir -p images/paper)

# Thesis image rules
images/thesis/ch2-diff-theta-histogram.png: \
	data/mechanical-heterotypic-tests-data-28-10-2021.csv \
	data/mechanical-homotypic-tests-data-28-10-2021.csv \
	figures/thesis/TT-18-2-analysis.py
	python3 figures/thesis/TT-18-2-analysis.py \
		data/mechanical-heterotypic-tests-data-28-10-2021.csv  \
		data/mechanical-homotypic-tests-data-28-10-2021.csv \
		diff-theta-histogram \
		$@

images/thesis/ch2-pearson-correlation.png: \
	data/mechanical-heterotypic-tests-data-28-10-2021.csv \
	data/mechanical-homotypic-tests-data-28-10-2021.csv \
	figures/thesis/TT-18-2-analysis.py
	python3 figures/thesis/TT-18-2-analysis.py \
		data/mechanical-heterotypic-tests-data-28-10-2021.csv  \
		data/mechanical-homotypic-tests-data-28-10-2021.csv \
		pearson-correlation \
		$@

images/thesis/ch4-extinction-frequencies.png: \
	figures/thesis/formatting.py \
	data/mc-extinction-exponential-data.csv \
	figures/thesis/plot-extinction-frequencies.py
	python3 figures/thesis/plot-extinction-frequencies.py \
		data/mc-extinction-exponential-data.csv \
		well-mixed \
		$@

images/thesis/ch4-extinction-frequencies-vertex.png: \
	figures/thesis/formatting.py \
	data/extinction-probability-regime-3-tests-data-01-10-2021.csv \
	figures/thesis/plot-extinction-frequencies.py
	python3 figures/thesis/plot-extinction-frequencies.py \
		data/extinction-probability-regime-3-tests-data-01-10-2021.csv \
		vertex \
		$@

images/thesis/ch4-extinction-frequencies-extrusion.png: \
	figures/thesis/formatting.py \
	data/extinction-probability-regime-3-tests-data-01-10-2021.csv \
	figures/thesis/plot-extinction-frequencies.py
	python3 figures/thesis/plot-extinction-frequencies.py \
		data/extinction-probability-regime-3-tests-data-01-10-2021.csv \
		extrusion \
		$@

images/thesis/ch4-effective-g1-durations.png: \
	figures/thesis/formatting.py \
	data/mc-g1-truncation-exponential-5-data.csv \
	data/mc-g1-truncation-uniform-5-data.csv \
	figures/thesis/plot-effective-g1-durations.py
	python3 figures/thesis/plot-effective-g1-durations.py \
		data/mc-g1-truncation-exponential-5-data.csv \
		data/mc-g1-truncation-uniform-5-data.csv \
		$@

images/thesis/ch4-effective-g1-duration-sample-sizes.png: \
	figures/thesis/formatting.py \
	data/mc-g1-truncation-exponential-5-data.csv \
	data/mc-g1-truncation-uniform-5-data.csv \
	figures/thesis/plot-effective-g1-duration-sample-sizes.py
	python3 figures/thesis/plot-effective-g1-duration-sample-sizes.py \
		data/mc-g1-truncation-exponential-5-data.csv \
		data/mc-g1-truncation-uniform-5-data.csv \
		$@

images/thesis/ch5-exponential-g1-proportion.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-exponential-g1-proportion.py
	python3 figures/thesis/plot-exponential-g1-proportion.py \
		$@

images/thesis/ch5-exponential-g1-proportion-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-exponential-g1-proportion-well-mixed.py \
	data/homotypic-g2-death-signal-exponential-data-26-08-2021.csv
	python3 figures/thesis/plot-exponential-g1-proportion-well-mixed.py \
		data/homotypic-g2-death-signal-exponential-data-26-08-2021.csv \
		n \
		$@

images/thesis/ch5-exponential-g1-proportion-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-exponential-g1-proportion-vertex.py \
	data/g1-proportion-tests-data-07-10-2021.csv
	python3 figures/thesis/plot-exponential-g1-proportion-vertex.py \
		data/g1-proportion-tests-data-07-10-2021.csv \
		$@

images/thesis/ch5-exponential-g1-proportion-vertex-no-extrusions.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-exponential-g1-proportion-vertex.py \
	data/g1-proportion-tests-data-07-10-2021.csv
	python3 figures/thesis/plot-exponential-g1-proportion-vertex.py \
		data/g1-proportion-tests-data-07-10-2021.csv \
		$@ \
		y

images/thesis/ch5-homotypic-exponential-proliferation-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-exponential-asymptotic-behaviour.py
	python3 figures/thesis/plot-homotypic-exponential-asymptotic-behaviour.py \
		$@

images/thesis/ch5-homotypic-exponential-estimated-proliferation-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-asymptotic-exponential.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
	data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv
	python3 figures/thesis/plot-homotypic-asymptotic-exponential.py \
		data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
		data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv \
		$@

images/thesis/ch5-quadratic-bounds.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-quadratic-bounds.py
	python3 figures/thesis/plot-quadratic-bounds.py \
		$@

images/thesis/ch5-uniform-parameter-space.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-uniform-parameter-space.py
	python3 figures/thesis/plot-uniform-parameter-space.py \
		$@

images/thesis/ch5-uniform-g1-proportion.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-uniform-g1-proportion.py
	python3 figures/thesis/plot-uniform-g1-proportion.py \
		$@

images/thesis/ch5-uniform-g1-proportion-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-uniform-g1-proportion-well-mixed.py \
	data/homotypic-model-2-uniform-g1-proportion-data-06-05-2021.csv
	python3 figures/thesis/plot-uniform-g1-proportion-well-mixed.py \
		data/homotypic-model-2-uniform-g1-proportion-data-06-05-2021.csv \
		$@

images/thesis/ch5-uniform-g1-proportion-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-uniform-g1-proportion-vertex.py \
	data/g1-proportion-tests-data-07-10-2021.csv
	python3 figures/thesis/plot-uniform-g1-proportion-vertex.py \
		data/g1-proportion-tests-data-07-10-2021.csv \
		$@

images/thesis/ch5-homotypic-uniform-proliferation-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-uniform-asymptotic-behaviour.py
	python3 figures/thesis/plot-homotypic-uniform-asymptotic-behaviour.py \
		$@

images/thesis/ch5-homotypic-uniform-estimated-proliferation-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-asymptotic-uniform.py \
	figures/thesis/plot_survival_difference.py \
	data/homotypic-model-2-uniform-asymptotic-I-data-02-10-2021.csv \
	data/homotypic-model-2-uniform-asymptotic-II-data-02-10-2021.csv \
	data/homotypic-model-2-uniform-asymptotic-I-tests-data-04-10-2021.csv \
	data/homotypic-model-2-uniform-asymptotic-II-tests-data-04-10-2021.csv
	python3 figures/thesis/plot-homotypic-asymptotic-uniform.py \
		data/homotypic-model-2-uniform-asymptotic-I-data-02-10-2021.csv \
		data/homotypic-model-2-uniform-asymptotic-II-data-02-10-2021.csv \
		data/homotypic-model-2-uniform-asymptotic-I-tests-data-04-10-2021.csv \
		data/homotypic-model-2-uniform-asymptotic-II-tests-data-04-10-2021.csv \
		$@

images/thesis/ch6-heterotypic-survival-difference-exponential-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-survival-difference-exponential-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-data-27-07-2021.csv
	python3 figures/thesis/plot-heterotypic-survival-difference-exponential-well-mixed.py \
		data/heterotypic-model-2-exponential-data-27-07-2021.csv \
		$@

images/thesis/ch6-heterotypic-survival-difference-exponential-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-survival-difference-exponential-vertex.py \
	figures/thesis/plot_survival_difference.py data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv
	python3 figures/thesis/plot-heterotypic-survival-difference-exponential-vertex.py  \
		data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
		$@

images/thesis/ch6-homotypic-survival-difference-exponential-A-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-exponential-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-data-27-07-2021.csv \
	data/heterotypic-model-2-exponential-control-data-04-10-2021.csv
	python3 figures/thesis/plot-homotypic-survival-difference-exponential-well-mixed.py \
		data/heterotypic-model-2-exponential-data-27-07-2021.csv \
		data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
		A \
		$@

images/thesis/ch6-homotypic-survival-difference-exponential-B-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-exponential-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-data-27-07-2021.csv \
	data/heterotypic-model-2-exponential-control-data-04-10-2021.csv
	python3 figures/thesis/plot-homotypic-survival-difference-exponential-well-mixed.py \
		data/heterotypic-model-2-exponential-data-27-07-2021.csv \
		data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
		B \
		$@

images/thesis/ch6-homotypic-survival-difference-exponential-A-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-exponential-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
	data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv
	python3 figures/thesis/plot-homotypic-survival-difference-exponential-vertex.py \
		data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
		data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv \
		A \
		$@

images/thesis/ch6-homotypic-survival-difference-exponential-B-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-exponential-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
	data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv
	python3 figures/thesis/plot-homotypic-survival-difference-exponential-vertex.py \
		data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
		data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv \
		B \
		$@

images/thesis/ch6-classification-of-competitive-interactions.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-classification-of-competitive-interactions.py
	python3 figures/thesis/plot-classification-of-competitive-interactions.py \
		$@

images/thesis/ch6-heterotypic-asymptotic-behaviour.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-asymptotic-behaviour.py
	python3 figures/thesis/plot-heterotypic-asymptotic-behaviour.py \
		$@

images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-asymptotic-exponential-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-asymptotic-I-data-02-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-II-data-02-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-III-data-02-08-2021.csv
	python3 figures/thesis/plot-heterotypic-asymptotic-exponential-well-mixed.py \
		data/heterotypic-model-2-exponential-asymptotic-I-data-02-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-II-data-02-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-III-data-02-08-2021.csv \
		$@

images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-vertex-random.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-asymptotic-exponential-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv
	python3 figures/thesis/plot-heterotypic-asymptotic-exponential-vertex.py \
		random \
		data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv \
		$@

images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-vertex-segregated.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-asymptotic-exponential-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv
	python3 figures/thesis/plot-heterotypic-asymptotic-exponential-vertex.py \
		segregated \
		data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv \
		$@

images/thesis/ch6-heterotypic-competition-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-heterotypic-competition-regimes.py
	python3 figures/thesis/plot-heterotypic-competition-regimes.py \
		I-II \
		$@

images/thesis/ch6-heterotypic-competition-regimes-transformed.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-heterotypic-competition-regimes.py
	python3 figures/thesis/plot-competition-regimes-transformed.py \
		$@

images/thesis/ch7-constant-emission.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-constant-emission.py
	python3 figures/thesis/plot-constant-emission.py \
		$@

