.PHONY: all thesis_images paper_images

THESISIMAGES = 	 \
			 images/thesis/chapter-2/survival-frequency-difference-histograms.png \
			 images/thesis/chapter-2/pearson-correlation.png \
			 images/thesis/chapter-4/extinction-frequencies-well-mixed.png \
			 images/thesis/chapter-4/extinction-frequencies-vertex.png \
			 images/thesis/chapter-4/effective-g1-durations.png \
			 images/thesis/chapter-5/homotypic-survival-probability.png \
			 images/thesis/chapter-5/exponential-simulations.png \
			 images/thesis/chapter-5/homotypic-survival-frequency-well-mixed.png \
			 images/thesis/chapter-5/homotypic-survival-frequency-vertex.png \
			 images/thesis/chapter-5/homotypic-survival-frequency-vertex-no-extrusions.png \
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
			 images/thesis/chapter-6/classification-of-competitive-interactions.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes.png \
			 images/thesis/ch6-heterotypic-asymptotic-survival-frequency-exponential-well-mixed.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes-vertex-random.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes-vertex-segregated.png \
			 images/thesis/ch6-heterotypic-competition-regimes.png \
			 images/thesis/chapter-6/competition-regimes-transformed.png \
			 images/thesis/chapter-7/constant-emission.png \

PAPERIMAGES = \
		 images/paper/homotypic-proliferation-regimes-diagram.png \
		 images/paper/homotypic-proliferation-regimes-estimated.png \
		 images/paper/heterotypic-proliferation-regimes-diagram.png \
		 images/paper/heterotypic-proliferation-regimes-estimated-well-mixed.png \
		 images/paper/heterotypic-proliferation-regimes-estimated-vertex-random.png \
		 images/paper/heterotypic-proliferation-regimes-estimated-vertex-segregated.png \
		 images/paper/classification-of-competitive-interactions.png \
		 images/paper/competition-regimes.png \
		 images/paper/competition-regimes-transformed.png \
		 images/paper/heterotypic-survival-difference-well-mixed.png \
		 images/paper/heterotypic-survival-difference-vertex.png \
		 images/paper/homotypic-survival-difference-A-well-mixed.png \
		 images/paper/homotypic-survival-difference-B-well-mixed.png \
		 images/paper/homotypic-survival-difference-A-vertex.png \
		 images/paper/homotypic-survival-difference-B-vertex.png \
		 images/paper/homotypic-survival-frequency-well-mixed.png \
		 images/paper/homotypic-survival-frequency-vertex.png

all: thesis_images paper_images

thesis_images: $(THESISIMAGES)

paper_images: $(PAPERIMAGES)

# create directories
$(shell mkdir -p images/thesis/chapter-2)
$(shell mkdir -p images/thesis/chapter-4)
$(shell mkdir -p images/thesis/chapter-5)
$(shell mkdir -p images/thesis/chapter-6)
$(shell mkdir -p images/thesis/chapter-7)
$(shell mkdir -p images/paper)

# Thesis image rules
images/thesis/chapter-2/survival-frequency-difference-histograms.png: \
	data/mechanical-cell-competition-heterotypic-data.csv \
	data/mechanical-cell-competition-homotypic-data.csv \
	figures/thesis/mechanical-cell-competition-analysis.py
	python3 figures/thesis/mechanical-cell-competition-analysis.py \
		data/mechanical-cell-competition-heterotypic-data.csv  \
		data/mechanical-cell-competition-homotypic-data.csv \
		diff-theta-histogram \
		$@

images/thesis/chapter-2/pearson-correlation.png: \
	data/mechanical-cell-competition-heterotypic-data.csv \
	data/mechanical-cell-competition-homotypic-data.csv \
	figures/thesis/mechanical-cell-competition-analysis.py
	python3 figures/thesis/mechanical-cell-competition-analysis.py \
		data/mechanical-cell-competition-heterotypic-data.csv  \
		data/mechanical-cell-competition-homotypic-data.csv \
		pearson-correlation \
		$@

images/thesis/chapter-4/extinction-frequencies-well-mixed.png: \
	figures/thesis/formatting.py \
	data/extinction-frequencies-well-mixed-data.csv \
	figures/thesis/plot-extinction-frequencies.py
	python3 figures/thesis/plot-extinction-frequencies.py \
		data/extinction-frequencies-well-mixed-data.csv\
		well-mixed \
		$@

images/thesis/chapter-4/extinction-frequencies-vertex.png: \
	figures/thesis/formatting.py \
	data/extinction-frequencies-vertex-data.csv \
	figures/thesis/plot-extinction-frequencies.py
	python3 figures/thesis/plot-extinction-frequencies.py \
		data/extinction-frequencies-vertex-data.csv \
		extrusion \
		$@

images/thesis/chapter-4/effective-g1-durations.png: \
	figures/thesis/formatting.py \
	data/effective-g1-durations-exponential-data.csv \
	data/effective-g1-durations-uniform-data.csv \
	figures/thesis/plot-effective-g1-durations.py
	python3 figures/thesis/plot-effective-g1-durations.py \
		data/effective-g1-durations-exponential-data.csv \
		data/effective-g1-durations-uniform-data.csv \
		$@

images/thesis/chapter-5/homotypic-survival-probability.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-probability.py
	python3 figures/thesis/plot-homotypic-survival-probability.py \
		$@

images/thesis/chapter-5/exponential-simulations.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-exponential-simulations.py \
	data/eta-0.2-simulation.csv \
	data/eta-0.05-simulation.csv
	python3 figures/thesis/plot-exponential-simulations.py \
		data/eta-0.2-simulation.csv \
		data/eta-0.05-simulation.csv \
		$@

images/thesis/chapter-5/homotypic-survival-frequency-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-frequency-well-mixed.py \
	data/homotypic-survival-frequency-well-mixed-data.csv
	python3 figures/thesis/plot-homotypic-survival-frequency-well-mixed.py \
		data/homotypic-survival-frequency-well-mixed-data.csv \
		n \
		$@

images/thesis/chapter-5/homotypic-survival-frequency-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-frequency-vertex.py \
	data/homotypic-survival-frequency-vertex-data.csv
	python3 figures/thesis/plot-homotypic-survival-frequency-vertex.py \
		data/homotypic-survival-frequency-vertex-data.csv \
		$@

images/thesis/chapter-5/homotypic-survival-frequency-vertex-no-extrusions.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-frequency-vertex.py \
	data/homotypic-survival-frequency-vertex-data.csv
	python3 figures/thesis/plot-homotypic-survival-frequency-vertex.py \
		data/homotypic-survival-frequency-vertex-data.csv \
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

images/thesis/chapter-6/classification-of-competitive-interactions.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-classification-of-competitive-interactions.py
	python3 figures/thesis/plot-classification-of-competitive-interactions.py \
		$@

images/thesis/chapter-6/heterotypic-proliferation-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-proliferation-regimes.py
	python3 figures/thesis/plot-heterotypic-proliferation-regimes.py \
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

images/thesis/chapter-6/heterotypic-proliferation-regimes-vertex-random.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-proliferation-regimes-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv
	python3 figures/thesis/plot-heterotypic-proliferation-regimes-vertex.py \
		random \
		data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv \
		$@

images/thesis/heterotypic-proliferation-regimes-vertex-segregated.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-proliferation-regimes-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv
	python3 figures/thesis/plot-heterotypic-proliferation-regimes-vertex.py \
		segregated \
		data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv \
		$@

images/thesis/ch6-heterotypic-competition-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-heterotypic-competition-regimes.py
	python3 figures/thesis/plot-heterotypic-competition-regimes.py \
		I-II \
		$@

images/thesis/chapter-6/competition-regimes-transformed.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-heterotypic-competition-regimes.py
	python3 figures/thesis/plot-competition-regimes-transformed.py \
		$@

images/thesis/chapter-7/constant-emission.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-constant-emission.py
	python3 figures/thesis/plot-constant-emission.py \
		$@

# Paper image rules
images/paper/homotypic-proliferation-regimes-diagram.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-proliferation-regimes-diagram.py
	python3 figures/paper/plot-homotypic-proliferation-regimes-diagram.py \
		$@

images/paper/homotypic-proliferation-regimes-estimated.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-proliferation-regimes-estimated.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
	data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv
	python3 figures/paper/plot-homotypic-proliferation-regimes-estimated.py \
		data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
		data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv \
		$@

images/paper/heterotypic-proliferation-regimes-diagram.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-proliferation-regimes-diagram.py
	python3 figures/paper/plot-heterotypic-proliferation-regimes-diagram.py \
		$@

images/paper/heterotypic-proliferation-regimes-estimated-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-proliferation-regimes-estimated-well-mixed.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-asymptotic-I-data-02-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-II-data-02-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-III-data-02-08-2021.csv
	python3 figures/paper/plot-heterotypic-proliferation-regimes-estimated-well-mixed.py \
		data/heterotypic-model-2-exponential-asymptotic-I-data-02-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-II-data-02-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-III-data-02-08-2021.csv \
		$@

images/paper/heterotypic-proliferation-regimes-estimated-vertex-random.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv
	python3 figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
		random \
		data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv \
		$@

images/paper/heterotypic-proliferation-regimes-estimated-vertex-segregated.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
	data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv
	python3 figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
		segregated \
		data/heterotypic-model-2-exponential-asymptotic-I-tests-data-06-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-II-tests-data-09-08-2021.csv \
		data/heterotypic-model-2-exponential-asymptotic-III-tests-data-09-08-2021.csv \
		$@

images/paper/classification-of-competitive-interactions.png: \
	figures/paper/formatting.py \
	figures/paper/plot-classification-of-competitive-interactions.py
	python3 figures/paper/plot-classification-of-competitive-interactions.py \
		$@

images/paper/competition-regimes.png: \
	figures/paper/formatting.py \
	figures/paper/competition_regimes.py \
	figures/paper/plot-competition-regimes.py
	python3 figures/paper/plot-competition-regimes.py \
		I-II \
		$@

images/paper/competition-regimes-transformed.png: \
	figures/paper/formatting.py \
	figures/paper/competition_regimes.py \
	figures/paper/plot-competition-regimes-transformed.py
	python3 figures/paper/plot-competition-regimes-transformed.py \
		$@

images/paper/heterotypic-survival-difference-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-survival-difference-well-mixed.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-data-27-07-2021.csv
	python3 figures/paper/plot-heterotypic-survival-difference-well-mixed.py \
		data/heterotypic-model-2-exponential-data-27-07-2021.csv \
		$@

images/paper/heterotypic-survival-difference-vertex.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-survival-difference-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv
	python3 figures/paper/plot-heterotypic-survival-difference-vertex.py  \
		data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
		$@

images/paper/homotypic-survival-difference-%-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-difference-well-mixed.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-data-27-07-2021.csv \
	data/heterotypic-model-2-exponential-control-data-04-10-2021.csv
	python3 figures/paper/plot-homotypic-survival-difference-well-mixed.py \
		data/heterotypic-model-2-exponential-data-27-07-2021.csv \
		data/heterotypic-model-2-exponential-control-data-04-10-2021.csv \
		$* \
		$@

images/paper/homotypic-survival-difference-%-vertex.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-difference-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
	data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv
	python3 figures/paper/plot-homotypic-survival-difference-vertex.py \
		data/heterotypic-model-2-exponential-tests-data-04-08-2021.csv \
		data/heterotypic-model-2-exponential-control-tests-data-04-08-2021.csv \
		$* \
		$@

images/paper/homotypic-survival-frequency-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-frequency-well-mixed.py \
	data/homotypic-g2-death-signal-exponential-data-26-08-2021.csv
	python3 figures/paper/plot-homotypic-survival-frequency-well-mixed.py \
		data/homotypic-g2-death-signal-exponential-data-26-08-2021.csv \
		n \
		$@

images/paper/homotypic-survival-frequency-vertex.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-frequency-vertex.py \
	data/g1-proportion-tests-data-07-10-2021.csv
	python3 figures/paper/plot-homotypic-survival-frequency-vertex.py \
		data/g1-proportion-tests-data-07-10-2021.csv \
		$@ y
