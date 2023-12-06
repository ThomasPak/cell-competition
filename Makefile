.PHONY: all thesis_images paper_images clean

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
			 images/thesis/chapter-5/homotypic-proliferation-regimes.png \
			 images/thesis/chapter-5/homotypic-proliferation-regimes-estimated.png \
			 images/thesis/chapter-5/quadratic-bounds.png \
			 images/thesis/chapter-5/parameter-regions-uniform.png \
			 images/thesis/chapter-5/homotypic-survival-probability-uniform.png \
			 images/thesis/chapter-5/uniform-simulations.png \
			 images/thesis/chapter-5/homotypic-survival-frequency-well-mixed-uniform.png \
			 images/thesis/chapter-5/homotypic-survival-frequency-vertex-uniform.png \
			 images/thesis/chapter-5/homotypic-proliferation-regimes-uniform.png \
			 images/thesis/chapter-5/homotypic-proliferation-regimes-estimated-uniform.png \
			 images/thesis/chapter-6/heterotypic-survival-difference-well-mixed.png \
			 images/thesis/chapter-6/heterotypic-survival-difference-vertex.png \
			 images/thesis/chapter-6/homotypic-survival-difference-A-well-mixed.png \
			 images/thesis/chapter-6/homotypic-survival-difference-B-well-mixed.png \
			 images/thesis/chapter-6/homotypic-survival-difference-A-vertex.png \
			 images/thesis/chapter-6/homotypic-survival-difference-B-vertex.png \
			 images/thesis/chapter-6/classification-of-competitive-interactions.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes-well-mixed.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes-vertex-random.png \
			 images/thesis/chapter-6/heterotypic-proliferation-regimes-vertex-segregated.png \
			 images/thesis/chapter-6/competition-regimes.png \
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
		 images/paper/supplementary/heterotypic-survival-difference-well-mixed.png \
		 images/paper/supplementary/heterotypic-survival-difference-vertex.png \
		 images/paper/supplementary/homotypic-survival-difference-A-well-mixed.png \
		 images/paper/supplementary/homotypic-survival-difference-B-well-mixed.png \
		 images/paper/supplementary/homotypic-survival-difference-A-vertex.png \
		 images/paper/supplementary/homotypic-survival-difference-B-vertex.png \
		 images/paper/supplementary/homotypic-survival-frequency-well-mixed.png \
		 images/paper/supplementary/homotypic-survival-frequency-vertex.png \
		 images/paper/supplementary/homotypic-simulation-eta-0.2.png \
		 images/paper/supplementary/homotypic-simulation-eta-0.05.png \
		 images/paper/supplementary/heterotypic-simulation.png

all: thesis_images paper_images

clean:
	rm -rf images

thesis_images: $(THESISIMAGES)

paper_images: $(PAPERIMAGES)

# create directories
$(shell mkdir -p images/thesis/chapter-2)
$(shell mkdir -p images/thesis/chapter-4)
$(shell mkdir -p images/thesis/chapter-5)
$(shell mkdir -p images/thesis/chapter-6)
$(shell mkdir -p images/thesis/chapter-7)
$(shell mkdir -p images/paper/supplementary)

# Thesis image rules
images/thesis/chapter-2/survival-frequency-difference-histograms.png: \
	data/mechanical-grid-heterotypic-data.csv \
	data/mechanical-grid-homotypic-data.csv \
	figures/thesis/mechanical-cell-competition-analysis.py
	python3 figures/thesis/mechanical-cell-competition-analysis.py \
		data/mechanical-grid-heterotypic-data.csv  \
		data/mechanical-grid-homotypic-data.csv \
		diff-theta-histogram \
		$@

images/thesis/chapter-2/pearson-correlation.png: \
	data/mechanical-grid-heterotypic-data.csv \
	data/mechanical-grid-homotypic-data.csv \
	figures/thesis/mechanical-cell-competition-analysis.py
	python3 figures/thesis/mechanical-cell-competition-analysis.py \
		data/mechanical-grid-heterotypic-data.csv  \
		data/mechanical-grid-homotypic-data.csv \
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

images/thesis/chapter-5/homotypic-proliferation-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-proliferation-regimes.py
	python3 figures/thesis/plot-homotypic-proliferation-regimes.py \
		$@

images/thesis/chapter-5/homotypic-proliferation-regimes-estimated.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-proliferation-regimes-estimated.py \
	figures/thesis/plot_survival_difference.py \
	data/homotypic-proliferation-regimes-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-vertex-data.csv
	python3 figures/thesis/plot-homotypic-proliferation-regimes-estimated.py \
		data/homotypic-proliferation-regimes-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-vertex-data.csv \
		$@

images/thesis/chapter-5/quadratic-bounds.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-quadratic-bounds.py
	python3 figures/thesis/plot-quadratic-bounds.py \
		$@

images/thesis/chapter-5/parameter-regions-uniform.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-parameter-regions-uniform.py
	python3 figures/thesis/plot-parameter-regions-uniform.py \
		$@

images/thesis/chapter-5/homotypic-survival-probability-uniform.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-probability-uniform.py
	python3 figures/thesis/plot-homotypic-survival-probability-uniform.py \
		$@

images/thesis/chapter-5/uniform-simulations.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-uniform-simulations.py \
	data/D-simulation.csv \
	data/E-simulation.csv \
	data/E-low-eta-simulation.csv
	python3 figures/thesis/plot-uniform-simulations.py \
		data/D-simulation.csv \
		data/E-simulation.csv \
		data/E-low-eta-simulation.csv \
		$@

images/thesis/chapter-5/homotypic-survival-frequency-well-mixed-uniform.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-frequency-well-mixed-uniform.py \
	data/homotypic-survival-frequency-uniform-well-mixed-data.csv
	python3 figures/thesis/plot-homotypic-survival-frequency-well-mixed-uniform.py \
		data/homotypic-survival-frequency-uniform-well-mixed-data.csv \
		$@

images/thesis/chapter-5/homotypic-survival-frequency-vertex-uniform.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-frequency-vertex-uniform.py \
	data/homotypic-survival-frequency-vertex-data.csv
	python3 figures/thesis/plot-homotypic-survival-frequency-vertex-uniform.py \
		data/homotypic-survival-frequency-vertex-data.csv \
		$@

images/thesis/chapter-5/homotypic-proliferation-regimes-uniform.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-proliferation-regimes-uniform.py
	python3 figures/thesis/plot-homotypic-proliferation-regimes-uniform.py \
		$@

images/thesis/chapter-5/homotypic-proliferation-regimes-estimated-uniform.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-proliferation-regimes-estimated-uniform.py \
	figures/thesis/plot_survival_difference.py \
	data/homotypic-proliferation-regimes-uniform-cross-section-I-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-uniform-cross-section-II-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-uniform-cross-section-I-vertex-data.csv \
	data/homotypic-proliferation-regimes-uniform-cross-section-II-vertex-data.csv
	python3 figures/thesis/plot-homotypic-proliferation-regimes-estimated-uniform.py \
		data/homotypic-proliferation-regimes-uniform-cross-section-I-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-uniform-cross-section-II-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-uniform-cross-section-I-vertex-data.csv \
		data/homotypic-proliferation-regimes-uniform-cross-section-II-vertex-data.csv \
		$@

images/thesis/chapter-6/heterotypic-survival-difference-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-survival-difference-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-survival-difference-well-mixed-data.csv
	python3 figures/thesis/plot-heterotypic-survival-difference-well-mixed.py \
		data/heterotypic-survival-difference-well-mixed-data.csv \
		$@

images/thesis/chapter-6/heterotypic-survival-difference-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-survival-difference-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-survival-difference-vertex-data.csv
	python3 figures/thesis/plot-heterotypic-survival-difference-vertex.py  \
		data/heterotypic-survival-difference-vertex-data.csv \
		$@

images/thesis/chapter-6/homotypic-survival-difference-A-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-survival-difference-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-well-mixed-data.csv
	python3 figures/thesis/plot-homotypic-survival-difference-well-mixed.py \
		data/heterotypic-survival-difference-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-well-mixed-data.csv \
		A \
		$@

images/thesis/chapter-6/homotypic-survival-difference-B-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-survival-difference-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-well-mixed-data.csv
	python3 figures/thesis/plot-homotypic-survival-difference-well-mixed.py \
		data/heterotypic-survival-difference-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-well-mixed-data.csv \
		B \
		$@

images/thesis/chapter-6/homotypic-survival-difference-A-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-survival-difference-vertex-data.csv \
	data/homotypic-proliferation-regimes-vertex-data.csv
	python3 figures/thesis/plot-homotypic-survival-difference-vertex.py \
		data/heterotypic-survival-difference-vertex-data.csv \
		data/homotypic-proliferation-regimes-vertex-data.csv \
		A \
		$@

images/thesis/chapter-6/homotypic-survival-difference-B-vertex.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-homotypic-survival-difference-vertex.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-survival-difference-vertex-data.csv \
	data/homotypic-proliferation-regimes-vertex-data.csv
	python3 figures/thesis/plot-homotypic-survival-difference-vertex.py \
		data/heterotypic-survival-difference-vertex-data.csv \
		data/homotypic-proliferation-regimes-vertex-data.csv \
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

images/thesis/chapter-6/heterotypic-proliferation-regimes-well-mixed.png: \
	figures/thesis/formatting.py \
	figures/thesis/plot-heterotypic-proliferation-regimes-well-mixed.py \
	figures/thesis/plot_survival_difference.py \
	data/heterotypic-proliferation-regimes-cross-section-I-well-mixed-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-II-well-mixed-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-III-well-mixed-data.csv
	python3 figures/thesis/plot-heterotypic-proliferation-regimes-well-mixed.py \
		data/heterotypic-proliferation-regimes-cross-section-I-well-mixed-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-II-well-mixed-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-III-well-mixed-data.csv \
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

images/thesis/chapter-6/heterotypic-proliferation-regimes-vertex-segregated.png: \
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

images/thesis/chapter-6/competition-regimes.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-competition-regimes.py
	python3 figures/thesis/plot-competition-regimes.py \
		I-II \
		$@

images/thesis/chapter-6/competition-regimes-transformed.png: \
	figures/thesis/formatting.py \
	figures/thesis/competition_regimes.py \
	figures/thesis/plot-competition-regimes-transformed.py
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
	data/homotypic-proliferation-regimes-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-vertex-data.csv
	python3 figures/paper/plot-homotypic-proliferation-regimes-estimated.py \
		data/homotypic-proliferation-regimes-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-vertex-data.csv \
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
	data/heterotypic-proliferation-regimes-cross-section-I-well-mixed-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-II-well-mixed-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-III-well-mixed-data.csv
	python3 figures/paper/plot-heterotypic-proliferation-regimes-estimated-well-mixed.py \
		data/heterotypic-proliferation-regimes-cross-section-I-well-mixed-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-II-well-mixed-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-III-well-mixed-data.csv \
		$@

images/paper/heterotypic-proliferation-regimes-estimated-vertex-random.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv
	python3 figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
		random \
		data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv \
		$@

images/paper/heterotypic-proliferation-regimes-estimated-vertex-segregated.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
	data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv
	python3 figures/paper/plot-heterotypic-proliferation-regimes-estimated-vertex.py \
		segregated \
		data/heterotypic-proliferation-regimes-cross-section-I-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-II-vertex-data.csv \
		data/heterotypic-proliferation-regimes-cross-section-III-vertex-data.csv \
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

# Supplementary material of paper
images/paper/supplementary/heterotypic-survival-difference-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-survival-difference-well-mixed.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-survival-difference-well-mixed-data.csv
	python3 figures/paper/plot-heterotypic-survival-difference-well-mixed.py \
		data/heterotypic-survival-difference-well-mixed-data.csv \
		$@

images/paper/supplementary/heterotypic-survival-difference-vertex.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-survival-difference-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-survival-difference-vertex-data.csv
	python3 figures/paper/plot-heterotypic-survival-difference-vertex.py  \
		data/heterotypic-survival-difference-vertex-data.csv \
		$@

images/paper/supplementary/homotypic-survival-difference-%-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-difference-well-mixed.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-survival-difference-well-mixed-data.csv \
	data/homotypic-proliferation-regimes-well-mixed-data.csv
	python3 figures/paper/plot-homotypic-survival-difference-well-mixed.py \
		data/heterotypic-survival-difference-well-mixed-data.csv \
		data/homotypic-proliferation-regimes-well-mixed-data.csv \
		$* \
		$@

images/paper/supplementary/homotypic-survival-difference-%-vertex.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-difference-vertex.py \
	figures/paper/plot_survival_difference.py \
	data/heterotypic-survival-difference-vertex-data.csv \
	data/homotypic-proliferation-regimes-vertex-data.csv
	python3 figures/paper/plot-homotypic-survival-difference-vertex.py \
		data/heterotypic-survival-difference-vertex-data.csv \
		data/homotypic-proliferation-regimes-vertex-data.csv \
		$* \
		$@

images/paper/supplementary/homotypic-survival-frequency-well-mixed.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-frequency-well-mixed.py \
	data/homotypic-survival-frequency-well-mixed-data.csv
	python3 figures/paper/plot-homotypic-survival-frequency-well-mixed.py \
		data/homotypic-survival-frequency-well-mixed-data.csv \
		n \
		$@

images/paper/supplementary/homotypic-survival-frequency-vertex.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-survival-frequency-vertex.py \
	data/homotypic-survival-frequency-vertex-data.csv
	python3 figures/paper/plot-homotypic-survival-frequency-vertex.py \
		data/homotypic-survival-frequency-vertex-data.csv \
		$@ y

images/paper/supplementary/homotypic-simulation-eta-%.png: \
	figures/paper/formatting.py \
	figures/paper/plot-homotypic-simulation.py \
	data/eta-%-simulation.csv
	python3 figures/paper/plot-homotypic-simulation.py \
		data/eta-$*-simulation.csv \
		$@

images/paper/supplementary/heterotypic-simulation.png: \
	figures/paper/formatting.py \
	figures/paper/plot-heterotypic-simulation.py \
	data/heterotypic-simulation.csv
	python3 figures/paper/plot-heterotypic-simulation.py \
		data/heterotypic-simulation.csv \
		$@
