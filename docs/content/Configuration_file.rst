Configuration file
==================

Below, we provide a template of the configutation file required for `optimalTAD`. While this template is needed for both 'optimization' and 'visualization' modes, `[visualization]` section is targeted for 'visualization' mode and doesn't affect TAD set optimization process. 

.. code-block:: INI

	# .ini basic configuration file

	[basic]
	mode = run

	[run]
	hic = ./testdata/LacZ1_2L.hdf5
	chipseq = ./testdata/LacZ1_2L.bedgraph  
	np = 3
	resolution = 20000
	stepsize = 0.05
	gamma_max = 4
	hic_format = txt.gz
	truncation = True
	log2_hic = True
	log2_chip = True
	zscore_chip = True
	empty_row_imputation = False

	[chromosomes]
	set_chromosomes = chr2L # should be separated by commas

	[hic]
	shrinkage_min = 0.5
	shrinkage_max = 1024
	balance = True

	[stair]
	index_min = -4
	index_max = 5
	acetyl_min = -3 
	acetyl_max = 5

	[output]
	path_to_amplitude_figure = output/figures/StairAmplitude
	path_to_amplitude_file = output/amplitudes.csv
	path_to_best_stair_data = output/stair.csv
	path_to_stair_data = output/all_stairs.csv
	path_to_best_stair_figure = output/figures/BestStairs
	figure_postfix = png
	figure_dpi = 300

	[visualization]
	samplename = LacZ_1
	region = chr2L:0-6,000,000
	path_to_chip = ./testdata/LacZ1_2L.bedgraph
	rnaseq = False
	fontsize = 12
	nticks = 5
	cmap = coolwarm
	hic_text = Hi-C
	chip_text = H3
	rnaseq_text = RNA-Seq
	vline_linewidth = 1
	vline_linestyle = dashed
	tad_linewidth = 1.2
	tad_linestyle = --
	filename = output/figures/TADandChIPseq.png
	dpi = 300











