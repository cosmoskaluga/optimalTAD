Optimization mode
===============================

To launch the algorithm  type the following at the command line:

.. code:: bash

	python3 -m optimalTAD run [-h] [--hic HIC [HIC ...]] [--chipseq CHIPSEQ [CHIPSEQ ...]] 
			[--np NP] [--resolution RESOLUTION] [--stepsize STEPSIZE] 
			[--gamma_max GAMMA_MAX] [--hic_format HIC_FORMAT] 
			[--empty_row_imputation] [--truncation] [--log2_hic] 
			[--log2_chip] [--zscore_chip]


"""""""""""""""""""""""""""""""
``-h, --help``
"""""""""""""""""""""""""""""""
Help message

"""""""""""""""""""""""""""""""
``--hic``
"""""""""""""""""""""""""""""""
Iteratively corrected Hi-C matrices in .hdf5 or .cool format

"""""""""""""""""""""""""""""""
``--chipseq``
"""""""""""""""""""""""""""""""
Epigenetic data (ChIP-seq in .bedgraph or .bw format)

"""""""""""""""""""""""""""""""
``--np``
"""""""""""""""""""""""""""""""
Number of processors (=1)

"""""""""""""""""""""""""""""""
``--resolution``
"""""""""""""""""""""""""""""""
Resolution of Hi-C matrices (=1)

"""""""""""""""""""""""""""""""
``--stepsize``
"""""""""""""""""""""""""""""""
Step size to increment gamma parameter in Armatus (=0.05)

"""""""""""""""""""""""""""""""
``--gamma_max``
"""""""""""""""""""""""""""""""
Max value of the gamma parameter (=4)

"""""""""""""""""""""""""""""""
``--hic_format``
"""""""""""""""""""""""""""""""
Hi-C matrices input format for armatus (=txt.gz)

"""""""""""""""""""""""""""""""
``--empty_row_imputation``
"""""""""""""""""""""""""""""""
Empty line imputation (=False)

"""""""""""""""""""""""""""""""
``--truncation``
"""""""""""""""""""""""""""""""
Truncation of a Hi-C-matrix (=False)

"""""""""""""""""""""""""""""""
``--log2_hic``
"""""""""""""""""""""""""""""""
log2 transformation of Hi-C matrix (=False)

"""""""""""""""""""""""""""""""
``--log2_chip``
"""""""""""""""""""""""""""""""
log2 transformation of ChIP-seq values (=False)

"""""""""""""""""""""""""""""""
``--zscore_chip``
"""""""""""""""""""""""""""""""
Z-score transformation of ChIP-seq values(=False)

