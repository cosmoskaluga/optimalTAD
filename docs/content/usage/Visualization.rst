Visualization mode
==================

Hi-C data with the obtained optimal TAD set can be visualized using the function below:

.. code:: bash

	python3 -m optimalTAD visualize [-h] [--samplename SAMPLENAME] 
						[--region REGION] [--resolution RESOLUTION] 
						[--chipseq CHIPSEQ] [--log2_chip] 
						[--zscore_chip] [--rnaseq RNASEQ]


"""""""""""""""""""""""""""""""
``-h, --help``
"""""""""""""""""""""""""""""""
Help message

"""""""""""""""""""""""""""""""
``--samplename``
"""""""""""""""""""""""""""""""
Samplename of Hi-C data (for example, LacZ_1)

"""""""""""""""""""""""""""""""
``--region``
"""""""""""""""""""""""""""""""
Genome region to plot (for example, chr2L:1,000,000-5,000,000)

"""""""""""""""""""""""""""""""
``--resolution``
"""""""""""""""""""""""""""""""
Resolution of Hi-C matrix (=1)

"""""""""""""""""""""""""""""""
``--chipseq``
"""""""""""""""""""""""""""""""
Path to epigenetic data

"""""""""""""""""""""""""""""""
``--log2_chip``
"""""""""""""""""""""""""""""""
log2 transformation of epigenetic data

"""""""""""""""""""""""""""""""
``--zscore_chip``
"""""""""""""""""""""""""""""""
Z-score transformation of epigenetic data

"""""""""""""""""""""""""""""""
``--rnaseq``
"""""""""""""""""""""""""""""""
Add additional track to the plot



