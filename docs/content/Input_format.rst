Input Format
=============

Both Hi-C and ChIP-Seq profiles are required for optimalTAD running. 

Hi-C data
----------
Hi-C contact matrices should be iteratively corrected and stored in the following formats:

- `.cool (.mcool)`. Cool is a widely usable binary Hi-C storage format introduced by Mirny lab.

.. code-block:: python

	import cooler
	c = cooler.Cooler('Filename.cool')
	c.info

.. code-block:: INI

	{'bin-size': 20000,
	 'bin-type': 'fixed',
	 'creation-date': '2020-06-02T14:03:32.896753',
	 'format': 'HDF5::Cooler',
	 'format-url': 'https://github.com/mirnylab/cooler',
	 'format-version': 3,
	 'generated-by': 'cooler-0.8.7',
	 'genome-assembly': 'unknown',
	 'metadata': {},
	 'nbins': 6024,
	 'nchroms': 7,
	 'nnz': 1903526,
	 'storage-mode': 'symmetric-upper',
	 'sum': 16003605}

- `.hdf5`. The algoritm also supports .hdf5 matrices, however the specific structure of these files is required. Individual chromosome Hi-C data must be stored in keys of the same names ('2L', '2R', etc). Chromosome names must be listed in the `chromosomeLabels` key and matrix resolution value must be indicated in the `resolution` key. Here is an example:

.. code-block:: python

	import h5py
	f = h5py.File("Filename.hdf5", 'r')
	f.keys()


.. code-block:: INI

	<KeysViewHDF5 ['chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chromosomeIndex', 'chromosomeLabels', 'chromosomeStarts', 'genome', 'positionIndex', 'resolution']>



ChIP-seq data
-------------

- `.bedgraph`. optimalTAD supports a classical bedgraph format consisting of 4 columns: `chromName`, `chromStart`, `chromEnd`, `dataValue`. Columns must be separated by a space (' '). 

.. code-block:: python

	import pandas as pd
	data = pd.read_csv("Filename.bedgraph", sep = ' ', header = None, names=['Chr', 'Start', 'End', 'Score'])
	data.head()

.. code-block:: INI
	
	Chr	Start	End	Score
	0	chr2L	0	20000	1.904566
	1	chr2L	20000	40000	2.963382
	2	chr2L	40000	60000	2.944759
	3	chr2L	60000	80000	4.394352
	4	chr2L	80000	100000	3.742936	


- `.bw` (`.BigWig` is also accepted) 


