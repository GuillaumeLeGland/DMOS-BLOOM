# DMOS_BLOOM

Code and data for the Dynamical Model of Organic Sulfur production and consumption during phytoplankton blooms (DMOS-BLOOM)

The version DMOS_BLOOM_v1, currently the only available version, was used to write the article submitted to Limnology and Oceanography in early November 2022
by Guillaume Le Gland, Marta Masdeu-Navarro, Martí Galí, Sergio M. Vallina, Matti Gralka, Flora Vincent, Otto Cordero, Assaf Vardi and Rafel Simó.

DMOS-BLOOM is designed to run on both open-source GNU-Octave (version 4.4.1 and more recent) and MATLAB (version R2010b and more recent).
The execution has been tested on Windows with a 2.5 GHz Intel i5-3210M processor, on Linux Ubuntu with a 2.4 GHz Intel Xeon E5645 processor, and on Linux Debian with a 2.6 GHz Intel Xeon E5-2640 processor.

The main module is DMOS_BLOOM.m. The model options can be modified in DMOS_BLOOM_keys.m and DMOS_BLOOM_parameters.m.
The data necessary to run DMOS_BLOOM, mainly observations of DMSP and DMS concentrations and microbial abundances, are located in the DATA folder.

Please contact Guillaume Le Gland (legland@icm.csic.es) if you cannot download, run or understand the model.
