# A multi-site, multi-participant magnetoencephalography resting-state dataset to study dementia: The BioFIND dataset

In BioFIND dataset, we introduced a MEG dataset with 324 individuals, half patients with MCI and half age- and sex-matched healthy controls. Their brain activity was recorded while resting with eyes closed, using a 306-channel MEG scanner at one of two sites (Madrid or Cambridge), enabling tests of generalization across sites. A T1-weighted MRI is provided to assist source localisation. The MEG and MRI data are formatted according to international BIDS standards. For more details please see (). Data are freely available from DPUK in ().


This repository contains all MATLAB scripts (and functions) including pre-processing, co-registration, source-localisation, feature extraction and machine learning steps to reproduce results of BioFIND data paper ().

Before downloading the repository:
To produce results succesfully, it needs that you download the updated OSL from (https://ohba-analysis.github.io/osl-docs/) which consequently includes SPM (https://www.fil.ion.ucl.ac.uk/spm/). Depending how you download OSL, you need to rename the OSL'root directory (after unziping/uncompressing) to osl and put it in your bwd directory. This should be a folder called osl with the following contents:

osl/

├── GLEAN/

├── HMM-MAR/

├── layouts/

├── MEG-ROI-nets/

├── ohba-external/

├── osl-core/

├── parcellations/

├── spm12/

├── std_masks/

└── version.txt

Then create a directory called "code" within bwd. So we will have a directory structure like this:

bwd/

├── osl/

├── code/


and download all scripts and function of this repository to the code directory.

The code directory content contains: 
preproc_beamform_ROI.m – includes all steps from reading raw FIFF file to ROI extraction.
feature_extraction_test.m –  includes all steps from different features extraction to call an SVM classifier.
participants-imputed.tsv - is participants.tsv file which missing values of covariates were imputed using the mean.

Work in progress 

rik.henson@mrc-cbu.cam.ac.uk
