# A multi-site, multi-participant magnetoencephalography resting-state dataset to study dementia: The BioFIND dataset

In BioFIND dataset, we introduced a MEG dataset with 324 individuals, half patients with MCI and half age- and sex-matched healthy controls. Their brain activity was recorded while resting with eyes closed, using a 306-channel MEG scanner at one of two sites (Madrid or Cambridge), enabling tests of generalization across sites. A T1-weighted MRI is provided to assist source localisation. The MEG and MRI data are formatted according to international BIDS standards. For more details please see the data paper in (). Data are freely available from DPUK in ().

This repository contains all MATLAB scripts (and functions) including pre-processing, co-registration, source-localisation, feature extraction and machine learning steps to reproduce results of data paper ().

Before downloading the repository:
To re-produce results succesfully, it is required that you first download and set up OSL toolbox from (https://ohba-analysis.github.io/osl-docs/) which consequently includes SPM toolbox (https://www.fil.ion.ucl.ac.uk/spm/). Depending on how you download OSL, you need to rename the OSL's root directory (after uncompression) to "osl" and put it in your current directory. This should be a folder called osl with the following contents:

osl/
.
.
.

├── MEG-ROI-nets/

├── osl-core/

├── parcellations/

├── spm12/
.
.
.

Then, create a directory called "code" at same level of osl means within your current directory and download all of this repository in code folder. So we will have a directory structure like this:

bwd/

├── osl/

├── code/


Repository content : 
preproc_beamform_ROI.m – includes all steps from reading raw FIFF file to ROI extraction.
feature_extraction_test.m –  includes the steps from different features extraction and call an SVM classifier and get nested cross-validation results (table 3, data paper).
participants-imputed.tsv - is MCIControls BIDS participants.tsv file in which missing values of covariates were imputed using the mean.

Work in progress 

rik.henson@mrc-cbu.cam.ac.uk
