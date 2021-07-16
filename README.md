# A multi-site, multi-participant magnetoencephalography resting-state dataset to study dementia: The BioFIND dataset
https://doi.org/10.48532/007000

In BioFIND dataset, we introduce a MEG dataset with 324 individuals, approximately half of whom are patients with MCI and the rest healthy controls. Brain activity was recorded while resting with eyes closed, using a 306-channel MEG scanner at one of two sites (Madrid or Cambridge), enabling tests of generalization across sites. A T1-weighted MRI is provided to assist source localisation. The MEG and MRI data are formatted according to the international BIDS standard. For more details please see the data paper in (PREPRINT LINK AVILABLE here https://www.medrxiv.org/content/10.1101/2021.05.19.21257330v1). 

Data are available to analyse on the DPUK platform: https://portal.dementiasplatform.uk/apply. They are stored in XNAT and BIDS format. The script "Xnat2Bids" can be run to convert from XNAT to BIDS format described in the data paper even though the BIDS format is already availlable.

This repository contains all MATLAB scripts (and functions) including pre-processing, co-registration, source-localisation, feature extraction and machine learning steps to reproduce results of data paper (https://www.medrxiv.org/content/10.1101/2021.05.19.21257330v1).

The BioFIND repository content: 

1. preproc_beamform_ROI.m – includes all steps from reading raw FIFF file to ROI extraction (Step 1)

2 .feature_extraction_test.m – includes all steps from feature extraction to calling an SVM classifier for cross-validated results (table 3, data paper). (Step 2)

Note: It is necessary to have MATLAB's "Digital Signal Processing" and "Statistics and Machine Learning" Toolboxes installed to run feature_extraction_test.m successfully which all are provided by DPUK.

3. participants-imputed.tsv - is MCIControl's participants.tsv file in which missing values of covariates were imputed using the mean (see data paper for more info about participants.tsv)

4. repeated_CV.m - a function for training an SVM model and produce cross-validation accuracy for "Nrun" permutations. (For more info about arguments see the function's description)
5. There is also instructions.md showing how to use data through Linux server on DPUK.

If you want to analyze data through DPUK you do not need to read below! All codes and data are available in https://portal.dpuk.ukserp.ac.uk/. Please first apply for data using https://portal.dementiasplatform.uk/apply and get necessary credentials.

Before downloading the repository:
It is required that you first download and set up OSL toolbox from (https://ohba-analysis.github.io/osl-docs/) which includes the SPM toolbox (https://www.fil.ion.ucl.ac.uk/spm/). Depending on how you download OSL, you need to rename the OSL's root directory to "osl" and put it in your current directory. This should be a folder called osl with at least the following contents:

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

Then, create a directory called "code" within your current directory (at same level of osl) and download all of this repository within code folde, so it looks like this:

bwd/

├── osl/

├── code/
.
.
.

Delshad Vaghari <dlshadv.ee@gmail.com>
Rik Henson <rik.henson@mrc-cbu.cam.ac.uk> 
