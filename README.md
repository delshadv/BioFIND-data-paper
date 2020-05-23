# A multi-site, multi-participant magnetoencephalography resting-state dataset to study dementia: The BioFIND dataset

Early detection of dementia, including Alzheimer’s Disease (AD), is vital to reduce the burden of disease and for developing effective treatments. Neuroimaging can detect early brain changes, such as hippocampal atrophy in prodromal AD presenting with Mild Cognitive Impairment (MCI). However, selecting and validating the most informative imaging features by machine-learning requires many participants (cases). While large publically-available datasets of people with dementia or prodromal disease exist for Magnetic Resonance Imaging (MRI), comparable datasets are missing for Magnetoencephalography (MEG). MEG offers advantages in its millisecond resolution, revealing physiological changes in brain oscillations or connectivity, before structural changes are evident with MRI.
In BioFIND dataset, we introduced a MEG dataset with 324 individuals, half patients with MCI and half age- and sex-matched healthy controls. Their brain activity was recorded while resting with eyes closed, using a 306-channel MEG scanner at one of two sites (Madrid or Cambridge), enabling tests of generalization across sites. A T1-weighted MRI is provided to assist source localisation. The MEG and MRI data are formatted according to international BIDS standards, and freely available from DPUK. 

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

The code directory content is: 
preproc_beamform_ROI.m – includes all steps from reading raw FIFF file to ROI extraction.
feature_extraction_test.m –  includes all steps from different features extraction to call an SVM classifier.
participants-imputed.tsv - is participants.tsv file which missing values of these covariates were imputed using the mean.

Work in progress 

rik.henson@mrc-cbu.cam.ac.uk
