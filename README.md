# A multi-site, multi-participant magnetoencephalography resting-state dataset to study dementia: The BioFIND dataset

Early detection of dementia, including Alzheimer’s Disease (AD), is vital to reduce the burden of disease and for developing effective treatments. Neuroimaging can detect early brain changes, such as hippocampal atrophy in prodromal AD presenting with Mild Cognitive Impairment (MCI). However, selecting and validating the most informative imaging features by machine-learning requires many participants (cases). While large publically-available datasets of people with dementia or prodromal disease exist for Magnetic Resonance Imaging (MRI), comparable datasets are missing for Magnetoencephalography (MEG). MEG offers advantages in its millisecond resolution, revealing physiological changes in brain oscillations or connectivity, before structural changes are evident with MRI. We introduce a MEG dataset with 324 individuals, half patients with MCI and half age- and sex-matched healthy controls. Their brain activity was recorded while resting with eyes closed, using a 306-channel MEG scanner at one of two sites (Madrid or Cambridge), enabling tests of generalization across sites. A T1-weighted MRI is provided to assist source localisation. The MEG and MRI data are formatted according to international BIDS standards, and freely available from DPUK. 

This repository contains all MATLAB codes and scripts including preprocessing, source-localisation, feature-extraction and machine learning codes to reproduce results and analysis of BioFIND data paper ()

To run the scripts succesfully, it needs to have an updated OSL (https://ohba-analysis.github.io/osl-docs/) and
consequently SPM (https://www.fil.ion.ucl.ac.uk/spm/) be installed on your PC.

The main scripts are: 
feature_extraction_test.m –  includes all steps from feature extraction to call an SVM classifier.
preproc_beamform_ROI.m – includes all steps from reading raw FIFF file to ROI extractions.

Work in progress 

rik.henson@mrc-cbu.cam.ac.uk
