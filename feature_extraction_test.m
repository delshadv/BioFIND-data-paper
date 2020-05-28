
% STEP 2
%% Features Extraction

% A multi-site, multi-participant magnetoencephalography resting-state dataset
% to study dementia: The BioFIND dataset
%
% This script calculates connectivity measure, and relative power
% for reproducing BioFIND data paper results.
% It requires ROI data (see "preproc_beamform_ROI.m" fisr)
%
% Note: Please make sure you are using the latest version of OSL and SPM
% within your OSL directory
% Note: Please copy all files from BioFIND repository (codes directory)

% Henson R.N 2020, Vaghari D 2020

%% Define Paths and variables

% Assuming you are currently in a directory which includes BioFIND data,
% OSL and code directories as described in readme.md

bwd = pwd;
wd  = fullfile(bwd,'code'); % Assuming you've already downloaded BioFind
% repository in the "code" directory.
cd (wd)

% Setup OSL
addpath(fullfile(bwd,'osl','osl-core'))
osl_startup
osl_check_installation

% BIDS and Processed directories
bidspth = fullfile(bwd,'BioFIND','MCIControls'); %BIDS Path
BIDS   = spm_BIDS(bidspth); % (If it does not work with OSL's SPM, so copy last version of spm_BIDS)
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

% Define Confounds and Covariate matrix
participants = spm_load(fullfile(wd,'participants-imputed.tsv'));
group_num    = grp2idx(participants.group);
site_num     = grp2idx(participants.site);
sex_num      = grp2idx(participants.sex);
Age          = participants.age;
mri_num      = grp2idx(participants.sImaging);
rec_time     = participants.Recording_time;
Move1        = participants.Move1;
Move2        = participants.Move2;
A_for_Cov = [site_num sex_num zscore([Age Move1 Move2 rec_time]) ones(size(Age))];

%% Cross-Validation setting

kFolds = 10; % Number of folds for nested cross-validation (inner and outer)

%% Relative power in source space %%

features = [];
parfor sub=1:length(subs)
    
    alphaRp = []; totalp = []; locAlpha = []; locLow = [];locHigh = [];
    
    infile = fullfile(processed_pth,subdir{sub},'beffdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',3); % Data must be in ROI montage
    
    % Estimation of power spectral density
    p = []; fP=[]; g=[]; y=[];
    g=D(:,:,:);
    g(:,:,D.badtrials)=[];
    
    for e = 1:size(g,3)
        %[p(e,:,:),fP] = pwelch(g(:,:,e)',500,450,1000,D.fsample);
        [p(e,:,:),fP] = periodogram(g(:,:,e)',[],1000,D.fsample);
    end
    
    p = squeeze(mean(p,1)); p=p';
    
    locAlpha = find(fP>=8 & fP<=12);
    locLow   = find(fP==0.5);
    locHigh  = find(fP==48);
    
    totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp = sum((p(:,locAlpha)./totalp),2);
    
    
    features(sub,:) = alphaRp;
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
A = A(mri_num==1,:); % Remove subjects without T1 MRI
rng(1) % For reproducibility
accuracy = nestedCVSP(A,kFolds)

%% Connectivity measure: Amplitude Envelope Correlation %%

features = [];
parfor sub = 1:length(subs)
    
    infile = fullfile(processed_pth,subdir{sub},'beffdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',3); % Data must be in ROI montage
    
    % Remove source leakage using sysmetric Orthogonalisation
    y = D(:,:,:);
    y = reshape(y,D.nchannels,D.nsamples*D.ntrials); % generalise to Nsource,Ntime*Ntrials
    
    y0 = ft_preproc_bandpassfilter(y, D.fsample, [6 14], 4, 'but', 'twopass', 'no');
    y0 = ROInets.remove_source_leakage(y0,'symmetric');

    y = reshape(y0,D.nchannels,D.nsamples,D.ntrials);
    
    Hen_lc_sep = [];
    for t=1:size(y,3)
        Hen_lc_sep(:,:,t) = hilbenv(squeeze(y(:,:,t)),1:D.nsamples,1,1);
    end
    
    Hen_lc_sep = reshape(Hen_lc_sep,D.nchannels,D.nsamples*D.ntrials);
    
    cm = corr(resample(Hen_lc_sep',1,D.fsample)); 
    features(sub,:) = (cm(find(triu(cm,1))))';
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
A = A(mri_num==1,:); % Remove subjects without T1 MRI
rng(1) % For reproducibility
accuracy = nestedCVCon(A,kFolds)

%% Relative power in sensor space %%

features = [];
parfor sub=1:length(subs)
    
    alphaRp = []; totalp = []; locAlpha = []; locLow = [];locHigh = [];
    
    infile = fullfile(processed_pth,subdir{sub},'effdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',0); % Data must be in sensor montage
    
    % Estimation of power spectral density
    p = []; fP=[]; g=[]; y=[];
    chans = D.indchantype('MEGANY','GOOD'); % Retrieve all MEG channels
    g = D(chans,:,:);
    g(:,:,D.badtrials)=[];
    
    
    for e = 1:size(g,3)
        %[p(e,:,:),fP] = pwelch(g(:,:,e)',500,450,1000,D.fsample);
        [p(e,:,:),fP] = periodogram(g(:,:,e)',[],1000,D.fsample);
    end
    
    p = squeeze(mean(p,1)); p=p';
    
    locAlpha = find(fP>=8 & fP<=12);
    locLow   = find(fP==0.5);
    locHigh  = find(fP==48);
    
    totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp = sum((p(:,locAlpha)./totalp),2);
    
    
    features(sub,:) = alphaRp;
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
rng(1) % to reproduce paper results
