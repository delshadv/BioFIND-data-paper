
% STEP 2
%% Features Extraction

% A multi-site, multi-participant magnetoencephalography resting-state dataset
% to study dementia: The BioFIND dataset
%
% This script calculates connectivity measure, and relative power
% for reproducing BioFIND data paper results.
% It requires ROI data (see "preproc_beamform_ROI.m" first)
%
% Note: Please make sure you are using the latest version of OSL and SPM
% within your OSL directory

% Henson R.N 2020, Vaghari D 2020

%% Define Paths and variables

% Assumed you are currently in a directory which includes BioFIND data,
% OSL and code directories as described in readme.md

bwd = pwd;
addpath(fullfile(bwd,'code')); % Assuming you've already downloaded BioFind
% repository in "code" directory.

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
processed_pth= fullfile(bwd,'Processed');

% Define Confounds and Covariate matrix
participants = spm_load(fullfile('code','participants-imputed.tsv'));
group_num    = grp2idx(participants.group);
site_num     = grp2idx(participants.site);
sex_num      = grp2idx(participants.sex);
mri_num      = grp2idx(participants.sImaging);

A_for_Cov = [site_num sex_num zscore([participants.age...
    participants.Move1 participants.Move2...
    participants.Recording_time]) ones(size(participants.age))];

% Cross-Validation setting
kFolds = 10; % Number of folds for cross-validation (inner and outer loops)
Nrun   = 100; % Number of runs to repeat cross-validation

%% Relative power in source space  %%

features = [];
parfor sub=1:length(subs)
    
    alphaRp = []; betaRp = []; totalp = []; locAlpha = []; locLow = [];locHigh = [];
    
    infile = fullfile(processed_pth,subdir{sub},'beffdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',3); % Data must be in ROI montage
    
    % Remove bad badtrials
    p = []; fP=[]; g=[];
    g = D(:,:,:);
    g(:,:,D.badtrials)=[];
    
    % Estimation of power spectral density 
    for e = 1:size(g,3)
        %[p(e,:,:),fP] = pwelch(g(:,:,e)',500,250,1000,D.fsample);
        [p(e,:,:),fP] = periodogram(g(:,:,e)',[],1000,D.fsample);
    end
    
    p = squeeze(mean(p,1)); p=p';
    
    locAlpha = find(fP>=8 & fP<=12);
    locBeta  = find(fP>12 & fP<=30);
    locLow   = find(fP==0.5);
    locHigh  = find(fP==48);
    
    totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp = sum((p(:,locAlpha)./totalp),2);
    betaRp  = sum((p(:,locBeta)./totalp),2);
    
    tmp = [alphaRp betaRp];
    features(sub,:) = tmp(:);
    tmp = [];
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
A = A(mri_num==1,:); % Remove subjects without T1 MRI
rng(1) % For reproducibility
%accuracy = repeated_CV(A,CVratio,kFolds,Nrun);
accuracy = repeated_CV_matlab(A,kFolds,Nrun);

mean(accuracy)
std(accuracy)

%% Relative power in sensor space %%

features = [];
parfor sub=1:length(subs)
    
    alphaRp = []; betaRp = []; totalp = []; locAlpha = []; locLow = [];locHigh = [];
    
    infile = fullfile(processed_pth,subdir{sub},'effdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',0); % Data must be in sensor montage
    
    % Remove bad badtrials
    p = []; fP=[]; g=[];
    chans = D.indchantype('MEGANY','GOOD'); % Retrieve all MEG channels
    g = D(chans,:,:);
    g(:,:,D.badtrials)=[];
    
    % Estimation of power spectral density
    for e = 1:size(g,3)
        %[p(e,:,:),fP] = pwelch(g(:,:,e)',500,250,1000,D.fsample);
        [p(e,:,:),fP] = periodogram(g(:,:,e)',[],1000,D.fsample);
    end
    
    p = squeeze(mean(p,1)); p=p';
    
    locAlpha = find(fP>=8 & fP<=12);
    locBeta  = find(fP>12 & fP<=30);
    locLow   = find(fP==0.5);
    locHigh  = find(fP==48);
    
    totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp = sum((p(:,locAlpha)./totalp),2);
    betaRp  = sum((p(:,locBeta)./totalp),2);
    
    tmp = [alphaRp betaRp];
    features(sub,:) = tmp(:);
    tmp = [];
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num]; % For all subjects
rng(1) % For reproducibility
accuracy = repeated_CV_matlab(A,kFolds,Nrun);

mean(accuracy)
std(accuracy)

A = A(mri_num==1,:); % Remove subjects without T1 MRI
rng(1) % For reproducibility
accuracy = repeated_CV_matlab(A,kFolds,Nrun);

mean(accuracy)
std(accuracy)

%% Connectivity measure: Amplitude Envelope Correlation %%

features  = [];
freqbands = {[8 13],[13 30]};

for ii = 1:length(freqbands)
    f=[];
    parfor sub = 1:length(subs)
        
        g = []; y = []; y0 = []; y1 = []; cm = [];
        infile = fullfile(processed_pth,subdir{sub},'beffdspmeeg');
        D = spm_eeg_load(infile);
        D = D.montage('switch',3); % Data must be in ROI montage
        
        % Remove bad badtrials
        g = D(:,:,:);
        g(:,:,D.badtrials)=[];
        
        % Remove source leakage using sysmetric Orthogonalisation
        y = reshape(g,D.nchannels,D.nsamples*size(g,3));  % generalise to Nsource,Ntime*Ntrials(good ones)
        y0 = ROInets.remove_source_leakage(y,'symmetric');
        
        % Filter to desired freq band
        y1 = ft_preproc_bandpassfilter(y0, D.fsample, freqbands{ii}, 4, 'but');
        y = reshape(y1,D.nchannels,D.nsamples,size(g,3));
        
        % Extract envelops
        Hen_lc_sep = [];
        for t=1:size(y,3)
            Hen_lc_sep(:,:,t) = hilbenv(squeeze(y(:,:,t)),1:D.nsamples,1,1);
        end
        
        Hen_lc_sep = reshape(Hen_lc_sep,D.nchannels,D.nsamples*size(g,3));
        
        % Calculate correlation
        cm = corr(resample(Hen_lc_sep',1,D.fsample));
        f(sub,:) = (cm(find(triu(cm,1))))';
        
    end
    features = [features f];
end
% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
A = A(mri_num==1,:); % Remove subjects without T1 MRI
rng(1) % For reproducibility
%accuracy = repeated_CV(A,CVratio,kFolds,Nrun);
accuracy = repeated_CV_matlab(A,kFolds,Nrun);

mean(accuracy)
std(accuracy)
