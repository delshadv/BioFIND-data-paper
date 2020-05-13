
%% Features Extraction
% A multi-site, multi-participant magnetoencephalography resting-state dataset
% to study dementia: The BioFIND dataset

% This script calculates connectivity measure, and relative power for reproducing BioFIND data
% paper results. It requires ROI data (see "preproc_beamform_ROI.m")

% Note: Please make sure you are using the latest version of OSL and SPM
% within your OSL directory

% Henson R.N 2020, Vaghari D 2020

%% Setup OSL
addpath('/imaging/dv01/Myscripts/BioFIND_datapaper')% Path to codes directory
addpath('/imaging/dv01/Toolboxes/osl/osl-core') % Path to OSL directory
osl_startup
osl_check_installation

%% Adding BIDS Paths, defining BIDS variables

% Please do all analysis in a sperate directory from your BIDS data. Here, we
% call it processed_pth
processed_pth= '/imaging/dv01/Processed_sss_new';
bidspth = '/imaging/dv01/BioFIND/MCIControls'; % BIDS Path
BIDS   = spm_BIDS(bidspth); % (If it does not work with OSL's SPM, so copied last version of spm_BIDS)
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);

% Define Confound, Covariate matrix

participants = spm_load('/imaging/dv01/Myscripts/participants-imputed.tsv');
group_num    = grp2idx(participants.group)-1;
site_num     = grp2idx(participants.site)-1;
sex_num      = grp2idx(participants.sex)-1;
Age          = participants.age;
mri_num      = grp2idx(participants.sImaging);
rec_time     = participants.Recording_time;
Edu          = participants.Edu_years;
Move1        = participants.Move1;
Move2        = participants.Move2;

A_for_Cov = [site_num sex_num zscore([Age Move1 Move2 rec_time]) ones(size(Age))];

%% Connectivity measure: Amplitude Envelope Correlation %%

features = [];
parfor sub = 1:length(subs)
    
    infile = fullfile(processed_pth,subdir{sub},'beffdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',3); % Data must be in ROI montage
    
    % Remove source leakage using sysmetric Orthogonalisation
    y = D(:,:,:);
    y = reshape(y,D.nchannels,D.nsamples*D.ntrials); % generalise to Nsource,Ntime*Ntrials
    y0 = ROInets.remove_source_leakage(y,'symmetric');
    
    % Filter to desired frequency band
    y0 = ft_preproc_bandpassfilter(y0, D.fsample, [6 14], 4, 'but', 'twopass', 'no');
    y = reshape(y0,D.nchannels,D.nsamples,D.ntrials);
    
    Hen_lc_sep = [];
    for t=1:size(y,3)
        Hen_lc_sep(:,:,t) = hilbenv(squeeze(y(:,:,t)),1:D.nsamples,1,1);
    end
    Hen_lc_sep = reshape(Hen_lc_sep,D.nchannels,D.nsamples*D.ntrials);
    
    cm = corr(resample(Hen_lc_sep',1,D.fsample)); %+diag(nan(38,1));
    %         ca = [min(cm(:)) max(cm(:))];
    %         figure;
    %         imagesc(cm), caxis(ca); colorbar;
    features(sub,:) = (cm(find(triu(cm,1))))';
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M;
A = [aM group_num];
A = A(mri_num==1,:);
% to reproduce paper results, created by Classifier Learner app
[~, validationAccuracy] = trainClassifierC(A)

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
        [p(e,:,:),fP]=pwelch(g(:,:,e)',500,400,1000,D.fsample);
    end
    p = squeeze(mean(p,1)); p=p';
    
    locAlpha = find(fP>=8 & fP<=12);
    locLow   = find(fP==0.5);
    locHigh  = find(fP==150);
    
    totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp = sum((p(:,locAlpha)./totalp),2);
    
    features(sub,:) = alphaRp;
    
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M;
A = [aM group_num];
A = A(mri_num==1,:);

% to reproduce paper results, created by Classifier Learner app

[~, validationAccuracy] = trainClassifierPsource(A)



%% Relative power in sensor space %%

features = [];

parfor sub=1:length(subs)
    
    infile = fullfile(processed_pth,subdir{sub},'beffdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',0); % Data must be in sensor montage
    
    % Estimation of power spectral density
    p = []; fP=[]; g=[]; y=[];
    chans = D.indchantype('MEGANY','GOOD'); % Retrieve all MEG channels
    g = D(chans,:,:);
    g(:,:,D.badtrials)=[];
    
    
    for e = 1:size(g,3)
        [p(e,:,:),fP]=pwelch(g(:,:,e)',500,400,1000,D.fsample);
    end
    p = squeeze(mean(p,1)); p=p';
    
    locAlpha = find(fP>=8 & fP<=12);
    locLow   = find(fP==0.5);
    locHigh  = find(fP==150);
    
    totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp = sum((p(:,locAlpha)./totalp),2);
    
    features(sub,:) = alphaRp;
end

% Machine Learning
M   = features;
aM = M - A_for_Cov*pinv(A_for_Cov)*M;
A = [aM group_num];
% to reproduce paper results, created by Classifier Learner app
[~, validationAccuracy] = trainClassifierPsensor(A)
