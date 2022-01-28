
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

addpath(fullfile(bwd,'code'))

% Define Confounds and Covariate matrix
participants = spm_load(fullfile('code','participants-imputed.tsv'));
group_num    = grp2idx(participants.group);
site_num     = grp2idx(participants.site);
sex_num      = grp2idx(participants.sex);
mri_num      = grp2idx(participants.sImaging);

covars = {[],[site_num zscore([participants.Move1 participants.Move2...
    participants.Recording_time PreT number_bad_epochs]) ones(size(participants.age))],...
    [site_num sex_num zscore([participants.Move1 participants.Move2 participants.age ...
    participants.Recording_time PreT number_bad_epochs edu]) ones(size(participants.age))]};

% Cross-Validation setting
kFolds = 10; % Number of folds for cross-validation (inner and outer loops)
Nrun   = 100; % Number of runs to repeat cross-validation
    
%% Relative power in sensor space %%

features = [];
parfor sub=1:nsub
    
    
    infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'effdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',0); % Data must be in sensor montage
    
    % Remove bad badtrials
    p1 = []; p2 = []; fP = [];
    chans1 = D.indchantype('MEGMAG','GOOD'); % Retrieve all GRADs channels
    chans2 = D.indchantype('MEGPLANAR','GOOD'); % Retrieve all GRADs channels
    g1 = D(chans1,:,:);
    g1(:,:,D.badtrials)=[];
    g2 = D(chans2,:,:);
    g2(:,:,D.badtrials)=[];
    
    % Estimation of power spectral density
    for e = 1:size(g1,3)
        %[p(e,:,:),fP] = pwelch(g(:,:,e)',500,250,1000,D.fsample);
        [p1(e,:,:),fP] = periodogram(g1(:,:,e)',[],1000,D.fsample);
        [p2(e,:,:),~] = periodogram(g2(:,:,e)',[],1000,D.fsample);
    end
    
    p1 = squeeze(mean(p1,1)); p1 = p1';
    p2 = squeeze(mean(p2,1)); p2 = p2';

    p1 = p1./sum(p1,1);
    p2 = p2./sum(p2,1);
    
    locAlpha = find(fP>=8 & fP<=12);
    locBeta  = find(fP>12 & fP<=30);
    locLow   = find(fP==2);
    locHigh  = find(fP==48);
    
    %totalp  = sum(p(:,locLow:locHigh),2);
    alphaRp1 = sum((p1(:,locAlpha)),2);
    betaRp1  = sum((p1(:,locBeta)),2);
    
    alphaRp2 = sum((p2(:,locAlpha)),2);
    betaRp2  = sum((p2(:,locBeta)),2);
    
    tmp = [alphaRp1' alphaRp2' betaRp1' betaRp2']; 
    features(sub,:) = tmp(:);
    tmp = [];
    
    
end

M   = features;

for nn = 1 : length(covars)
    
    A = [];
    A_for_Cov = covars{nn};
    % Machine Learning
    
    if isempty(A_for_Cov)
        A = [M group_num];
        rng(1) % For reproducibility
        accuracy = repeated_CV_matlab(A,kFolds,Nrun);
        mean(accuracy)
        std(accuracy)
    else
        aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
        A = [aM group_num]; % For all subjects
        rng(1) % For reproducibility
        accuracy = repeated_CV_matlab(A,kFolds,Nrun);
        mean(accuracy)
        std(accuracy)
        
        if nn == 2
            A = A(mri_num==1,:); % Remove subjects without T1 MRI
            rng(1) % For reproducibility
            accuracy = repeated_CV_matlab(A,kFolds,Nrun);
            mean(accuracy)
            std(accuracy)
        end
        
        
    end
    
end


%% Relative power in source space  %%
    
features = [];
parfor sub=1:nsub
    
    
    infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'b_normeffdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',4); % Data must be in ROI montage
    
    % Remove bad badtrials
    p = []; fP=[];
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
    locLow   = find(fP==2);
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
A_for_Cov = covars{2};

aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
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
    parfor sub = 1:nsub
        
        infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'b_normeffdspmeeg');
        D = spm_eeg_load(infile);
        D = D.montage('switch',4); % Data must be in ROI montage
        
        % Remove bad badtrials
        g = D(:,:,:);
        g(:,:,D.badtrials)=[];
        
        % Remove source leakage using sysmetric Orthogonalisation
        y = reshape(g,D.nchannels,D.nsamples*size(g,3));  % generalise to Nsource,Ntime*Ntrials(good ones)
        y0 = ROInets.remove_source_leakage(y,'closest');
        
        % Filter to desired freq band
        y1 = ft_preproc_bandpassfilter(y0, D.fsample, freqbands{ii}, 4, 'but');
        y = reshape(y1,D.nchannels,D.nsamples,size(g,3));
        
        % Extract envelops
        Hen_lc_sep = [];
        for t=1:size(y,3)
            Hen_lc_sep(:,:,t) = hilbenv(squeeze(y(:,:,t)),1:D.nsamples,1,1);
        end
        
        Hen_lc_sep = reshape(Hen_lc_sep,D.nchannels,D.nsamples*size(g,3));
        %Hen_lc_sep = abs(Hen_lc_sep)
        
        % Calculate correlation
        cm = corr(resample(Hen_lc_sep',1,D.fsample));
        f(sub,:) = (cm(find(triu(cm,1))))';
        
    end
    features = [features f];
end
% Machine Learning
M   = features;
A_for_Cov = covars{2};
aM = M - A_for_Cov*pinv(A_for_Cov)*M; % Regress out
A = [aM group_num];
A = A(mri_num==1,:); % Remove subjects without T1 MRI
rng(1) % For reproducibility
accuracy = repeated_CV_matlab(A,kFolds,Nrun);

mean(accuracy)
std(accuracy)
