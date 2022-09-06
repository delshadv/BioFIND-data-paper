% Step 1
%% Pre-Processing, Beamforming, ROI Extraction %%

% A multi-site, multi-participant magnetoencephalography resting-state
% dataset to study dementia: The BioFIND dataset
%
% This script contains preprocessing, co-registration, source-localisation
% and ROI extraction steps for reproducing BioFIND data paper results

% Note: Please make sure you are using the latest version of OSL and SPM
% within your OSL directory

% Henson R.N 2020, Vaghari D 2020

%% Define Paths ands variables

% Assumed you are currently in a directory which includes BioFIND data,
% OSL and code directories as described in readme.md

%restoredefaultpath
bwd = pwd;
wd  = fullfile(bwd,'code'); % Assumed you've already downloaded BioFind
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
procpth = fullfile(bidspth,'derivatives','meg_derivatives'); % If want maxfiltered files

%% Create Processed directory
% Please do all analysis in a separate directory from BIDS
% Here we call it "Processed"

processed_pth= fullfile(bwd,'Processed');

if ~exist(processed_pth,'dir')
    
    mkdir('Processed');
    cd ('Processed')
    for s=1:nsub
        mkdir(sprintf('sub-Sub%04d',s))
    end
else
end

cd (processed_pth)

%% PreProcess- Part 1 (Convert, Downsample, Filter)

parfor sub = 1:nsub
    
    % Read event & json file to extract desirable length of MEG Recordings
    
    tmp = spm_jsonread(fullfile(procpth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_proc-sss_meg.json']));
    event_file = spm_load(fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_events.tsv']));
    onset = (event_file.onset*tmp.SamplingFrequency)+1;
    
    % offset = onset + event_file.duration*tmp.SamplingFrequency;
    offset = (onset + 120 *tmp.SamplingFrequency)-1; % we put 120 seconds due to min length of raw data
    
    % Converting
    
    S = [];
    S.outfile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'spmeeg');
    S.dataset = fullfile(procpth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_proc-sss_meg.fif']);
    S.mode = 'epoched';
    S.channels = {'EOG', 'ECG', 'MEGMAG', 'MEGPLANAR'}; % EEG was removed
    S.checkboundary = 0;
    S.trl = [onset offset 0];
    try
        S.conditionlabels = event_file.stim_type;
    catch
        S.conditionlabels = event_file.trial_type;
    end
    D = spm_eeg_convert(S);
    
    % Set channel types and bad channels
    S = [];
    S.D    = D;
    S.task = 'bidschantype';
    S.save = 1;
    S.filename = fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_channels.tsv']);
    D = spm_eeg_prep(S);
    D = chantype(D,indchantype(D,'MEGMAG'),'MEGMAG');
    D = chantype(D,indchantype(D,'MEGGRADPLANAR'),'MEGPLANAR');
    D.save
    
    % Downsampling the data
    S = [];
    S.D = D;
    S.method = 'resample';
    S.fsample_new = 500;
    D = spm_eeg_downsample(S);
    delete(S.D)
    
    % High-pass filter
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 0.5; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    delete(S.D)
    
    
    % Low-pass filter
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'low';
    S.freq = 150; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    delete(S.D)
    
end


%% PreProcess- Part 2 - Epoching, OSL Artifacts detection

parfor sub = 1:nsub
    
    infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'ffdspmeeg');
    D = spm_eeg_load(infile);
    
    EpochLength = 2 * D.fsample; % in samples
    t = [1:EpochLength:(D.nsamples-EpochLength)]';
    nt = length(t);
    
    S = [];
    S.D = D;
    S.trl = [t t+EpochLength-1 zeros(nt,1)];
    S.conditionlabels = repmat({'EYES_CLOSED'},1,nt);
    S.prefix='e';
    S.bc = 0;
    D = spm_eeg_epochs(S);
    D.save;
    % OSL artifact detection
    D = osl_detect_artefacts(D,'modalities',unique(D.chantype(D.indchantype('MEGANY'))),'badchannels',false);
    D.save;
    
end


%% RUN Standard COREGISTRATION using SPM

% Copy T1w files to processed directory;

parfor sub = 1:nsub
    bidsT1 = fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','anat',[sprintf('sub-Sub%04d',sub) '_ses-meg1_T1w.nii.gz']);
    T1file = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),[sprintf('sub-Sub%04d',sub) '_ses-meg1_T1w.nii'])
    if exist(bidsT1,'file') & ~exist(T1file,'file')
        copyfile(bidsT1,[T1file '.gz']);
        gunzip([T1file '.gz'])
        delete([T1file '.gz'])
    end
end

% Coregisteration

UseHPs = 1; % Use headpoints

parfor sub = 1:nsub
    
    
    infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'effdspmeeg');
    D = spm_eeg_load(infile);
    D = D.montage('switch',0);
    
    % Remove headpoints around face (y>0 & z<0), particularly important if MRI de-faced
    fid = D.fiducials;
    fid.pnt(find(fid.pnt(:,2)>0 & fid.pnt(:,3)<0),:)=[];
    D=fiducials(D,fid); D.save
    
    %D = rmfield(D,'inv');
    
    T1file = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),[sprintf('sub-Sub%04d',sub) '_ses-meg1_T1w.nii']);
    if exist(T1file,'file')
        D = spm_eeg_inv_mesh_ui(D,1,T1file,2);
        
        V = spm_vol(T1file);
        fids = spm_jsonread(fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','anat',[sprintf('sub-Sub%04d',sub) '_ses-meg1_T1w.json']));
        
        mrifid = [];
        megfid = D.fiducials;
        mrifid.fid.label = {'NAS';'LPA';'RPA'};
        nas = V.mat*[fids.AnatomicalLandmarkCoordinates.NAS; 1]; nas = nasion(1:3)';
        lpa    = V.mat*[fids.AnatomicalLandmarkCoordinates.LPA; 1];    lpa    = lpa(1:3)';
        rpa    = V.mat*[fids.AnatomicalLandmarkCoordinates.RPA; 1];    rpa    = rpa(1:3)';
        mrifid.fid.pnt = [nas; lpa; rpa];
        mrifid.pnt = D.inv{1}.mesh.fid.pnt;
        D = spm_eeg_inv_datareg_ui(D, 1, megfid, mrifid, UseHPs);
        if UseHPs, D.inv{1}.comment = 'With headpoints'; else, D.inv{1}.comment = 'SPM fids only'; end
        D.save;
        
    else
        D = spm_eeg_inv_mesh_ui(D,1,1,2);
        mrifid = [];
        for f=1:3
            mrifid.fid.label{f} = D.inv{1}.mesh.fid.fid.label{f};
            mrifid.fid.pnt(f,:) = D.inv{1}.mesh.fid.fid.pnt(f,:);
        end
        megfid = D.fiducials;
        megfid.fid.label{find(strcmp([megfid.fid.label],'LPA'))} = 'lpa';
        megfid.fid.label{find(strcmp([megfid.fid.label],'RPA'))} = 'rpa';
        megfid.fid.label{find(strcmp([megfid.fid.label],'NAS'))} = 'nas';
        D = spm_eeg_inv_datareg_ui(D, 1, megfid, mrifid, 0);
        D.save;
    end
    
    D.inv{1}.forward(1).voltype = 'Single Shell';
    D = spm_eeg_inv_forward(D);
    D.save;
    
    
end

%% Beamforming and Extract ROIs using OSL

p           = parcellation(fullfile('fMRI_parcellation_ds8mm.nii.gz'));
mni_coords  = p.template_coordinates;

parfor sub = 1:nsub
    
    infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'effdspmeeg');
    D = spm_eeg_load(infile);
    D = osl_filter(D,[0.5 48]);
    
    S=[];
    S.D = D;
    S.datatype ='neuromag';
    S.do_plots = false;
    [D,~] = osl_normalise_sensor_data(S);
    % Run Beamforming
    
    S = struct;
    S.modalities        = {'MEG','MEGPLANAR'};
    S.fuse              = 'all';
    S.timespan          = [0 Inf];
    S.pca_order         = 50;
    S.type              = 'Scalar';
    S.inverse_method    = 'beamform';
    S.prefix            = 'b_norm';
    
    D = osl_inverse_model(D,mni_coords,S);
    
    % Select montage
    D = D.montage('switch',2); % Unweighted inverse model
    
    % Extract ROI (Region Of Interest)
    D = ROInets.get_node_tcs(D,p.to_matrix(p.binarize),'pca');
    
    % Save the data
    D.save;
    
    
end

%% Table 1 - pre-requirements

participants = spm_load(fullfile(wd,'participants-imputed.tsv'));
group_num    = grp2idx(participants.group);
mri_num      = grp2idx(participants.sImaging);

%% Number of bad trials

number_bad_epochs = nan(length(subs),1);
for sub = 1:length(subs)
    
    infile = fullfile(processed_pth,sprintf('sub-Sub%04d',sub),'effdspmeeg');
    D = spm_eeg_load(infile);
    
    number_bad_epochs(sub,:) = length(D.badtrials);
    
end

nanmean(number_bad_epochs(group_num==1)) % Control
nanstd(number_bad_epochs(group_num==1)) % Control
nanmean(number_bad_epochs(group_num==2)) % MCI
nanstd(number_bad_epochs(group_num==2)) % MCI

[h,p,ci,stats] = ttest2(number_bad_epochs(group_num==1),number_bad_epochs(group_num==2))

%% Recording duration

rec_dur = nan(length(subs),1);
for sub = 1:length(subs)
    
    % Read event & json file to calculate offset
    
    tmp = spm_jsonread(fullfile(procpth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_proc-sss_meg.json']));
    event_file = spm_load(fullfile(bidspth,sprintf('sub-Sub%04d',sub),'ses-meg1','meg',[sprintf('sub-Sub%04d',sub) '_ses-meg1_task-Rest_events.tsv']));
    rec_dur(sub) = event_file.duration;
    
end
[p,h,stats] = ranksum(rec_dur(group_num==1),rec_dur(group_num==2))
nanmedian(rec_dur(group_num==1)) % Control
iqr(rec_dur(group_num==1)) % Control
nanmedian(rec_dur(group_num==2)) % MCI
iqr(rec_dur(group_num==2)) % MCI


%% Pre Task

PreT  = participants.Pre_task;
[p,h,stats] = ranksum(PreT(group_num==1),PreT(group_num==2))
nanmedian(PreT(group_num==1)) % Control
iqr(PreT(group_num==1)) % Control
nanmedian(PreT(group_num==2)) % MCI
iqr(PreT(group_num==2)) % MCI

%% Mean and SD of head translation

move1 = participants.Move1; %Mean of head trans
move2 = participants.Move2; %STD of head trans

nanmean(move1(group_num==1)) % Control
nanstd(move1(group_num==1)) % Control
nanmean(move1(group_num==2)) % MCI
nanstd(move1(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(move1(group_num==1),move1(group_num==2))

nanmean(move2(group_num==1)) % Control
nanstd(move2(group_num==1)) % Control
nanmean(move2(group_num==2)) % MCI
nanstd(move2(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(move2(group_num==1),move2(group_num==2))

%% Recording hour and recording year

tday = participants.Recording_time;
tyear = participants.Recording_year;

nanmean(tday(group_num==1)) % Control
nanstd(tday(group_num==1)) % Control
nanmean(tday(group_num==2)) % MCI
nanstd(tday(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(tday(group_num==1),tday(group_num==2))

nanmean(tyear(group_num==1)) % Control
nanstd(tyear(group_num==1)) % Control
nanmean(tyear(group_num==2)) % MCI
nanstd(tyear(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(tyear(group_num==1),tyear(group_num==2))

%% MMSE, Education, Year

edu = participants.Edu_years;
age = participants.age;
MMSE = participants.MMSE;

nanmean(edu(group_num==1)) % Control
nanstd(edu(group_num==1)) % Control
nanmean(edu(group_num==2)) % MCI
nanstd(edu(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(edu(group_num==1),edu(group_num==2))

nanmean(age(group_num==1)) % Control
nanstd(age(group_num==1)) % Control
nanmean(age(group_num==2)) % MCI
nanstd(age(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(age(group_num==1),age(group_num==2))

nanmean(MMSE(group_num==1)) % Control
nanstd(MMSE(group_num==1)) % Control
nanmean(MMSE(group_num==2)) % MCI
nanstd(MMSE(group_num==2)) % MCI
[h,p,ci,stats] = ttest2(MMSE(group_num==1),MMSE(group_num==2))

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
