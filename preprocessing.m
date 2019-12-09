
%% Setup OSL including SPM

addpath('/imaging/dv01/Toolboxes/osl/osl-core') % Path to OSL directory
osl_startup
osl_check_installation


%% Definition of BIDS


bidspth = '/imaging/dv01/BIDS'; %Path to BIDS directory
BIDS   = spm_BIDS(bidspth);
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);
procpth = fullfile(bidspth,'derivatives','meg_derivatives'); % If want maxfiltered files


%% Pre-Processing- (Convertion, Preparation, Filtering(LP+HP), downsampling)


parfor sub = 1:length(subs)
    
    %tmp = spm_BIDS(BIDS,'metadata','ses','meg1','sub',subs{sub}); % This doesn't work for subjects with no MRI
    tmp = spm_jsonread(fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_meg.json']));
    event_file = spm_load(fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_events.tsv']));
    
    onset = event_file.onset*tmp.SamplingFrequency;
    offset = onset + event_file.duration*tmp.SamplingFrequency;
    
    % Converting %
    S = [];
    S.outfile = sprintf('/imaging/dv01/Processed/%s/spmeeg',subdir{sub});
    S.dataset = fullfile(procpth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_proc-tsss_meg.fif']);
    S.mode = 'epoched';
    S.channels = {'EOG', 'ECG', 'MEGMAG', 'MEGPLANAR'}; %EEG was removed
    S.checkboundary = 0;
    S.trl = [onset offset 0];
    try
        S.conditionlabels = event_file.stim_type;
    catch
        S.conditionlabels = event_file.trial_type;
    end
    D = spm_eeg_convert(S);
    
    % Set channel types and bad channels %
    S = [];
    S.D    = D;
    S.task = 'bidschantype';
    S.save = 1;
    S.filename = fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_channels.tsv']);
    D = spm_eeg_prep(S);
    
    D = chantype(D,indchantype(D,'megmag'),'MEG');
    D = chantype(D,indchantype(D,'megplanar'),'MEGPLANAR');
    D.save
    
    
    % Downsampling the data %
    S = [];
    S.D = D;
    S.method = 'resample';
    S.fsample_new = 500;
    D = spm_eeg_downsample(S);
    
    
    % High-pass filter %
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 0.5; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    
    % Low-pass filter %
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'low';
    S.freq = 150; % Cutoff frequency
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    
    
    % Epoch data % IF NEEDED
    time=2; %Epoch length in sec
    EpochLength = time * D.fsample; % in samples
    t = [1:EpochLength:(D.nsamples-EpochLength)]';
    nt = length(t);
    
    S = [];
    S.D       = D;
    S.trl = [t t+EpochLength-1 zeros(nt,1)];
    S.conditionlabels = repmat({'EYES_CLOSED'},1,nt);
    S.prefix='e';
    S.bc = 0;
    D = spm_eeg_epochs(S);
    
end

%% Automatic Artifact Detection using OSL


parfor sub = 1:length(subs)
    
    % OSL Artefact Detection
    infile = sprintf('/imaging/dv01/Processed/%s/effdspmeeg',subdir{sub}); % Can be epoched or continuous data
    D = spm_eeg_load(infile);
    D = osl_detect_artefacts(D,'modalities',unique(D.chantype(D.indchantype('MEGANY'))));
    % or D = osl_detect_artefacts(D,'badchannels',false); to only remove bad epochs
    D.save;
end


