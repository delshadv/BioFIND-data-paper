
%% Pre-Processing, Beamforming, ROI Extraction %%

% A multi-site, multi-participant magnetoencephalography resting-state dataset
% to study dementia: The BioFIND dataset

% This script contains preprocessing, co-registration, Source-localisation
% and ROI extraction steps for reproducing BioFIND data paper results

% Note: Please make sure you are using the latest version of OSL and SPM
% within your OSL directory

% Henson R.N 2020, Vaghari D 2020

%% Setup OSL

addpath('/imaging/dv01/Toolboxes/osl/osl-core') % Path to OSL directory
osl_startup
osl_check_installation

%% Adding BIDS Paths, defining BIDS variables

% Please do all analysis in a sperate directory from your BIDS data. Here, we
% call it processed_pth
processed_pth= '/imaging/dv01/Processed_sss_new';
bidspth = '/imaging/dv01/BioFIND/MCIControls'; %BIDS Path;
BIDS   = spm_BIDS(bidspth); % (If it does not work with OSL's SPM, so copied last version of spm_BIDS)
subs   = spm_BIDS(BIDS,'subjects', 'task', 'Rest');
nsub   = numel(subs);
subdir = cellfun(@(s) ['sub-' s], subs, 'UniformOutput',false);
procpth = fullfile(bidspth,'derivatives','meg_derivatives'); % If want maxfiltered files

%% PreProcess- Part 1 (Convert, Downsample, Filter)

parfor sub = 1:length(subs)
    
    % Read event & json file to extract desirable length of MEG Recordings
    
    tmp = spm_jsonread(fullfile(procpth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_proc-sss_meg.json']));
    event_file = spm_load(fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_events.tsv']));
    onset = (event_file.onset*tmp.SamplingFrequency)+1;
    
    % offset = onset + event_file.duration*tmp.SamplingFrequency;
    offset = (onset + 120 *tmp.SamplingFrequency)-1; % we put 120 seconds due to min length of raw data
    
    % Converting
    
    S = [];
    S.outfile = fullfile(processed_pth,subdir{sub},'spmeeg');
    S.dataset = fullfile(procpth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_proc-sss_meg.fif']);
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
    
    % Set channel types and bad channels (Does not work with OSL's SPM, so copied last version of spm_eeg_prep)
    S = [];
    S.D    = D;
    S.task = 'bidschantype';
    S.save = 1;
    S.filename = fullfile(bidspth,subdir{sub},'ses-meg1','meg',[subdir{sub} '_ses-meg1_task-Rest_channels.tsv']);
    D = spm_eeg_prep(S);
    D = chantype(D,indchantype(D,'megmag'),'MEGMAG');
    D = chantype(D,indchantype(D,'megplanar'),'MEGPLANAR');
    D.save
    
    % Downsampling the data %
    S = [];
    S.D = D;
    S.method = 'resample';
    S.fsample_new = 500;
    D = spm_eeg_downsample(S);
    delete(S.D)
    
    % High-pass filter %
    S = [];
    S.D = D;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 0.5; % Cutoff frequency    %%%% should need to use lower freq to have freq=1?
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

parfor sub = 1:length(subs)
    
    infile = fullfile(processed_pth,subdir{sub},'ffdspmeeg');
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

parfor sub=1:length(subs)
    bidsT1 = fullfile(bidspth,subdir{sub},'ses-meg1','anat',[subdir{sub} '_ses-meg1_T1w.nii.gz']);
    T1file = fullfile(processed_pth,subdir{sub},[subdir{sub} '_ses-meg1_T1w.nii'])
    if exist(bidsT1,'file') & ~exist(T1file,'file')
        copyfile(bidsT1,[T1file '.gz']);
        gunzip([T1file '.gz'])
        delete([T1file '.gz'])
    end
end

% Coregisteration

parfor sub=1:length(subs)
    
    infile = fullfile(processed_pth,subdir{sub},'effdspmeeg');
    D = spm_eeg_load(infile);
    
    T1file = fullfile(processed_pth,subdir{sub},[subdir{sub} '_ses-meg1_T1w.nii']);
    if exist(T1file,'file')
        D = spm_eeg_inv_mesh_ui(D,1,T1file,2);
        
        V = spm_vol(T1file);
        fids = spm_jsonread(fullfile(bidspth,subdir{sub},'ses-meg1','anat',[subdir{sub} '_ses-meg1_T1w.json']));
        
        mrifid = [];
        megfid = D.fiducials;
        mrifid.fid.label = {'Nasion';'LPA';'RPA'};
        nasion = V.mat*[fids.AnatomicalLandmarkCoordinates.Nasion; 1]; nasion = nasion(1:3)';
        lpa    = V.mat*[fids.AnatomicalLandmarkCoordinates.LPA; 1];    lpa    = lpa(1:3)';
        rpa    = V.mat*[fids.AnatomicalLandmarkCoordinates.RPA; 1];    rpa    = rpa(1:3)';
        mrifid.fid.pnt = [nasion; lpa; rpa];
        D = spm_eeg_inv_datareg_ui(D, 1, megfid, mrifid, 0);
        D.inv{1}.comment = 'SPM fids only';
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
        megfid.fid.label{find(strcmp([megfid.fid.label],'Nasion'))} = 'nas';
        D = spm_eeg_inv_datareg_ui(D, 1, megfid, mrifid, 0);
        D.save;
    end
    
    D.inv{1}.forward(1).voltype = 'Single Shell';
    D = spm_eeg_inv_forward(D);
    fid = D.fiducials;
    fid.pnt(find(fid.pnt(:,2)>0 & fid.pnt(:,3)<0),:)=[];
    D=fiducials(D,fid);
    D.save
    
end

%% Beamforming and Extract ROIs using OSL

p           = parcellation(fullfile('fMRI_parcellation_ds8mm.nii.gz'));
mni_coords  = p.template_coordinates;

parfor sub = 1:length(subs) 
    
    infile = fullfile(processed_pth,subdir{sub},'effdspmeeg');
    D = spm_eeg_load(infile);
    
    % Run Beamforming

    S = struct;
    S.modalities        = {'MEG','MEGPLANAR'};
    S.fuse              = 'all';
    S.timespan          = [0 Inf];
    S.pca_order         = 50;
    S.type              = 'Scalar';
    S.inverse_method    = 'beamform';
    S.prefix            = 'b';
    
    D = osl_inverse_model(D,mni_coords,S);

    % Select montage  
    D = D.montage('switch',1);
    
    % Extract ROI (Region Of Interest)
    D = ROInets.get_node_tcs(D,p.to_matrix(p.binarize),'pca');
          
    % Save the data
    D.save;
    
end

