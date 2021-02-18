

cd S: % 
mkdir ('BioFIND','MCIControls')
mkdir ('BioFIND','TravelBrains')
path_xnat = '/imaging/dv01/Github/Xnat/projects/biofind/';
path_bids = '/imaging/dv01/Github/BioFIND-BioFIND/';
%% MCIControls

cd BioFIND/MCIControls % all need to be windows based format e.g. / => \
mkdir ('derivatives','meg_derivatives')
mkdir ('sub-emptyroom')

% Raw data
xnatfile = fullfile(path_xnat,sprintf('Resources/NO LABEL/'));
bidsfile = fullfile(path_bids,sprintf('BioFIND/MCIControls/'));

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% create directories for raw data

for s = 1:324
    
    mkdir (fullfile(sprintf('sub-Sub%04d',s),'ses-meg1','anat'));
    mkdir (fullfile(sprintf('sub-Sub%04d',s),'ses-meg1','meg'));
    
end

% copy related content for above


for s = 1:324

    % move anat dir
    xnatfile = fullfile(path_xnat,sprintf('subjects/sub-Sub%04d/experiments/sub-Sub%04d_ses-meg1/scans/anat/resources/NIFTI/',s,s));
    bidsfile = fullfile(path_bids,sprintf('BioFIND/MCIControls/sub-Sub%04d/ses-meg1/anat/',s));

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
    % move meg dir
    xnatfile = fullfile(path_xnat,sprintf('subjects/sub-Sub%04d/experiments/sub-Sub%04d_ses-meg1/scans/meg/resources/FIF/',s,s));
    bidsfile = fullfile(path_bids,sprintf('BioFIND/MCIControls/sub-Sub%04d/ses-meg1/meg/',s));

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT

end

% Derivatives (MaxFilter)

xnatfile = fullfile(path_xnat,sprintf('Resources/Derivatives/'));
bidsfile = fullfile(path_bids,sprintf('BioFIND/MCIControls/derivatives/meg_derivatives/'));

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

cd derivatives/meg_derivatives
% create directories for raw data

for s = 1:324
    
    mkdir (fullfile(sprintf('sub-Sub%04d',s),'ses-meg1','meg'));

end

% copy related content for above

for s = 1:324
        
    % move meg dir
    xnatfile = fullfile(path_xnat,sprintf('subjects/sub-Sub%04d/experiments/sub-Sub%04d_ses-meg1_derivatives/scans/meg/resources/FIF/',s,s));
    bidsfile = fullfile(path_bids,sprintf('BioFIND/MCIControls/derivatives/meg_derivatives/sub-Sub%04d/ses-meg1/meg/',s));

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT

end


% Empty Rooms
% will bee specified later

%% TravelBrains
% will bee specified later





