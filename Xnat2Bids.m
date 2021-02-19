

cd S: % 
mkdir ('BioFIND','MCIControls')
mkdir ('BioFIND','TravelBrains')

nsub = 324;

bwd = pwd;
path_xnat = '/imaging/dv01/Github/Xnat/projects/biofind/'; % will be specified later: bwd/projects/biofind
path_bids = '/imaging/dv01/Github/BioFIND-BioFIND/'; % will be specified later : bwd/BioFINF
%% MCIControls

mkdir ('BioFIND','MCIControls','derivatives','meg_derivatives')
mkdir ('BioFIND','MCIControls''sub-emptyroom')

% Raw data
xnatfile = fullfile(path_xnat,'Resources','NO LABEL');
bidsfile = fullfile(path_bids,'BioFIND','MCIControls');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% create directories for raw data

for s = 1:nsub
    
    mkdir (fullfile(sprintf('sub-Sub%04d',s),'ses-meg1','anat'));
    mkdir (fullfile(sprintf('sub-Sub%04d',s),'ses-meg1','meg'));
    
end

% copy related content for above


for s = 1:nsub

    % move anat dir
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-Sub%04d',s),'experiments',sprintf('sub-Sub%04d_ses-meg1',s),'scans','anat','resources','NIFTI');
    bidsfile = fullfile(path_bids,'BioFIND','MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','anat');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
    % move meg dir
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-Sub%04d',s),'experiments',sprintf('sub-Sub%04d_ses-meg1',s),'scans','meg','resources','FIF');
    bidsfile = fullfile(path_bids,'BioFIND','MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
end

% Derivatives (MaxFilter)

xnatfile = fullfile(path_xnat,'Resources','Derivatives');
bidsfile = fullfile(path_bids,'BioFIND','MCIControls','derivatives','meg_derivatives');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

cd derivatives/meg_derivatives
% create directories for raw data

for s = 1:nsub
    
    mkdir (fullfile('derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-meg1','meg'));

end

% copy related content for above

for s = 1:nsub
        
    % move meg dir
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-Sub%04d',s),'experiments',sprintf('sub-Sub%04d_ses-meg1_derivatives',s),'scans','meg','resources','FIF');
    bidsfile = fullfile(path_bids,'BioFIND','MCIControls','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-meg1','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT

end


% Empty Rooms
% will bee specified later

%% TravelBrains
% will bee specified later




