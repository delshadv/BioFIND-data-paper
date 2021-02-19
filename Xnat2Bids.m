

cd S: % 
mkdir ('BioFIND','MCIControls')
mkdir ('BioFIND','TravelBrains')

nsub = 324;

bwd = pwd;
path_xnat = '/imaging/dv01/Github/Xnat/projects/biofind/'; % will be specified later: bwd/projects/biofind
path_bids = '/imaging/dv01/Github/BioFIND-BioFIND/'; % will be specified later : bwd/BioFIND
%% MCIControls

mkdir (fullfile(path_bids,'BioFIND','MCIControls','derivatives','meg_derivatives'))
mkdir (fullfile(path_bids,'BioFIND','MCIControls','sub-emptyroom'))

% Raw data
xnatfile = fullfile(path_xnat,'Resources','NO LABEL');
bidsfile = fullfile(path_bids,'BioFIND','MCIControls');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% create directories for raw data

for s = 1:nsub
    
    mkdir (fullfile(path_bids,'BioFIND','MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','anat'));
    mkdir (fullfile(path_bids,'BioFIND','MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','meg'));
    
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

% create directories for maxfiltered data

for s = 1:nsub
    
    mkdir (fullfile(path_bids,'BioFIND','MCIControls','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-meg1','meg'));

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

mkdir (fullfile(path_bids,'BioFIND','TravelBrains','derivatives','meg_derivatives'))
nsub = 7;

% Raw data
xnatfile = fullfile(path_xnat,'Resources','TravelBrains_Participants');
bidsfile = fullfile(path_bids,'BioFIND','TravelBrains');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT


for s = 1:nsub
    
    mkdir (fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','anat'));
    mkdir (fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','meg'));
    mkdir (fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','anat'));
    mkdir (fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','meg'));
    
end



for s = 1:nsub

    % move anat dir for CBU
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-TravelBrains%04d',s),'experiments',sprintf('sub-TravelBrains%04d_ses-megCBU',s),'scans','anat','resources','NIFTI');
    bidsfile = fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','anat');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
    % move anat dir for CTB
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-TravelBrains%04d',s),'experiments',sprintf('sub-TravelBrains%04d_ses-megCTB',s),'scans','anat','resources','NIFTI');
    bidsfile = fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','anat');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
    % move meg dir for CBU
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-TravelBrains%04d',s),'experiments',sprintf('sub-TravelBrains%04d_ses-megCBU',s),'scans','meg','resources','FIF');
    bidsfile = fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
    % move meg dir for CTB
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-TravelBrains%04d',s),'experiments',sprintf('sub-TravelBrains%04d_ses-megCTB',s),'scans','meg','resources','FIF');
    bidsfile = fullfile(path_bids,'BioFIND','TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
  
end

% Derivatives (MaxFilter)

xnatfile = fullfile(path_xnat,'Resources','Derivatives_TravelBrains');
bidsfile = fullfile(path_bids,'BioFIND','Travelbrains','derivatives','meg_derivatives');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% create directories for maxfiltered data

for s = 1:nsub
    
    mkdir (fullfile(path_bids,'BioFIND','Travelbrains','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-meg1','meg'));

end

% copy related content for above

for s = 1:nsub
        
    % move meg dir
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-TravelBrains%04d',s),'experiments',sprintf('sub-TravelBrains%04d_ses-megCBU_derivatives',s),'scans','meg','resources','FIF');
    bidsfile = fullfile(path_bids,'BioFIND','Travelbrains','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-megCBU','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
     % move meg dir
    xnatfile = fullfile(path_xnat,'subjects',sprintf('sub-TravelBrains%04d',s),'experiments',sprintf('sub-TravelBrains%04d_ses-megCTB_derivatives',s),'scans','meg','resources','FIF');
    bidsfile = fullfile(path_bids,'BioFIND','Travelbrains','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-megCTB','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT

end
