

cd 'S:\BioFIND - BioFIND' % 
mkdir ('BioFIND','MCIControls')
mkdir ('BioFIND','TravelBrains')

nsub = 324;

bwd = pwd;
path_xnat = 'S:\BioFIND - BioFIND\images'; 
path_bids = 'S:\BioFIND - BioFIND\BioFIND';
%% MCIControls

mkdir (fullfile(path_bids,'MCIControls','derivatives','meg_derivatives'))
mkdir (fullfile(path_bids,'MCIControls','sub-emptyroom'))

% Raw data
xnatfile = fullfile(path_xnat,'Resources','44390');
bidsfile = fullfile(path_bids,'MCIControls');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% create directories for raw data

for s = 1:nsub
    
    mkdir (fullfile(path_bids,'MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','anat'));
    mkdir (fullfile(path_bids,'MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','meg'));
    
end

% copy related content for above


for s = 1:nsub
    
    try
    % move anat dir
    xnatfile = fullfile(path_xnat,sprintf('sub-Sub%04d',s),sprintf('sub-Sub%04d_ses-meg1',s),'anat');
    bidsfile = fullfile(path_bids,'MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','anat');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    catch
        warning('NO MRI')
    end
    % move meg dir
    xnatfile = fullfile(path_xnat,sprintf('sub-Sub%04d',s),sprintf('sub-Sub%04d_ses-meg1',s),'meg');
    bidsfile = fullfile(path_bids,'MCIControls',sprintf('sub-Sub%04d',s),'ses-meg1','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
end

% Derivatives (MaxFilter)

xnatfile = fullfile(path_xnat,'Resources','Derivatives');
bidsfile = fullfile(path_bids,'MCIControls','derivatives','meg_derivatives');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% create directories for maxfiltered data

for s = 1:nsub
    
    mkdir (fullfile(path_bids,'MCIControls','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-meg1','meg'));

end

% copy related content for above

for s = 1:nsub
        
    % move meg dir
    xnatfile = fullfile(path_xnat,sprintf('sub-Sub%04d',s),sprintf('sub-Sub%04d_ses-meg1_derivatives',s),'meg');
    bidsfile = fullfile(path_bids,'MCIControls','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-meg1','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT

end


% Empty Rooms

for s = 2009:2019
    try
    % move meg dir
    xnatfile = fullfile(path_xnat,'sub-emptyroom',sprintf('sub-emptyroom_ses-%04dCBU',s),'meg');
    bidsfile = fullfile(path_bids,'MCIControls','sub-emptyroom',sprintf('ses-%04dCBU',s),'meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    catch
        warning('NO SITE FOR THIS YEAR')
    end
    try
    % move meg dir
    xnatfile = fullfile(path_xnat,'sub-emptyroom',sprintf('sub-emptyroom_ses-%04dCTB',s),'meg');
    bidsfile = fullfile(path_bids,'MCIControls','sub-emptyroom',sprintf('ses-%04dCTB',s),'meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    catch
        warning('NO SITE FOR THIS YEAR')
    end
end

%% TravelBrains

mkdir (fullfile(path_bids,'TravelBrains','derivatives','meg_derivatives'))
nsub = 7;

% Raw data
xnatfile = fullfile(path_xnat,'Resources','TravelBrains_Participants');
bidsfile = fullfile(path_bids,'TravelBrains');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT


for s = 1:nsub
    
    mkdir (fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','anat'));
    mkdir (fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','meg'));
    mkdir (fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','anat'));
    mkdir (fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','meg'));
    
end



for s = 1:nsub
    
    try
    % move anat dir for CBU
    xnatfile = fullfile(path_xnat,sprintf('sub-TravelBrains%04d',s),sprintf('sub-TravelBrains%04d_ses-megCBU',s),'anat');
    bidsfile = fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','anat');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    catch
        warning('NO MRI')
    end
    
    try
    % move anat dir for CTB
    xnatfile = fullfile(path_xnat,sprintf('sub-TravelBrains%04d',s),sprintf('sub-TravelBrains%04d_ses-megCTB',s),'anat');
    bidsfile = fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','anat');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    catch
        warning('NO MRI')
    end
    % move meg dir for CBU
    xnatfile = fullfile(path_xnat,sprintf('sub-TravelBrains%04d',s),sprintf('sub-TravelBrains%04d_ses-megCBU',s),'meg');
    bidsfile = fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCBU','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
    % move meg dir for CTB
    xnatfile = fullfile(path_xnat,sprintf('sub-TravelBrains%04d',s),sprintf('sub-TravelBrains%04d_ses-megCTB',s),'meg');
    bidsfile = fullfile(path_bids,'TravelBrains',sprintf('sub-Sub%04d',s),'ses-megCTB','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
  
end

% Derivatives (MaxFilter)

xnatfile = fullfile(path_xnat,'Resources','Derivatives_TravelBrains');
bidsfile = fullfile(path_bids,'Travelbrains','derivatives','meg_derivatives');

movefile(fullfile(xnatfile,'*'),bidsfile); % could be copyfile if needed to keep XNAT

% copy the content, No need to create empty folders

for s = 1:nsub
        
    % move meg dir
    xnatfile = fullfile(path_xnat,sprintf('sub-TravelBrains%04d',s),sprintf('sub-TravelBrains%04d_ses-megCBU_derivatives',s),'meg');
    bidsfile = fullfile(path_bids,'Travelbrains','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-megCBU','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT
    
     % move meg dir
    xnatfile = fullfile(path_xnat,sprintf('sub-TravelBrains%04d',s),sprintf('sub-TravelBrains%04d_ses-megCTB_derivatives',s),'meg');
    bidsfile = fullfile(path_bids,'Travelbrains','derivatives','meg_derivatives',sprintf('sub-Sub%04d',s),'ses-megCTB','meg');

    movefile(fullfile(xnatfile,'sub*'),bidsfile); % could be copyfile if needed to keep XNAT

end
