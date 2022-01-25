
% Paths below for FIFFACCESS toolbox (http://kimmouutela.yolasite.com/meg-pd.php)
addpath /imaging/local/meg_misc/
addpath /neuro/meg_pd_1.2/
addpath /neuro/bin/util

clear

bwd = '/imaging/henson/users/dv01/BioFIND/MCIControls';

hwd = '/imaging/henson/users/rh01/Projects/AD_MCI/Biofind/BIDS'

dosubs = [1:324]+1; %+1 because of TSV header below!

Nsub = length(dosubs); 

[subid,group,site,sex,age,MMSE,sMRI,Converter,Recording_year,Recording_time,Edu_years,Move1,Move2,PreTask] = textread(fullfile(bwd,'participants.tsv'),'%s%s%s%s%s%s%s%s%s%s%s%s%s%s');

cal_num = textread('CTBcal.txt','%d'); % CTB site has two different calibration files

cd(hwd)

owd = fullfile(hwd,'derivatives','meg_derivatives'); 

sites = {'CTB';'CBU'};

groups = {'control','patient'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Maxfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxfstr = '!/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2.12 '
%%tSSSstr = ' -st 10 -corr 0.98'; 
%tSSSstr = ' -st 10 -corr 0.9'; 
tSSSstr = '';
dsstr = ''; % dsstr = ' -ds 4';  % downsample to 250Hz - can cause MF crash?
incEEG = 0;
OverWrite = 1;
Nbuf = 600; % Number of buffers per session 

autobad_flag = 1;

% Set-up directories because doesn't always work in parfor loop 
for subnum = dosubs
    sub_dir{subnum} = fullfile(owd,subid{subnum});
    try eval(sprintf('!mkdir %s',sub_dir{subnum})); end
    out_dir{subnum} = fullfile(sub_dir,'ses-meg1');
    try eval(sprintf('!mkdir %s',out_dir{subnum})); end
    out_dir{subnum} = fullfile(sub_dir,'ses-meg1','meg');
    try eval(sprintf('!mkdir %s',out_dir{subnum})); end
end


parfor subnum = dosubs

    cal_ver  = '';
    
    if strcmp(site{subnum},'CTB')
%        cal_sub = cal_sub+1;
        cal_ver = sprintf('_%d',cal_num(subnum-1))
    end
        
    base_name = sprintf('%s_ses-meg1_task-Rest',subid{subnum});
    raw_file = fullfile(bwd,subid{subnum},'ses-meg1','meg',[base_name '_meg.fif'])
   
    basestr = sprintf(' -ctc %s/ct_sparse_%s.fif -cal %s/sss_cal_%s%s.dat',hwd,site{subnum},hwd,site{subnum},cal_ver);
    basestr = [basestr ' -linefreq 50 -hpisubt amp -force'];
    
    if ~exist(raw_file,'file')
        error(raw_file,' does not exist')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit sphere (since better than MaxFilter does)
    
    if exist(fullfile(sub_dir{subnum},'fittmp.txt'),'file'); delete(fullfile(sub_dir{subnum},'fittmp.txt'),'file'); end
    [orig,rad,fit] = meg_fit_sphere(raw_file,sub_dir{subnum},[base_name '_hpi.txt'],incEEG);
    delete(fullfile(sub_dir{subnum},'fittmp.txt'))
    
    origstr = sprintf(' -origin %d %d %d -frame head',orig(1),orig(2),orig(3))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Bad channel detection (this email says important if doing tSSS later https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=NEUROMEG;d3f363f3.1205)
    
    badfile = fullfile(sub_dir{subnum},[base_name '_bad.txt']); %sprintf('bad_autobad%d.txt',autobad_flag)); 
    
    if autobad_flag
        outfile = fullfile(sub_dir{subnum},'bad'); logfile = fullfile(sub_dir{subnum},'bad.log');
        badstr  = sprintf(' -autobad %d -badlimit %d',1800,7); % 1800s is 30mins - ie enough for all do_sessions?
        
        if ~exist(logfile,'file') | OverWrite
            filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);
            posfile = fullfile(sub_dir{subnum},[base_name '_headpos.txt']);
            compstr = sprintf(' -headpos -hpistep 10 -hp %s',posfile);
            finstr = [maxfstr filestr origstr basestr badstr compstr sprintf(' -v | tee %s.log',outfile)]
            rik_eval(finstr);
            delete(sprintf('%s.fif',outfile));
        end
        
        % Pull out bad channels from logfile:
        delete(badfile); rik_eval(sprintf('!echo '' '' > %s', badfile));
        rik_eval(sprintf('!cat %s.log | sed -n -e ''/Detected/p'' -e ''/Static/p'' | cut -f 5- -d '' '' >> %s',outfile,badfile));
        delete(logfile)
    else
         error('Manual bad channel files not available yet')
    end
    
    tmp=dlmread(badfile,' '); % Nbuf = size(tmp,1)-1; Doesn't work if subset of channels bad for short time
    tmp=reshape(tmp,1,prod(size(tmp)));
    tmp=tmp(tmp>0); % Omit zeros (padded by dlmread):
    
    % Mark bad based on threshold (currently ~5% of buffers (assuming 500 buffers)):
    [frq,allbad] = hist(tmp,unique(tmp));
    badchans = allbad(frq>0.05*Nbuf);
    if isempty(badchans) badstr = '';
    else
        badstr = sprintf(' -bad %s',num2str(badchans))
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. tSSS and trans to first file
    
    outfile = fullfile(sub_dir{subnum},[base_name '_proc-sss']); %sprintf('sss_autobad%d',autobad_flag));
    if ~exist(sprintf('%s.fif',outfile),'file') | ~exist([outfile '.log'],'file') | OverWrite
        transtr = '';
        skipstr = '';
        compstr = ''
                
        filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);
        finstr = [maxfstr filestr basestr badstr tSSSstr skipstr compstr origstr transtr dsstr sprintf(' -v | tee %s.log',outfile)]
        rik_eval(finstr);
    end
end

return
