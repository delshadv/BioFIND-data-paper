
% Just re-run for 3 Cambridge replacements (but could be re-run on all

% Notes:   Mvcomp turned off for all because fails for some subjects -
% could re-run just a mvcomp option to estimate motion even if not corrected?

addpath /imaging/local/meg_misc/
addpath /neuro/meg_pd_1.2/
addpath /neuro/bin/util

clear

bwd = '/imaging/dv01/BioFIND/MCIControls';

hwd = '/imaging/rh01/Projects/AD_MCI/Biofind/BIDS'

%dosubs = [166:168]+1; % 3 CBU replacements of OHBA D&F, +1 because of TSV header below
dosubs = [1:324]+1; % new CTB ones, +1 because of TSV header below

Nsub = length(dosubs); 

[subid,group,site,sex,age,MMSE,sMRI,Converter,Recording_year,Recording_time,Edu_years,Move1,Move2] = textread(fullfile(bwd,'participants.tsv'),'%s%s%s%s%s%s%s%s%s%s%s%s%s');

cd(hwd)

owd = fullfile(hwd,'derivatives','meg_derivatives'); % should this include "derivatives" directory before "meg_derivatives"?

sites = {'CTB';'CBU'};

groups = {'control','patient'};

CTB = load('/imaging/rh01/Projects/AD_MCI/Biofind/Biofind_Info_CTB_new.mat');

try CTB = CTB.CTB; end

%% Do just once - renaming files from CTB
% 
% for s = 1:length(dosubs)
%     ss = dosubs(s)-1;
%     if (ss >= 244 & ss <= 276) | ss == 324 % new controls and last patient! (see "Biofind Batches 1 & 2.xls")
%         CTB.Biofind.Madrid.calfiles(ss) = 2;
%     elseif ss >= 277 & ss <= 323
%         CTB.Biofind.Madrid.calfiles(ss) = 1;
%     end
% %    if strcmp(group{s},'control'), cgrp = 'CN'; elseif strcmp(group{s},'patient'), cgrp = 'MCI'; else error('!'); end
% %    eval(sprintf('!mv /imaging/rh01/Projects/AD_MCI/Biofind/NewCTB/MEG_raw/MAD-%s-%04d_resting.fif %s',cgrp,s-1,
% end
% 
% save('/imaging/rh01/Projects/AD_MCI/Biofind/Biofind_Info_CTB_new.mat','CTB');


%FailedMvComp = {'CN015','CN019','CN020','MCI003','MCI004','MCI005','MCI006','MCI007','MCI020','MCI021','MCI022'};
%CBUFailedMvComp = [87]; % Added later. Note done on subject numbers not original names. Subject 87 had missed cHPI for first lot of scan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Maxfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxfstr = '!/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2.12 '
%%tSSSstr = ' -st 10 -corr 0.98'; 
%tSSSstr = ' -st 10 -corr 0.9'; 
tSSSstr = '';
dsstr = ''; % dsstr = ' -ds 4';  % downsample to 250Hz - can cause MF crash?
incEEG = 0;
TransTarget = '';  % Safer not to align sessions since different experiments (trans default anyway)
TransDefaultFlag = 0;
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

% if ~autobad_flag 
%     mad_bad = read_bad_madrid('MAD'); % hacky using Matlab's import (don't seem to resemble bad channels detected by maxfilter...?)
%     cam_bad = read_bad_madrid('CAM'); % hacky using Matlab's import
% end
% 
% bad_agreement = [];

% Below just done once
% eval(sprintf('!cp /imaging/rh01/Projects/AD_MCI/Collaborations/Biofind/CTB/ct_sparse_Madrid.fif %s/ct_sparse_CTB.fif',owd));
% eval(sprintf('!cp /imaging/rh01/Projects/AD_MCI/Collaborations/Biofind/CTB/Madrid_sss_cal_1.dat %s/sss_cal_CTB_1.dat',owd));  % First 13 controls had cal_2.dat...
% eval(sprintf('!cp /imaging/rh01/Projects/AD_MCI/Collaborations/Biofind/CTB/Madrid_sss_cal_2.dat %s/sss_cal_CTB_2.dat',owd));  % ...rest had cal_1.dat?
% 
% eval(sprintf('!cp /neuro/databases/ctc/ct_sparse.fif %s/ct_sparse_CBU.fif',owd));
% eval(sprintf('!cp /neuro/databases/sss/sss_cal.dat   %s/sss_cal_CBU.dat',owd));  
% 
% eval(sprintf('!cp /imaging/rh01/Projects/AD_MCI/DandF/Oxford1/ct_sparse.fif %s/ct_sparse_OHBA.fif',owd));
% eval(sprintf('!cp /imaging/rh01/Projects/AD_MCI/DandF/Oxford1/sss_cal.dat  %s/sss_cal_OHBA.dat',owd)); 

%cal_sub = 0;

parfor subnum = dosubs

    cal_ver  = '';
    
    if strcmp(site{subnum},'CTB')
%        cal_sub = cal_sub+1;
        cal_ver = sprintf('_%d',CTB.Biofind.Madrid.calfiles(subnum-1))
    end
    
    % using maxfilter to cut data doesn't work, so done outside with mne_python
%     eve_file = fullfile(bwd,subid{subnum},'ses-meg1','meg',sprintf('%s_ses-meg1_task-Rest_events.tsv',subid{subnum}))
%     [eve_srt,eve_end,~,~] = textread(eve_file,'%s%s%s%s');  %csvread(eve_file,1,0,[1 0 1 1])
%    skipstr = sprintf(' -nosss -skip 0 %d %d 9999',str2num(eve_srt{2}),str2num(eve_srt{2})+str2num(eve_end{2}));
%     finstr = [maxfstr skipstr filestr sprintf(' -v | tee %s.log',outfile)]
%     rik_eval(finstr);
%     filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);
%    outfile = fullfile(sub_dir{subnum},'resting'); 
    
    base_name = sprintf('%s_ses-meg1_task-Rest',subid{subnum});
    raw_file = fullfile(bwd,subid{subnum},'ses-meg1','meg',[base_name '_meg.fif'])
   
    basestr = sprintf(' -ctc %s/ct_sparse_%s.fif -cal %s/sss_cal_%s%s.dat',hwd,site{subnum},hwd,site{subnum},cal_ver);
    basestr = [basestr ' -linefreq 50 -hpisubt amp -force'];
    
%     movfile = fullfile(sub_dir{subnum},[base_name '_move.txt']); % This file will record translations between runs
%     rik_eval(sprintf('!touch %s',movfile));
%     rik_eval(sprintf('!echo ''%s'' >> %s',datestr(now),movfile));   
    
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
%         if ~exist(badfile,'file') | OverWrite
%             if strcmp(site{subnum},'CTB') % below needed because not all (CAM) subjects have manual bad definitions
%                 if strcmp(group{subnum},'patient'), grp = 'MCI'; else grp = 'CN'; end
%                 [~,id] = modrem2(subnum,42);  %% HACKY!!!
%                 sid = find(strcmp(sprintf('MAD-%s-%04d',grp,id),mad_bad.MADID));
%                 if ~isempty(sid), man_bad = mad_bad.BadChan{sid}; else man_bad = []; end
%             elseif strcmp(site{subnum},'CBU') | strcmp(site{subnum},'OHBA')
%                 if strcmp(group{subnum},'patient'), grp = 'AD'; else grp = 'CC'; end
%                 [~,id] = modrem2(subnum,42);  %% HACKY!!!
%                 sid = find(strcmp(sprintf('sub-%s%06d',grp,id),cam_bad.MADID));
%                 if ~isempty(sid), man_bad = cam_bad.BadChan{sid}; else man_bad = []; end
%             end
%         end
%         
%          tmp=dlmread(fullfile(sub_dir{subnum},'bad_autobad1.txt'),' '); % Nbuf = size(tmp,1)-1; Doesn't work if subset of channels bad for short time
%          tmp=reshape(tmp,1,prod(size(tmp)));
%          tmp=tmp(tmp>0); % Omit zeros (padded by dlmread):
%          [frq,allbad] = hist(tmp,unique(tmp));
%          badchans = allbad(frq>0.05*Nbuf);
%          
%          bad_agreement(subnum,:) = [length(badchans) length(man_bad) length(intersect(badchans,man_bad))];
% 
%         if isempty(man_bad)
%             warning('No manual bad channels provided, so using previous auto-detected ones')
%             man_bad = badchans;
%         end
%         
%         fp = fopen(badfile,'w');
%         for n=1:Nbuf  % just to match autodetect version
%             fprintf(fp,'%s\r\n',num2str(man_bad));
%         end
%         fclose(fp);
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
        %             transfstfile = [outfile '.fif'];
        %             if ~isempty(TransTarget) & ~strcmp(do_sessions{ses},TransTarget)
        %                 outfile = fullfile(sub_dir{subnum},sprintf('sss_trans_%s',TransTarget))
        %                 transtr = sprintf(' -trans %s/sss.fif',fullfile(sub_dir{subnum},TransTarget));
        %             else
        transtr = '';
        %                outfile = fullfile(sub_dir{subnum},sprintf('sss'))
        %            end
        
        skipstr = '';
%        compstr = sprintf(' -movecomp inter')
%         if strcmp(site{subnum},'CTB')
%             if ismember(ctb_dir,FailedMvComp)
%                 %                  skipstr = ' -skip 375 895' % doesn't help
%                 %                   compstr = sprintf(' -movecomp'); % doesn't help
%                 compstr = ''
%             end
%         else 
%             if ismember(subnum,CBUFailedMvComp)
                compstr = ''
%             end                
%         end
%         
        filestr = sprintf(' -f %s -o %s.fif',raw_file,outfile);
        finstr = [maxfstr filestr basestr badstr tSSSstr skipstr compstr origstr transtr dsstr sprintf(' -v | tee %s.log',outfile)]
        rik_eval(finstr);
    end
    
    %         if ~isempty(TransTarget) & ~strcmp(do_sessions{ses},TransTarget)
    %             %fprintf(fp,'\nTransfirst %s to %s: ',do_sessions{ses},TransTarget);
    %             rik_eval(sprintf('!echo ''Trans %s to %s'' >> %s',do_sessions{ses},TransTarget,movfile));
    %             rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',outfile,movfile));
    %         end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. trans to default helmet space (align across subjects)
    
    transdeffile = fullfile(sub_dir{subnum},sprintf('transdef_autobad%d',autobad_flag));
    if TransDefaultFlag
        if ~exist([transdeffile '.fif'],'file') | ~exist([transdeffile '.log'],'file') | OverWrite
            transtr = sprintf(' -trans default -origin %d %d %d -frame head -force',orig+[0 -13 6])
            filestr = sprintf(' -f %s.fif -o %s.fif',outfile,transdeffile);
            finstr = [maxfstr filestr transtr sprintf(' -v | tee %s.log',transdeffile)]
            rik_eval(finstr);
        end
        
        rik_eval(sprintf('!echo ''Transdef '' >> %s',movfile));
        %fprintf(fp,'\nTransDef %s: ',do_sessions{ses});
        rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',transdeffile,movfile));
    end
end

return



%% Collect up results
olddir = ''; % have to edit *.txt files in latest run to get code below to work, so use old versions (results should be same)
%olddir = 'MF0.98ST';
trdef = []; orig = []; bad = []; move = []; move2 = []; trl_info = []; time_info = [];
for subnum = 1:Nsub
    sub_dir = fullfile(mwd,olddir,sprintf('sub-Sub%02d',subnum));
    
    tmp = textread(fullfile(sub_dir,'trans_move.txt'),'%s');
    trdef(subnum,1) = str2num(tmp{4}); %tmp{end-1});
    
    %            tmp = load(fullfile(sub_dir,sprintf('%s_org.txt',subjects{o}{g}{sub})));
    tmp = load(fullfile(sub_dir,'fittmp.txt'));
    orig(subnum,:) = round(tmp(1:3)*1000);
    
    %    tmp = dlmread(fullfile(sub_dir,'bad.txt'),' '); % Nbuf = size(tmp,1)-1; Doesn't work if subset of channels bad for short time
    tmp = dlmread(fullfile(sub_dir,sprintf('bad_autobad%d.txt',autobad_flag)),' '); % Nbuf = size(tmp,1)-1; Doesn't work if subset of channels bad for short time
    tmp=reshape(tmp,1,prod(size(tmp))); tmp=tmp(tmp>0);
    [frq,allbad] = hist(tmp,unique(tmp));
    bad(subnum,1) = length(allbad(frq>0.05*Nbuf));
    
     S.logfname = fullfile(sub_dir,'sss.log');
     S.toplot = 0;
     [md,mr,sd,sr] = plot_maxfilter_move(S);
     move(subnum,:) = [md mr sd sr];

     if exist(fullfile(sub_dir,'headpos.txt'))
         tmp = textread(fullfile(sub_dir,'headpos.txt'),'%s');
         pos = [];
         for t=15:10:length(tmp)
             pos(end+1,:) = [str2num(tmp{t}) str2num(tmp{t+1}) str2num(tmp{t+2})];  % q4-q6 are T only? (page 77 of MF manual)
         end
         tmp=[];
         for t=1:(size(pos,1)-1)
             tmp(t) = sqrt(sum([pos(t+1,:)-pos(1,:)].^2))*1000;
         end
         move2(subnum,:) = [mean(tmp) std(tmp)];
     else
         move2(subnum,:) = [NaN NaN];
     end
    
    trl_file = fullfile(bwd,'BIDS',sprintf('sub-Sub%02d',subnum),'meg',sprintf('sub-Sub%02d_task-Rest_events_Madrid.tsv',subnum));
    if exist(trl_file)
        [ons,dur,~,tname] = textread(trl_file,'%s%s%s%s');
        ons = str2num(strvcat(ons{2:end}));
        dur = str2num(strvcat(dur{2:end}));
        if length(unique(dur))>1, error;
        elseif dur~=3.999; error;
        end
        trl_info(subnum,:) = [ons(1) size(ons,1) ons(end)];
    else
        trl_info(subnum,:) = [NaN NaN NaN];
    end
    
    % How extract time and date from a FIF file?
    raw_file = fullfile(bwd,'BIDS',sprintf('sub-Sub%02d',subnum),'meg',sprintf('sub-Sub%02d_task-Rest_meg.fif',subnum))
    eval(sprintf('!/imaging/local/software/utils/fiff_date_time %s > %s/rawdatetime.txt',raw_file,sub_dir))
    tmp = textread(sprintf('%s/rawdatetime.txt',sub_dir),'%s');    
    time_info(subnum,1) = str2num(tmp{end-4})+str2num(tmp{end-2})/60; % Time of Day
    time_info(subnum,2) = str2num(tmp{end-10}); % Day
    time_info(subnum,3) = str2num(tmp{end-8}); % Month
    time_info(subnum,4) = str2num(tmp{end-6})-2000; % Year
end

Gnames = {'ConCTB','PatCTB','ConCBU','PatCBU','PatOHBA'}; % check matches above!

X = zeros(Nsub,length(groups)+length(sites)); 
c=0; Gind={};
for s = 1:length(sites)
    for g = 1:length(groups)
        c=c+1;
        l = find(strcmp(group,groups{g}) & strcmp(site,sites{s}));
        Gind{end+1} = l;
        X(l,c)=1;
    end
end
sum(X)
X(:,5)=[]; % No OHBA patients
Oind = {};
for c=[1:4 6]
    Oind{end+1} = Gind{c};
end
Gind{4} = [Gind{4}; Gind{6}];  % !!!Add OHBA to CBU for some analyses

% X = zeros(Nsub,length(groups)+length(sites)); 
% for g = 1:length(groups)
%     l = find(strcmp(group,groups{g}));
%     X(l,g)=1;
% end
% for s = 1:length(sites)
%     l = find(strcmp(site,sites{s}));
%     X(l,s+length(groups))=1;
% end

%X = [X age/100 sexn ones(42*4,1)]; % Hacky assumption of how data organised!!! Could include age/sex
figure,imagesc(X);

valid_move = find(~isnan(move(:,1)));
invalid_move = find(isnan(move(:,1)));
% Problem with SVD is lose units - stick with translations?
% [u,s,v]=spm_svd(move(valid_move,:)); s(1)^2/sum(s.^2);  %s = full(diag(s)), 
% move_summary = NaN(Nsub,1);
% move_summary(valid_move) = -full(u(:,1));  % care minus!!!

cd(bwd)
y = [bad move(:,1:2) trdef time_info(:,1) trl_info]; Ylabel = {'Bad','MMove','SDMove','Position','TimeOfDay','DataOnset','NumEpochs','DataOffset'};
%y = cat(1,sf{:}); Ylabel = {'MagSF','GrdSF'};
for c=1:size(y,2)
    valid = find(~isnan(y(:,c)));
    yy = y(valid,c); XX = X(valid,:);
    Ylabel{c}
    [t,F,p,df,R2,cR2,B,r,aR2,iR2] = glm(yy,XX,[1 -1 1 -1 0]'); [F,2*p]  % Con vs Pat
    [t,F,p,df,R2,cR2,B,r,aR2,iR2] = glm(yy,XX,[1 1 -1 -1 0]'); [F,2*p]  % Mad vs Cam
    [t,F,p,df,R2,cR2,B,r,aR2,iR2] = glm(yy,XX,[1 -1 -1 1 0]'); [F,2*p]  % Interaction
    yyy=[];cind=[];my = [];sd=[]; sv={}; si={};
    for g = 1:size(XX,2)
        my(g) = mean(yy(find(XX(:,g))));
        sd(g) = std(yy(find(XX(:,g))));
        [sv{g},si{g}] = sort(yy(find(XX(:,g))));
        yyy = [yyy; yy(find(XX(:,g)))];
%        yyy{g} = yy(find(XX(:,g)));
        cind = [cind; ones(size(find(XX(:,g))))*g];
    end
    figure,boxplot(yyy,cind);
    my
    %sd
%    figure,distributionPlot(yyy)
    set(gca,'XTick',[1:5]);
    set(gca,'XTickLabel',Gnames);
    ylabel(gca,Ylabel{c},'FontSize',18)
    set(gca,'FontSize',18);
    axis([0.5 5.5 0 max(yy)])
%    eval(sprintf('print -dpng violin_%s.png',Ylabel{c}))
    eval(sprintf('print -dpng boxplot_%s.png',Ylabel{c}))
end

excsubs = [];
excsubs = [excsubs find(move(:,1)>=6 | move(:,2)>=6)]  % 166 ??
%excsubs = [excsubs find(bad>=20)]  % Too many if manual Madrid??
excsubs = [excsubs find(trdef>=30)]  % 159 only if do trans def analysis?


for c=1:size(move,2)
    move(isnan(move(:,c)),c) = nanmean(move(:,c));
end

mfconfs = [bad move(:,1:2) trdef];

cd(mwd)
save('MEG_confs','mfconfs')



