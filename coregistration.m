%%  C O R E G I S T R A T I O N 

%% Copy and unzip Structral data

copy_MRIs = 1; % set it to 0, if you ran this bit before

if copy_MRIs
    
    parfor sub=1:length(subs)
        bidsT1 = fullfile(bidspth,subdir{sub},'ses-meg1','anat',[subdir{sub} '_ses-meg1_T1w.nii.gz']);
        T1file = fullfile('/imaging/dv01/Processed',subdir{sub},[subdir{sub} '_ses-meg1_T1w.nii'])
        if exist(bidsT1,'file') & ~exist(T1file,'file')
            copyfile(bidsT1,[T1file '.gz']);
            gunzip([T1file '.gz'])
            delete([T1file '.gz'])
        end
        
        
    end
    
end
%% Run COREGISTRATION using SPM- Fiducial only


parfor sub=1:length(subs)
    
    infile = sprintf('/imaging/dv01/Processed/%s/ffdspmeeg',subdir{sub});
    D = spm_eeg_load(infile);
    
    T1file = fullfile('/imaging/dv01/Processed',subdir{sub},[subdir{sub} '_ses-meg1_T1w.nii']);
    if exist(T1file,'file') % If a subject does not have sMRI
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
    D.save;
    
end



