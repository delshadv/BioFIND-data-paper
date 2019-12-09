%% Beamforming using M , G and fusion of them

p          = parcellation('/imaging/dv01/Toolboxes/osl/parcellations/fMRI_parcellation_ds8mm.nii.gz');
mni_coords = p.template_coordinates;

modalities = {{'MEG','MEGPLANAR'}, {'MEG'}, {'MEGPLANAR'}}; % Might be possible to store multiple montages per modality

parfor sub = 1:length(subs)
    
    infile = sprintf('/imaging/dv01/Processed/%s/ffdspmeeg',subdir{sub}); %Choose to use epoched or continuous data
    
    D = spm_eeg_load(infile);
    D = osl_filter(D,[2 45]); % If you need filtration before source-localisation
    
    SS = struct;
    SS.timespan          = [0 Inf];
    SS.pca_order         =  100;
    SS.type              = 'Scalar';
    SS.inverse_method    = 'beamform';
    SS.prefix            = '';
    
    for m = 1:length(modalities)
        
        S = []; S.D = D;
        
        if length(modalities{m})>1
            SS.fuse = 'all';
            SS.modalities = modalities{m};
            SS.dirname = [infile '_' strcat(modalities{m}{:}) '_BFn'];
            S.outfile = [infile '_' strcat(modalities{m}{:})];
        else
            SS.fuse = 'no';
            SS.modalities = modalities{m}{1};
            SS.dirname = [infile '_' modalities{m}{1} '_BFn'];
            S.outfile = [infile '_' modalities{m}{1}];
        end
        
        DS = spm_eeg_copy(S);
        
        DS = osl_inverse_model(DS, mni_coords, SS);
        
        % Select montage %
        DS = DS.montage('switch',1); % Select 0/1 for using unnormalised/normalised weights
        DS.save;
        
        %Ectract Region of intrest (ROIs)
        DS = ROInets.get_node_tcs(DS,p1.to_matrix(p1.binarize),'pca');
        
        DS = DS.montage('switch',3);
        DS.save;
        
    end
    
end