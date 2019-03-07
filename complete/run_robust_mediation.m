%% Load Data

% Version 2 uses the updated fitting procedure, ERP_Fitting_v3, which does
% not use absolute value of R for fitting. Outliers are also removed in
% both lat and str variables (if an outlier is detected in lat it is also
% removed in str) because the data point is considered bad. 

clear all
close all
clc

addpath /home/dp01/matlab/NIFTI_20110610/
addpath /imaging/camcan/QueryFunction/QueryFun_v1/
addpath /home/dp01/matlab/lib/
addpath /imaging/dp01/scripts/av_nature/
addpath(genpath('/home/dp01/matlab/MediationToolbox/mediation_toolbox/'))
addpath(genpath('/home/dp01/matlab/RobustToolbox/robust_toolbox/'))
addpath(genpath('/home/dp01/matlab/SCN_Core_Support'))
loadspm 12

%% Make cortical white matter mask.

NII = load_nii('/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii');
NII.img(29:65, 29:59, 7:32) = 0;
save_nii(NII, '/imaging/dp01/templates/JHU-ICBM-labels-2mm_cortical_mask.nii')
view_nii(NII)

%%
[GM_fl, GM_index] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002/avpassive/<CCID>/anatomicals/GM/s122m*.nii');
[MK_fl, MK_index] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002/avpassive/<CCID>/anatomicals/MK/s12swdki_MK.nii');

Info = LoadSubIDs;

[TIV_fl, TIV_index] = CCQuery_QuickCheck('/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_structuralstats_00001/<MRIID>/structurals/*.mat');
TIV = zeros(708,1);
disp('Loading TIV')
for ii = 1:length(TIV_fl)
    if TIV_index(ii)
        load(TIV_fl{ii});
        TIV(ii) = S.TIV.mm3;
    end
end
TIV(TIV == 0) = NaN;

ERPvis  = load('/imaging/dp01/results/av_nature/Fits_AVpassiveVisual.mat');
ERPaud  = load('/imaging/dp01/results/av_nature/Fits_AVpassiveAudio.mat');

% Remove zeros (these are just subjects that had no data)
ind = ERPvis.FitOut.NstepsFit == 0;
ERPvis.FitOut(ind, :) = dataset(NaN);

% Remove zeros
ind = ERPaud.FitOut.NstepsFit == 0;
ERPaud.FitOut(ind, :) = dataset(NaN);

% Remove cases where amplitude was less than 0
ind = ERPaud.FitOut.Scaling < 0;
ERPaud.FitOut(ind, :) = dataset(NaN);

ind = ERPvis.FitOut.Scaling < 0;
ERPvis.FitOut(ind, :) = dataset(NaN);

VisLat = dp_removeoutliers(ERPvis.FitOut.Latency);
VisStr = dp_removeoutliers(ERPvis.FitOut.Stretch);
AudLat = dp_removeoutliers(ERPaud.FitOut.Latency);
AudStr = dp_removeoutliers(ERPaud.FitOut.Stretch);
AudAmp = dp_removeoutliers(ERPaud.FitOut.Scaling);
VisAmp = dp_removeoutliers(ERPvis.FitOut.Scaling);

[AudLat, AudStr, ind] = complete(AudLat, AudStr);
AudLat(ind == 0) = NaN;
AudStr(ind == 0) = NaN;

[VisLat, VisStr, ind] = complete(VisLat, VisStr);
VisLat(ind == 0) = NaN;
VisStr(ind == 0) = NaN;

noOLind = ~isnan(VisLat) & ~isnan(VisStr) & ~isnan(AudStr) & ~isnan(AudLat);

loadspm 12
spm_jobman('initcfg')

hand = Info.AddData.handedness;
gender = Info.GenderNum;


%% Run Mediation Analysis on Whole Brain - MK

% Load NII data
[MK_fl, MK_index] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002/avpassive/<CCID>/anatomicals/MK/s12sw*.nii');

include = MK_index;
ind(540) = 0; % bad subject data
ind(Info.AddData.DTIStripeExclude) = 0;
NIImask = load_nii('/imaging/dp01/templates/JHU-ICBM-LGN-V1_tractROIRS.nii');
MK_ff = find(include);

f1 = find(NIImask.img ~= 0,1,'first');
WMmean = nan(708,1);
WMstd = nan(708,1);
WMmin = nan(708,1);
for ii = 1:length(MK_ff)
    disp(str('Loading ', ii))
    NII = load_nii(MK_fl{MK_ff(ii)});
    WMmean(MK_ff(ii)) = mean(NII.img(NIImask.img ~= 0));
    WMstd(MK_ff(ii)) = std(NII.img(NIImask.img ~= 0));
    WMmin(MK_ff(ii)) = min(NII.img(NIImask.img ~= 0));
    WMmedian(MK_ff(ii)) = median(NII.img(NIImask.img ~= 0));
    WMvox(MK_ff(ii)) = NII.img(f1);
end
wmclean = dp_removeoutliers(WMmean);
include(isnan(wmclean)) = 0;


NIImask = load_nii('/imaging/dp01/templates/JHU-ICBM-labels-2mm.nii');
NIImask.img(NIImask.img > 0) = 1;
save_nii(NIImask,'/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii')
scn_map_image('/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii','/imaging/dp01/cc700meg-mf22-spm12/release002/avpassive/CC110033/anatomicals/MK/sswdki_MK.nii','write','/imaging/dp01/templates/JHU-ICBM-labels-2mm_maskRS.nii');
scn_map_image('/imaging/dp01/templates/JHU-ICBM-LGN-V1_tractROI.nii','/imaging/dp01/cc700meg-mf22-spm12/release002/avpassive/CC110033/anatomicals/MK/sswdki_MK.nii','write','/imaging/dp01/templates/JHU-ICBM-LGN-V1_tractROIRS.nii');

exclude = double(include);
exclude(exclude == 0) = NaN;

[Y,C,X, ind] = complete(VisLat, TIV, Info.Age, (1:708)', exclude);

mkdir /imaging/dp01/results/av_nature/vbm/mediation/
cd /imaging/dp01/results/av_nature/vbm/mediation/
M = MK_fl(ind);
M = char(M);
names = {'X' 'Y' 'M'};

mkdir /imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS
cd /imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS
medbrain = mediation_brain(X,Y,M,'covs',C,'names',names,'mask','/imaging/dp01/templates/JHU-ICBM-labels-2mm_maskRS.nii');
thresh = mediation_brain_corrected_threshold('fdr', 'mask', '/imaging/dp01/templates/JHU-ICBM-labels-2mm_maskRS.nii');
[clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh', thresh.fdr_p_thresh, 'size', [5 5 5],'overlay','/imaging/dp01/templates/MNI152_T1_2mm.nii');


%%

running = 1;

while done
    try
        dp_matlabpool_start(75)
        
        dp_mediation_parfor(MK_fl, MK_index, VisLat, TIV, Info.Age,...
            '/imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/VisLat/',...
            '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii');
        
        dp_mediation_parfor(GM_fl, GM_index, AudStr, TIV, Info.Age,...
            '/imaging/dp01/results/av_nature/vbm/mediation/GM_RobustParfor/AudStr/',...
            '/imaging/camcan/templates/HarvardOxford-combo-maxprob-thr25-2mm_nobrainstem.nii');
        
        dp_mediation_parfor(MK_fl, MK_index, AudStr, TIV, Info.Age,...
            '/imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/AudStr/',...
            '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii');
        
        dp_mediation_parfor(GM_fl, GM_index, VisLat, TIV, Info.Age,...
            '/imaging/dp01/results/av_nature/vbm/mediation/GM_RobustParfor/VisLat/',...
            '/imaging/camcan/templates/HarvardOxford-combo-maxprob-thr25-2mm_nobrainstem.nii');
        
        running = 0;
    catch ME
        fid = fopen('/home/dp01/parfor_errorlog.txt','a');
        fprintf(fid, 'Error occured at %s\t%s\n', datestr(now), ME.message);
        fclose(fid);
    end
end

matlabpool close


%% Recreate pathdata_sum
dp_mediation_create_pathdata_sum('/imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/AudStr/')
dp_mediation_create_pathdata_sum('/imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/VisLat/')




%% Compile results and display significant clusters MK

applymasks = {
    '/imaging/dp01/templates/JHU-ICBM-labels-2mm_cortical_mask.nii'
    '/imaging/dp01/templates/JHU-ICBM-labels-2mm_cortical_mask.nii'
    };

results_folders = {
    '/imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/VisLat/'
    '/imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/AudStr/'
    };

for ri = 1
    applymask = applymasks{ri};
    results_folder = results_folders{ri};
    
    cd(results_folder)
    !rm *.hdr *.img
    
    dp_mediation_compile_results(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii',...
        applymask);
    
    thresh = mediation_brain_corrected_threshold('fdr', 'mask', applymask);
    mediation_brain_corrected_threshold('fdr', 'mask', applymask)
    [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh',...
        thresh.fdr_p_thresh/4, 'size', 250, 'overlay','/imaging/dp01/templates/MNI152_T1_2mm.nii', 'prune');
    
    % Save Mask
    dp_mediation_create_clustermask(applymask, clpos_data, clneg_data)
    
    dp_mediation_compile_results_masked(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii',...
        [pwd '/thresholdmask.hdr'],'effectsize');
    
    dp_mediation_compile_results_masked(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii',...
        [pwd '/thresholdmask.hdr'],'paths');
    
    dp_mediation_compile_results_unmasked(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii' ,'effectsize');
end



%% GM
applymask = '/imaging/dp01/templates/HOA_nobrainstem_mask.nii';
% applymask = '/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii';

results_folders = {
    '/imaging/dp01/results/av_nature/vbm/mediation/GM_RobustParfor/AudStr/'
    '/imaging/dp01/results/av_nature/vbm/mediation/GM_RobustParfor/VisLat/'
    };

for ri = 1
    results_folder = results_folders{ri};

    cd(results_folder)
    !rm *.hdr *.img
    
    dp_mediation_compile_results(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        applymask,...
        applymask);
    
    thresh = mediation_brain_corrected_threshold('fdr', 'mask', applymask);
    [clpos, clneg, clpos_data, clneg_data] = mediation_brain_results('all', 'thresh',...
        thresh.fdr_p_thresh/4, 'size', 250, 'overlay','/imaging/dp01/templates/MNI152_T1_2mm.nii', 'prune');
    
    % Save Mask
    dp_mediation_create_clustermask(applymask, clpos_data, clneg_data)
    
    dp_mediation_compile_results_masked(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        applymask,...
        [pwd '/thresholdmask.hdr'],'effectsize');
    dp_mediation_compile_results_masked(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        applymask,...
        [pwd '/thresholdmask.hdr'],'effectsize');
    dp_mediation_compile_results_unmasked(results_folder,...
        '/imaging/dp01/results/av_nature/vbm/mediation/JHUall_OLS/',...
        applymask,'effectsize');
end


%% Time Left GM
% cd /imaging/dp01/results/av_nature/vbm/mediation/GM_RobustParfor/AudStr/
% 
% % Resume Analysis
% for li = 1:2
%     clear pathdata
%     NIImask = load_nii('/imaging/camcan/templates/HarvardOxford-combo-maxprob-thr25-2mm_nobrainstem.nii');
%     maskind = find(NIImask.img);
%     N = length(maskind);
%     pathdata_sum = zeros(N,5,2);
%     
%     D = dir('pathdata_w*.mat');
%     for ii = 1:length(D)
%         load(D(ii).name)
%         pathdata_sum = pathdata_sum + pathdata;
%     end
%     completen = pathdata_sum(:,1,1) ~= 0;
%     pcomp(li) = sum(completen)./N*100;
%     fprintf('%2.2f%% complete', pcomp(li))
%     pause(60)
% end
% hoursleft = (100-pcomp(2)) ./ ((pcomp(2)-pcomp(1))) / 60;
% 
% % % End Resume Analysis
% 
% 
% %% Time left in MK
% % cd /imaging/dp01/results/av_nature/vbm/mediation/MK_RobustParfor/AudStr/
% cd /imaging/dp01/results/av_nature/vbm/mediation/standard_smoothing/MK_Robust/AudStr/
% % Resume Analysis
% for li = 1:2
%     clear pathdata
%     NIImask = load_nii('/imaging/dp01/templates/JHU-ICBM-labels-2mm_mask.nii');
%     maskind = find(NIImask.img);
%     N = length(maskind);
%     pathdata_sum = zeros(N,5,2);
%     
%     D = dir('pathdata_w*.mat');
%     for ii = 1:length(D)
%         load(D(ii).name)
%         pathdata_sum = pathdata_sum + pathdata;
%     end
%     completen = pathdata_sum(:,1,1) ~= 0;
%     pcomp(li) = sum(completen)./N*100;
%     fprintf('%2.2f%% complete', pcomp(li))
%     pause(60)
% end
% hoursleft = (100-pcomp(2)) ./ ((pcomp(2)-pcomp(1))) / 60;
% % % End Resume Analysis
% 
