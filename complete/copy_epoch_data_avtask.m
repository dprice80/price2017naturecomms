clear all
close all
clc

addpath /imaging/camcan/QueryFunction/QueryFun_v1/
Info = LoadSubIDs;

DAT = [];
DAT.ListFound = 1;
DAT.SessionList = {
    %     'MEGmat'     '/imaging/camcan/cc700-meg/pipeline/release002/data_trans/aamod_meg_denoise_ICA_2_applytrajectory_00002/*<ccID>/passive/Mspm12_transdef_transrest_mf2pt2_passive_raw.mat'
    'MEGtask'    '/imaging/camcan/cc700-meg/pipeline/release002/data_trans/aamod_meg_denoise_ICA_2_applytrajectory_00002/<MEGID>_<ccID>/task/Mspm12_transdef_transrest_mf2pt2_task_raw.mat'
    };
DAT = CCQuery_CheckFiles(DAT);

ff = DAT.AllExistSubInd;

for ii = 1:length(ff)
    %     mkdir(['/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/' Info.SubCCIDc{ff(ii)}])
    %     unix(['ln -s ' DAT.FileNames.MEGmat{ff(ii)} ' /imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/' Info.SubCCIDc{ff(ii)} '/Mspm12_transdef_transrest_mf2pt2_passive_raw.mat'])
    mkdir(['/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/' Info.SubCCIDc{ff(ii)}])
    unix(['ln -s ' DAT.FileNames.MEGtask{ff(ii)} ' /imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/' Info.SubCCIDc{ff(ii)} '/Mspm12_transdef_transrest_mf2pt2_task_raw.mat'])
end



%
addpath /home/dp01/matlab/lib/
addpath /imaging/dp01/scripts/av_nature/batch_functions/
loadspm 12
spm('defaults', 'EEG');
Info = LoadSubIDs;

dp_matlabpool_start(96)

% Filter
DAT = [];
DAT.SessionList = {
    'MEGtask' '/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/<CCID>/Mspm12_transdef_transrest_mf2pt2_task_raw.mat'
    };
DAT = CCQuery_CheckFiles(DAT);
ff = DAT.AllExistSubInd;
FileNames = DAT.FileNames.MEGtask(ff);

parfor ii = 1:length(ff)
    filter_job(FileNames(ii));
end

AfterSoundCard = zeros(1,708);

%% Epoch

DAT = [];
DAT.SessionList = {
    'MEGtask' '/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/<CCID>/fMspm12_transdef_transrest_mf2pt2_task_raw.mat'
    };
DAT = CCQuery_CheckFiles(DAT);
ff = DAT.AllExistSubInd;

FileNames = DAT.FileNames.MEGtask(ff);

CCID = Info.SubCCIDc(ff);
dp_matlabpool_start(3)
MES = {};

soundcard = AfterSoundCard(ff);

parfor ii = 1:length(ff)
    
    if soundcard(ii) == 1
        trlshift = 13;
    else
        trlshift = 0;
    end
    
    try
        epoch_job_avtask(FileNames(ii));
    catch ME
        MES{ii} = ME;
    end
end


%%

[fn,fc,ff] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');


