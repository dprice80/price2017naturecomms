clear all
close all
clc

addpath /imaging/camcan/QueryFunction/QueryFun_v1/

Info = LoadSubIDs;

DAT = [];
DAT.SessionList = {
    'MEGmat' '/imaging/camcan/cc700-meg/pipeline/release002/data_trans/aamod_meg_denoise_ICA_2_applytrajectory_00002/*<ccID>/passive/Mspm12_transdef_transrest_mf2pt2_passive_raw.mat'
    };
DAT = CCQuery_CheckFiles(DAT);

ff = DAT.AllExistSubInd;

for ii = 1:length(ff)
    mkdir(['/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/' Info.SubCCIDc{ff(ii)}])
    unix(['ln -s ' DAT.FileNames.MEGmat{ff(ii)} ' /imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/' Info.SubCCIDc{ff(ii)} '/Mspm12_transdef_transrest_mf2pt2_passive_raw.mat'])
end



%% Filter and Epoch Data 

addpath /home/dp01/matlab/lib/
addpath /imaging/dp01/scripts/av_nature/batch_functions/
loadspm 12
spm('defaults', 'EEG');
Info = LoadSubIDs;

dp_matlabpool_start(96)

% Filter
DAT = [];
DAT.SessionList = {
    'MEGmat' '/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/Mspm12_transdef_transrest_mf2pt2_passive_raw.mat'
    };
DAT = CCQuery_CheckFiles(DAT);
ff = DAT.AllExistSubInd;
FileNames = DAT.FileNames.MEGmat(ff);

parfor ii = 1:length(ff)
    filter_job(FileNames(ii));
end


%% Find sound card change

DATsc = [];
DATsc.ListFound = 1;
DATsc.SelectFirstFile = 1;
DATsc.SessionList = {
    'meg_passive_raw' '/megdata/camcan/camcan/<MEGID>_<ccID>/*/passive_raw.fif'
    };
DATsc = CCQuery_CheckFiles(DATsc);
% Convert string to a number (in number of days from the birth of baby
% Jesus) and which files have created dates before that time.
BeforeSoundCard = datenum('12/08/2011') >= DATsc.FileInfo.meg_passive_raw.datenum;
AfterSoundCard = datenum('12/08/2011') < DATsc.FileInfo.meg_passive_raw.datenum;
% Create index before and after that matches the filecheck


%% Epoch
addpath /home/dp01/matlab/lib/
addpath /imaging/dp01/scripts/av_nature/batch_functions/
loadspm 12
spm('defaults', 'EEG');
Info = LoadSubIDs;

dp_matlabpool_start(96)
DAT = [];
DAT.SessionList = {
    'MEGmat' '/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/fMspm12_transdef_transrest_mf2pt2_passive_raw.mat'
    };
DAT = CCQuery_CheckFiles(DAT);
ff = DAT.AllExistSubInd;

FileNames = DAT.FileNames.MEGmat(ff);

CCID = Info.SubCCIDc(ff);

MES = {};

soundcard = AfterSoundCard(ff);

parfor ii = 1:length(ff)
    if soundcard(ii) == 1
        trlshift = 13;
    else
        trlshift = 0;
    end
    
    try
        epoch_job(FileNames(ii), trlshift);
    catch ME
        MES{ii} = ME;
    end
end

matlabpool close

%%

[fn,fc,ff] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');


