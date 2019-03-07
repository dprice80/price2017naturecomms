clear all
close all
clc

load('/imaging/dp01/results/av_nature/GrdsBothAV.mat')
addpath /home/dp01/matlab/lib
addpath /imaging/camcan/sandbox/projects/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/scripts/av_nature/
addpath /imaging/dp01/toolboxes/FastICA_25/
addpath /home/dp01/matlab/lib/fitgaussian/
addpath /home/dp01/matlab/export_fig/
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))
addpath /imaging/camcan/QueryFunction/QueryFun_v1/
loadspm 12

[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');

Info = LoadSubIDs;
% Select Audio 1 or Visual 2

fitwindows{1} = [0 140];
fitwindows{2} = [0 140];

lab = {'Aud' 'Vis'};
Ncomps = 1;
for av = 2
    % ERP Fits
    clear Fit FitTimeSeries
    close all
    BadData = zeros(size(GrdsBothAV,3),1);
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,:,av),3));
        for ii = 1:size(GrdsBothAV,3)
            x = squeeze(GrdsBothAV(pci,:,ii,av)); % use original data
            showplots = 0; % 3 = show only final fit
            F = ERPfit(x', y, linspace(-100, 500, length(y)), 50, fitwindows{av}, showplots, 1e-6);
            Fit(pci,ii,:) = [F.shift F.stretch F.amp F.error sqrt(F.R2) F.Nsteps];
            
            title(['Subject ', Info.SubCCIDc{ii}, ' index = ', num2str(ii), ' PCI = ', num2str(pci), 'R = ', num2str(Fit(pci,ii,5)), 'Amp = ', num2str(Fit(pci,ii,3),'%2.4f')])
            if showplots == 3
                export_fig(['/imaging/dp01/results/av_nature/ERPFittingChecks' lab{av} '-subid-' num2str(ii,'%2.3d') '.bmp'],'-bmp')
            end
            disp(str(pci, ' ', ii,' Fit = Delay: ',Fit(pci,ii,1),' Stretch: ',Fit(pci,ii,2), ' R2: ',Fit(pci,ii,5)^2,' Nsteps: ',Fit(pci,ii,6)))
        end
        
        % Save Data
        if av == 1
            SaveValsAs = '/imaging/dp01/results/av_nature/Fits_AVpassiveAudio_shotwin.mat';
        else
            SaveValsAs = '/imaging/dp01/results/av_nature/Fits_AVpassiveVisual_shortwin.mat';
        end
        
        FitOut = zeros(Ncomps,708,6);
        FitOut(pci,FindFiles,:) = squeeze(Fit);
        FitOutMat = FitOut;
        FitOutCell{pci} = mat2dataset(squeeze(FitOut(1,:,:)),'VarNames',{'Latency' 'Stretch' 'Scaling' 'RMSE' 'R' 'NstepsFit'},'ObsNames',Info.SubCCIDc);
    end
    
    FitOut = FitOutCell{1}; % if there are any scripts still using FitOut, then they will get the correct PC (first, rather than the last one).
    save(SaveValsAs,'FitOutCell','FitOut','FindFiles','y','FitOutMat')
    
end















