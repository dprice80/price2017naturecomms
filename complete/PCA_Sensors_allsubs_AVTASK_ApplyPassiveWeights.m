clear all
close all
clc

addpath /home/dp01/matlab/lib
addpath /home/dp01/matlab/lib/paruly/
addpath /imaging/camcan/sandbox/projects/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/scripts/proj1/
addpath /imaging/dp01/toolboxes/FastICA_25/
addpath /home/dp01/matlab/lib/fitgaussian/
addpath /home/dp01/matlab/export_fig/
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))

Info = LoadSubIDs;

loadspm 12

[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/<CCID>/mefMspm12_transdef_transrest_mf2pt2_task_raw.dat');

Dalls = zeros(306,601,4,length(FindFiles));
Files = FileList(FindFiles);

for ii = 1:length(FindFiles)
    disp(ii)
    D = spm_eeg_load(FileList{FindFiles(ii)});
    Di = mean(D(1:306,:,1:3), 3);
    Dalls(:,:,1,ii) = Di-repmat(mean(Di(:, 1:100), 2), [1 size(Di ,2)]);
    Di = D(1:306,:,4);
    Dalls(:,:,2,ii) = Di-repmat(mean(Di(:, 1:100), 2), [1 size(Di ,2)]);
end
clear Di

% Load Data

[FileList] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/<CCID>/efMspm12_transdef_transrest_mf2pt2_task_raw.mat');

%
D = spm_eeg_load(FileList{1});
Cindex = D.coor2D(D.indchantype('MEGPLANAR'));
Clabels = D.chanlabels(D.indchantype('MEGPLANAR'));
Cindex = Cindex(:,1:2:end);
Clabels = Clabels(1:2:end);
Age = Info.Age(FileCheck);

Dalls = Dalls(103:306,:,:,:);


%%
Nsubs = size(Dalls, 4);
xwin = [-50 400];
lx = [0 0;100 100;200 200; 300 300];
ly = [18 88;18 88; 18 88; 18 88];
    
for av = 1:2

    Dallrs = reshape(squeeze(Dalls(:,:,av,:)),size(Dalls,1),size(Dalls,2)*size(Dalls,4));

    % Load PCA components from passive session
    if av == 1
        load /imaging/dp01/camcan/meg/avpassive/PCAcomps_audio.mat % load weights from passive
    else
        load /imaging/dp01/camcan/meg/avpassive/PCAcomps_visual.mat
    end
    GrdsBoth = Dallrs'*ug;
    % Variance explained by each component
    lgsum = full(sum(diag(lg)));
    varexp = full(diag(lg)./lgsum);
    
    GrdsBoth = full(GrdsBoth)';
    GrdsBoth = reshape(GrdsBoth,10,601,Nsubs);
    ug = reshape(full(ug),2,102,10);
    ug = sqrt(sum(ug.^2,1));
    ug = reshape(ug,102,10);
    disp('Done')
    
    Nbins = size(Dalls, 4);
    
    % Average age bins into new matrix.
    % Otherwise just copy data
    % Correct for audio/visual stimulus delay (13 aud, 34 vis)+23ms already accounted for.
    cs = [-13 -34]+23;
    GrdsBoth = circshift(GrdsBoth, [0 cs(av) 0]);
    
    GrdsBothSG = GrdsBoth;
    AgeBins = Age;
    
    % Sorted With AgeBins
    figure(200)
    set(gcf,'position',[26 419 1380 660],'color','w');
    SortingVar = AgeBins;
    GaussConv = 5;

    for ii = 1:5
        pc = shiftdim(ug(:,ii)');
        
        [IM,f] = spm_eeg_plotScalpData(pc,Cindex,Clabels);
        close
        
        subplot(3,5,ii)
        IMmask = repmat(~isnan(IM),[1 1 3]);
        clims = [-quantile(IM(:), 0.999) quantile(IM(:), 0.999)];
        IMRGB = dp_Convert2DToColor(IM,clims,paruly);
        IMRGB(IMmask == 0) = 1;
        imagesc(flipdim(IMRGB,1));
        colormap(paruly)
        clims2 = [-quantile(IM(:), 0.999) quantile(IM(:), 0.999)];
        caxis(clims)
        axis off
        
        subplot(3,5,ii+5)
        pl = plot(linspace(-100,500,601),squeeze(mean(GrdsBothSG(ii,:,:),3)));
        xlim(xwin);
        set(pl,'color',[0 160 255]/255,'linewidth',2);
        xlabel('Time (ms)');ylabel('Amplitude (A.U)')
        set(gca, 'XTick',[0, 100, 200, 300, 400])
        ct = get(gca,'YTick');
        set(gca, 'YTick',[min(ct), 0, max(ct)])
        set(gca,'YtickLabel',{'-', '0', '+'})
        
        subplot(3,5,ii+10)
        Np = size(GrdsBothSG,2);
        t = linspace(-100,500,601)';
        imdat = squeeze(GrdsBothSG(ii,:,:))';
        [Sorted, SortInd] = sort(SortingVar);
        if GaussConv > 0
            h = fspecial('gaussian', size(imdat), GaussConv); % create a gaussian blurring kernel
            imdat = conv2(imdat(SortInd,:),h,'same'); % convolve that kernel with the gradient image.
        else
            imdat = imdat(SortInd,:); % convolve that kernel with the gradient image.
        end
        i1 = imagesc(t,Sorted,imdat);
        colormap(paruly)
        set(gca,'ydir','normal','xlim',xwin)
        caxis(quantile(imdat(:),[0.001 0.999]))
        set(gca, 'XTick',[0, 100, 200, 300, 400])
        ln = line(lx',ly');
        set(ln,'Color',[0 0 0],'LineStyle','--')
        xlabel('Time (ms)');ylabel('Age (years)')
    end
    %
        if av == 1
            export_fig('/imaging/dp01/results/av_nature/Figures/PCA_ERP_avtask_aud_Paper_PassiveWeights','-pdf')
        else
            export_fig('/imaging/dp01/results/av_nature/Figures/PCA_ERP_avtask_vis_Paper_PassiveWeights','-pdf')
        end
    GrdsBothAV(:,:,:,av) = GrdsBothSG;
end

save('/imaging/dp01/results/av_nature/GrdsBothAVtask_ApplyPassiveWeights.mat','GrdsBothAV')


%% Load Presaved data (if data is already saved, start from here)

clear all
close all
clc

load('/imaging/dp01/results/av_nature/GrdsBothAVtask_ApplyPassiveWeights.mat')
addpath /home/dp01/matlab/lib
addpath /imaging/camcan/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/scripts/av_nature/
addpath /imaging/dp01/toolboxes/FastICA_25/
addpath /home/dp01/matlab/lib/fitgaussian/
addpath /home/dp01/matlab/export_fig/
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))

[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avtask/<CCID>/mefMspm12_transdef_transrest_mf2pt2_task_raw.mat');

Info = LoadSubIDs;
% Select Audio 1 or Visual 2

Ncomps = 1;
for av = 1:2
    % ERP Fits
    clear Fit FitTimeSeries
    close all
    FitOut = zeros(Ncomps,708,6);
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,:,av),3));
        for ii = 1:size(GrdsBothAV,3)
            x = squeeze(GrdsBothAV(pci,:,ii,av)); % use original data
            showplots = 0; % 3 = show only final fit
            F = ERPfit(x', y, linspace(-100, 500, length(y)), 50, [0 400], 0, 1e-4);
            Fit(pci,ii,:) = [F.shift F.stretch F.amp F.error sqrt(F.R2) F.Nsteps];
            disp(str(pci, ' ', ii,' Fit = Delay: ',Fit(pci,ii,1),' Stretch: ',Fit(pci,ii,2), ' R2: ',Fit(pci,ii,5)^2,' Nsteps: ',Fit(pci,ii,6)))
        end
        
        % Save Data
        FitOut(pci,FindFiles,:) = squeeze(Fit(pci,:,:));
        FitOutMat = FitOut;
        FitOutCell{pci} = mat2dataset(squeeze(FitOut(1,:,:)),'VarNames',{'Latency' 'Stretch' 'Scaling' 'RMSE' 'R' 'NstepsFit'},'ObsNames',Info.SubCCIDc);
    end
    
    if av == 1
        SaveValsAs = '/imaging/dp01/results/av_nature/Fits_AVtaskAud_ApplyPassive.mat';
    else
        SaveValsAs = '/imaging/dp01/results/av_nature/Fits_AVtaskVis_ApplyPassive.mat';
    end
    
    FitOut = FitOutCell{1}; % if there are any scripts will using FitOut, then they will get the correct PC (rather than the last one).
    save(SaveValsAs,'FitOutCell','FitOut','FindFiles','y','FitOutMat')
end

