clear all
close all
clc

addpath /home/dp01/matlab/lib
addpath /imaging/camcan/sandbox/projects/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/scripts/av_nature/
addpath /imaging/dp01/scripts/av_nature/batch_functions/
addpath /imaging/dp01/toolboxes/FastICA_25/
addpath /home/dp01/matlab/lib/fitgaussian/
addpath /home/dp01/matlab/export_fig/
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))

loadspm 12

%% Load principal component

[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.dat');

dp_matlabpool_start(96)

Dalls = zeros(306,601,2,length(FindFiles));
Files = FileList(FindFiles);

parfor ii = 1:length(FindFiles)
    disp(['Loading: ' Files{ii}])
    CreateAverages(Files{ii})
end

ind = 0;
for ii = 1:length(FindFiles)
    ind = ind + 1;
    disp(['Loading into Dalls: ' fullfile(fileparts(Files{ii}),'trialaverages.mat')])
    load(fullfile(fileparts(Files{ii}),'trialaverages.mat'))
    Dalls(:,:,1,ind) = DmeanA(1:306,:);
    Dalls(:,:,2,ind) = DmeanV(1:306,:);
end

save('/imaging/dp01/results/av_nature/Dalls_avpassive_allchans.mat','Dalls','FileCheck','FileList','FindFiles')

% matlabpool close


%% Perform PCA and Save
clear all
close all
clc

addpath /home/dp01/matlab/lib
addpath /imaging/camcan/sandbox/projects/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/scripts/av_nature
addpath /imaging/dp01/toolboxes/FastICA_25/
addpath /home/dp01/matlab/lib/fitgaussian/
addpath /home/dp01/matlab/export_fig/
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))
addpath /home/dp01/matlab/lib/paruly/
addpath /home/dp01/matlab/lib/boundedline/
loadspm 12
Info = LoadSubIDs;

clear Dalls
load('/imaging/dp01/results/av_nature/Dalls_avpassive_allchans.mat')
Dalls = Dalls(103:306,:,:,:);

Nsubs = sum(FileCheck);

[FileList] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');

%
D = spm_eeg_load(FileList{1});
Age = Info.Age(FileCheck);

for av = 1:2
        
    Dallrs = reshape(squeeze(Dalls(:,:,av,:)),size(Dalls,1),size(Dalls,2)*size(Dalls,4));
    % Create another PCA dataset with both grad channels
    [ug,lg,GrdsBoth] = spm_svd(Dallrs);
    ug = full(ug(:,1:10));
    GrdsBoth = GrdsBoth(:,1:10);
    
    if av == 1
        save('/imaging/dp01/camcan/meg/avpassive/PCAcomps_audio.mat','ug','lg')
    end
    
    if av == 2
        save('/imaging/dp01/camcan/meg/avpassive/PCAcomps_visual.mat','ug','lg')
    end
    
    % Variance explained by each component
    lg = lg.^2;
    lgsum = full(sum(diag(lg)));
    varexp = full(diag(lg)./lgsum);
    disp(varexp(1:5))
    
    GrdsBoth = full(GrdsBoth)';
    GrdsBoth = reshape(GrdsBoth,10,601,Nsubs);
    ug = reshape(full(ug),2,102,10);
    ug = sqrt(sum(ug.^2,1));
    ug = reshape(ug,102,10);
    disp('Done')
    
    % Create Subgroups Based on age.
    
    
    clear Bins AgeBins S GrdsBothSG
    
    Nbins = Nsubs;
    
    % Average age bins into new matrix.
    
    % Otherwise just copy data
    GrdsBothSG = GrdsBoth;
    AgeBins = Age;
    
    
    % Sorted With AgeBins
    figure(200)
    set(gcf,'position',[26 419 1380 660],'color','w');
    SortingVar = AgeBins;
    GaussConv = 5;
    lx = [0 0;100 100;200 200; 300 300];
    ly = [18 88;18 88; 18 88; 18 88];
    timevec = linspace(-100,500,601);
    xlims = [-50 400];
    
    for ii = 1:5
        pc = shiftdim(ug(:,ii)');
        Cindex = D.coor2D(D.indchantype('MEGPLANAR'));
        Clabels = D.chanlabels(D.indchantype('MEGPLANAR'));
        Cindex = Cindex(:,1:2:end);
        Clabels = Clabels(1:2:end);
        [IM,f] = spm_eeg_plotScalpData(pc,Cindex,Clabels);close
        
        subplot(3,5,ii)
        IMmask = repmat(~isnan(IM),[1 1 3]);
        clims = [-quantile(IM(:), 0.999) quantile(IM(:), 0.999)];
        IMRGB = dp_Convert2DToColor(IM,clims,paruly);
        IMRGB(IMmask == 0) = 1;

        imagesc(flipdim(IMRGB,1));
        
   
        caxis(clims)
        axis off
        
        subplot(3,5,ii+5)
        pl = plot(timevec,squeeze(mean(GrdsBothSG(ii,:,:),3)));
        
        xlim(xlims);
        set(pl,'color',[0 160 255]/255,'linewidth',2);
        xlabel('Time (ms)');ylabel('Amplitude (A.U)')
        set(gca, 'XTick',[0, 100, 200, 300, 400])
        ct = get(gca,'YTick');
        set(gca, 'YTick',[min(ct), 0, max(ct)])
        set(gca,'YtickLabel',{'-', '0', '+'})
        
        subplot(3,5,ii+10)
        Np = size(GrdsBothSG,2);
        t = timevec';
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
        set(gca,'ydir','normal','xlim',xlims)
        caxis(quantile(imdat(:),[0.001 0.999]))
        set(gca, 'XTick',[0, 100, 200, 300, 400])
        ln = line(lx',ly');
        set(ln,'Color',[0 0 0],'LineStyle','--')
        xlabel('Time (ms)');ylabel('Age (years)')
    end
    
    if av == 2
        export_fig('/imaging/dp01/results/av_nature/Figures/PCA_ERP_Visual_Paper_parula.pdf','-pdf')
    else
        export_fig('/imaging/dp01/results/av_nature/Figures/PCA_ERP_Audio_Paper_parula.pdf','-pdf')
    end
    GrdsBothAV(:,:,:,av) = GrdsBothSG;
    
    
end

save('/imaging/dp01/results/av_nature/GrdsBothAV.mat','GrdsBothAV')


%% Load Presaved data (if data is already saved, start from here)
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

[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');


Info = LoadSubIDs;

% Select Audio 1 or Visual 2
lab = {'Aud' 'Vis'};
Ncomps = 1;
T0 = 50;


for av = 1:2
    % ERP Fits
    clear Fit FitTimeSeries
    close all
    BadData = zeros(size(GrdsBothAV,3),1);
    BaseLine = zeros(Ncomps,708);
    AdjTempM = zeros(Ncomps,708,size(GrdsBothAV,2));
    clear BL AdjTemp
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,:,av),3));
        for ii = 1:size(GrdsBothAV,3)
            x = squeeze(GrdsBothAV(pci,:,ii,av)); % use original data
            showplots = 0; % 3 = show only final fit
            F = ERPfit(x', y, linspace(-100, 500, length(y)), T0, [0 400], showplots, 1e-6);
            Fit(pci,ii,:) = [F.shift F.stretch F.amp F.error sqrt(F.R2) F.Nsteps];
            BL(pci,ii) = F.baseline;
            AdjTemp(pci,ii,:) = F.AdjustedTemplate;
            title(['Subject ', Info.SubCCIDc{ii}, ' index = ', num2str(ii), ' PCI = ', num2str(pci), 'R = ', num2str(Fit(pci,ii,5)), 'Amp = ', num2str(Fit(pci,ii,3),'%2.4f')])
            if showplots == 3
                export_fig(['/imaging/dp01/results/av_nature/ERPFittingChecks' lab{av} '-subid-' num2str(ii,'%2.3d') '.bmp'],'-bmp')
            end
            disp(str(pci, ' ', ii, ' Fit = Delay: ', Fit(pci,ii,1), ' Stretch: ', Fit(pci,ii,2), ' R2: ', Fit(pci,ii,5)^2, ' Nsteps: ', Fit(pci,ii,6)))
        end
        
        % Save Data
        if av == 1
            SaveValsAs = '/imaging/dp01/results/av_nature/Fits_AVpassiveAudio.mat';
        else
            SaveValsAs = '/imaging/dp01/results/av_nature/Fits_AVpassiveVisual.mat';
        end
        
        FitOut = zeros(Ncomps,708,6);
        
        BaseLine(pci,FindFiles) = BL;
        AdjTempM(pci,FindFiles,:) = AdjTemp;
        FitOut(pci,FindFiles,:) = squeeze(Fit);
        FitOutMat = FitOut;
        FitOutCell{pci} = mat2dataset(squeeze(FitOut(1,:,:)),'VarNames',{'Latency' 'Stretch' 'Scaling' 'RMSE' 'R' 'NstepsFit'},'ObsNames',Info.SubCCIDc);
    end
    
    FitOut = FitOutCell{1}; % if there are any scripts will using FitOut, then they will get the correct PC (rather than the last one).
    save(SaveValsAs,'FitOutCell','FitOut','FindFiles','y','FitOutMat','BaseLine','AdjTempM')
    
end















