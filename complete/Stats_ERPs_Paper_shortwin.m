clear all
close all
clc

addpath /imaging/camcan/QueryFunction/QueryFun_v1/
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))
addpath /home/dp01/matlab/lib/fisherstats
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath /home/dp01/matlab/lib/DelayProjAnalysis/
addpath /imaging/dp01/scripts/av_nature/
addpath /imaging/dp01/toolboxes/mni2fs/
addpath /home/dp01/matlab/lib/
addpath /home/dp01/matlab/export_fig/

 
% 
Info = LoadSubIDs;

ERPvis  = load('/imaging/dp01/results/av_nature/Fits_AVpassiveVisual_shortwin.mat');
ERPaud  = load('/imaging/dp01/results/av_nature/Fits_AVpassiveAudio_shotwin.mat');

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


%%

load('/imaging/dp01/results/av_nature/GrdsBothAV.mat')
[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');

age = Info.Age(FileCheck);

close all

age_low = [18 44 65];
age_high =[43 64 88];

cols = {'r' 'g' 'b'};

timex = linspace(-100,500,601);

for av = 1:2
    close all
    for ii = 1:length(age_low)
        tc = squeeze(mean(GrdsBothAV(1,:,age > age_low(ii) & age < age_high(ii),av),3));
        %     tc = tc ./ std(tc);
        plot(timex, tc, cols{ii},'LineWidth',1.5)
        hold on
    end
    xlim([-50 400])
    xlabel('Time (ms)')
    ylabel('Amplitude (A.U)')
    set(gcf, 'position', [560 720 547 223],'color','w')
    export_fig(sprintf('/imaging/dp01/results/av_nature/Figures_0to200ms/Age_Group_ERFs_%d',av),'-pdf')
end

%% Plot Robust Regression of Delay Parameters

AudStr = dp_removeoutliers(ERPaud.FitOut.Stretch);
AudAmp = dp_removeoutliers(ERPaud.FitOut.Scaling);
AudLat = dp_removeoutliers(ERPaud.FitOut.Latency);
VisStr = dp_removeoutliers(ERPvis.FitOut.Stretch);
VisLat = dp_removeoutliers(ERPvis.FitOut.Latency);
VisAmp = dp_removeoutliers(ERPvis.FitOut.Scaling);

[AudStr_temp, AudLat_temp, ind] = complete(AudStr, AudLat);
AudStr(ind == 0) = NaN;
AudLat(ind == 0) = NaN;
[VisStr_temp, VisLat_temp, ind] = complete(VisStr, VisLat);
VisStr(ind == 0) = NaN;
VisLat(ind == 0) = NaN;
clear *_temp

Age = Info.Age;

Style = [];
Style.Scatter.LineStyle = 'none';
Style.Scatter.MarkerSize = 12;
Style.Scatter.Marker = '.';
Style.Scatter.Color = [0 126 189]/255;
Style.Line.LineStyle = '-';
Style.Line.Color = [189 63 0]/255;
Style.Line.LineWidth = 3;
Style.RemoveOutliers = 0;
Style.RobustOpts = 'on';
Style.ShowEquation = true;
Style.eqnformat = 'y=%.0f [%.0f,%.0f] %+2.2fx [%2.2f,%2.2f]';

close all
fig(1,'color','w','position',[200 200 1050 640])

X = Age;

subplot(2,3,1)
y = AudLat;
mdl = dp_plotFitRegression(X,y,Style);
ci = mdl.coefCI;
fprintf('AudLat Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
ylabel('Constant Delay (ms)')
xlabel('Age (years)')
axis square
line('XData',[10 96],'ydata',[0 0],'LineStyle','--','LineWidth',2)
ylim([-25 25])
xlim([10 96])

subplot(2,3,2)
y = AudStr;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('AudStr Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
hold on
ylabel('Cumulative Delay (%)')
xlabel('Age (years)')
axis square
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
ylim([60 140])
xlim([10 96])

subplot(2,3,3)
y = AudAmp;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('AudAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])

subplot(2,3,4)
y = VisLat;
mdl = dp_plotFitRegression(X,y,Style);
ci = mdl.coefCI;
fprintf('VisLat Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Constant Delay (ms)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[0 0],'LineStyle','--','LineWidth',2)
ylim([-65 +65])
xlim([10 96])

subplot(2,3,5)
y = VisStr;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('VisStr Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Cumulative Delay (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
ylim([30 170])
xlim([10 96])

subplot(2,3,6)
y = VisAmp;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('VisAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])

export_fig('/imaging/dp01/results/av_nature/Figures_0to200ms/Age_vs_Delay_Robust_Passive_shotwin.pdf','-pdf')
!cp /imaging/dp01/results/av_nature/Figures_0to200ms/Age_vs_Delay_Robust_Passive_0_150ms.pdf /home/dp01/results/av_nature/Figures/
