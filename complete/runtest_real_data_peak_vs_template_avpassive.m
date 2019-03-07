%% Load Presaved data (if data is already saved, start from here)
%#ok<*SAGROW>

clear all
close all
clc

load('/imaging/dp01/results/av_nature/GrdsBothAV.mat')
addpath /home/dp01/matlab/lib
addpath /imaging/camcan/sandbox/projects/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/toolboxes/FastICA_25/
addpath /home/dp01/matlab/lib/fitgaussian/
addpath /home/dp01/matlab/export_fig/
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))
addpath /imaging/dp01/scripts/ERPfitting
addpath /imaging/camcan/QueryFunction/QueryFun_v1/
addpath /imaging/dp01/scripts/av_nature/code_tests/

[FileList, FileCheck, FindFiles] = CCQuery_QuickCheck('/imaging/dp01/cc700meg-mf22-spm12/release002_transdef/avpassive/<CCID>/efMspm12_transdef_transrest_mf2pt2_passive_raw.mat');

I = LoadSubIDs;

age = I.Age(FileCheck);


%% Standard

clear FitStruct;

Info = LoadSubIDs;
% Select Audio 1 or Visual 2
FitCentres = [50 50];
lab = {'Aud' 'Vis'};


Ncomps = 1;
for av = 1:2
    % ERP Fits
    BadData = zeros(size(GrdsBothAV,3),1);
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,:,av),3));
        for ii = 1:size(GrdsBothAV,3)
            x = squeeze(GrdsBothAV(pci,:,ii,av)); % use original data
            xmin = -100; % t minimum of data
            xmax = 500; % t maximum of data
            FitCentre = FitCentres(av);
            Xrange = [0 400]; % samples
            difflim = 1e-6; % tolerance limit for gradient ascent (correlation value)
            showplots = 0; % 3 = show only final fit
%             Fit(pci,av,ii,:) = ERPfitting_v1(x,y,xmin,xmax,FitCentre,101:500,showplots,difflim,[2 2]); % Fit = [delay stretch];

            FitStruct(pci,ii,av) = ERPfit(x', y, linspace(xmin, xmax, length(y)), FitCentre, Xrange, showplots, 1e-6); % Fit = [delay stretch];
            
            if showplots == 3
                export_fig(['/imaging/dp01/results/proj2trans/ERPFittingChecks' lab{av} '-subid-' num2str(ii,'%2.3d') '.bmp'],'-bmp')
            end
            disp(str(pci, ' ', ii,' Fit = Delay: ',FitStruct(pci,ii).shift,' Stretch: ',FitStruct(pci,ii).stretch, ' R2: ',FitStruct(pci,ii).R2^2,' Nsteps: ',FitStruct(pci,ii).Nsteps))
        end
    end
end

%% Young template
Ncomps = 1;

ageCutoff = 30;

age = I.Age(FileCheck);

templateIndices = age < ageCutoff;
testIndicesYoung = age > 0;

ageTest = age(testIndicesYoung);

for av = 1:2
    % ERP Fits
    close all
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,templateIndices,av),3));
        GrdsBothAVtest = GrdsBothAV(pci,:,testIndicesYoung,:);

        for ii = 1:size(GrdsBothAVtest,3)
            x = squeeze(GrdsBothAVtest(pci,:,ii,av)); % use original data
            xmin = -100; % t minimum of data
            xmax = 500; % t maximum of data
            FitCentre = FitCentres(av);
            Xrange = [0 400]; % samples
            difflim = 1e-6; % tolerance limit for gradient ascent (correlation value)
            showplots = 0; % 3 = show only final fit

            FitStructYoungTemplate(pci,ii,av) = ERPfit(x', y, linspace(xmin, xmax, length(y)), FitCentre, Xrange, showplots, 1e-6); % Fit = [delay stretch];
            
            disp(str(pci, ' ', ii,' Fit = Delay: ',FitStruct(pci,ii).shift,' Stretch: ',FitStruct(pci,ii).stretch, ' R2: ',FitStruct(pci,ii).R2^2,' Nsteps: ',FitStruct(pci,ii).Nsteps))
        end
    end
end


%% Old Template

Ncomps = 1;

ageCutoff = 60;

age = I.Age(FileCheck);

templateIndices = age > ageCutoff;
testIndicesOld = age > 0;

ageTest = age(testIndicesOld);

for av = 1:2
    % ERP Fits
    close all
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,templateIndices,av),3));
        GrdsBothAVtest = GrdsBothAV(pci,:,testIndicesOld,:);

        for ii = 1:size(GrdsBothAVtest,3)
            x = squeeze(GrdsBothAVtest(pci,:,ii,av)); % use original data
            xmin = -100; % t minimum of data
            xmax = 500; % t maximum of data
            FitCentre = FitCentres(av);
            Xrange = [0 400]; % samples
            difflim = 1e-6; % tolerance limit for gradient ascent (correlation value)
            showplots = 0; % 3 = show only final fit

            FitStructOldTemplate(pci,ii,av) = ERPfit(x', y, linspace(xmin, xmax, length(y)), FitCentre, Xrange, showplots, 1e-6); % Fit = [delay stretch];
            
            disp(str(pci, ' ', ii,' Fit = Delay: ',FitStruct(pci,ii).shift,' Stretch: ',FitStruct(pci,ii).stretch, ' R2: ',FitStruct(pci,ii).R2^2,' Nsteps: ',FitStruct(pci,ii).Nsteps))
        end
    end
end


save('/imaging/dp01/results/av_nature/Old_vs_Young_Templates.mat','FitStruct*','FileCheck','FindFiles')

%% Plot young and old templates

load /imaging/dp01/results/av_nature/Old_vs_Young_Templates.mat

AllFits(1,:,:,:) = FitStructYoungTemplate;
AllFits(2,:,:,:) = FitStructOldTemplate;


for agei = 1:2
    
    Age = I.Age;
    
    AudStr = nan(708,1);
    AudAmp = nan(708,1);
    AudLat = nan(708,1);
    VisStr = nan(708,1);
    VisAmp = nan(708,1);
    VisLat = nan(708,1);
    
    zeroind = find([AllFits(agei,1,:,1).amp] == 0);
    
    for ii = 1:length(zeroind)
        FitStructOldTemplate(:,zeroind(ii),:).stretch = NaN;
    end
    
    AudStr(FileCheck) = dp_removeoutliers([AllFits(agei,1,:,1).stretch]');
    AudAmp(FileCheck) = dp_removeoutliers([AllFits(agei,1,:,1).amp]');
    AudLat(FileCheck) = dp_removeoutliers([AllFits(agei,1,:,1).shift]');
    VisStr(FileCheck) = dp_removeoutliers([AllFits(agei,1,:,2).stretch]');
    VisAmp(FileCheck) = dp_removeoutliers([AllFits(agei,1,:,2).amp]');
    VisLat(FileCheck) = dp_removeoutliers([AllFits(agei,1,:,2).shift]');
    
    [AudStr_temp, AudLat_temp, ind] = complete(AudStr, AudLat);
    AudStr(ind == 0) = NaN;
    AudLat(ind == 0) = NaN;
    [VisStr_temp, VisLat_temp, ind] = complete(VisStr, VisLat);
    VisStr(ind == 0) = NaN;
    VisLat(ind == 0) = NaN;
   
    
    clear *_temp
    
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
    
    if agei == 1
        export_fig('/imaging/dp01/results/av_nature/Figures/ERF_Fitting_Young_Template','-pdf')
    else
        export_fig('/imaging/dp01/results/av_nature/Figures/ERF_Fitting_Old_Template','-pdf')
    end
    
end


%% First Test Peak Finding Script

close all


timex = linspace(-100,500,601);
peakO = [100, 200, 300, ]; % original peak locations
peakA = [-1, 0.5, -1];
peakT = peakO;
peakW = [25, 35, 100];
P(2).peakwin = {
    [-50 50] 
    [-50 50]
    [-50 50]
    };
peaksign = [-1 1 -1];

template = createERPTemplate(timex, peakA, peakW, peakT);
plmin = min(template);
plmax = max(template);

testpeaks = findpeaksSimple(template, P(2).peakwin, peakO, peaksign, timex);

plot(timex, template)
hold on

for ii = 1:length(testpeaks); 
   line('XData',[testpeaks(ii); testpeaks(ii)],'YData',[plmin; plmax])
end

title('Peak Locations')


%% Plot the mean ERFs 
close all
pci = 1

figure
plot(timex, squeeze(mean(GrdsBothAV(pci,:,:,1),3)));

figure
plot(timex, squeeze(mean(GrdsBothAV(pci,:,:,2),3)));


%% Peak finding

Ncomps = 1;

age = I.Age(FileCheck);

timex = linspace(-100, 500, 601);

P = [];
window = [-60 60];
P(1).peakzero = [54 106 198];
P(1).peaksign = [1 -1 1];
P(1).peakwin = {
    window
    window
    window
    };

window = [-60 60];
P(2).peakzero = [112 210];
P(2).peaksign = [-1 1];
P(2).peakwin = {
    window
    window
    };

for av = 1:2
    
    for pci = 1:Ncomps
        y = squeeze(mean(GrdsBothAV(pci,:,:,av),3));
        GrdsBothAVtest = GrdsBothAV(pci,:,:,:);

        for ii = 1:size(GrdsBothAVtest,3)
            x = squeeze(GrdsBothAVtest(pci,:,ii,av)); % use original data
            P(av).Peaks(ii, pci, :) = findpeaksSimple(x, P(av).peakwin, P(av).peakzero, P(av).peaksign, timex);
            disp(str(pci, ' ', ii,' Fit = Delay: ',FitStruct(pci,ii).shift,' Stretch: ',FitStruct(pci,ii).stretch, ' R2: ',FitStruct(pci,ii).R2^2,' Nsteps: ',FitStruct(pci,ii).Nsteps))
        end
    end
end



% Plot Peak Latencies
close all

Settings.Scatter.LineStyle = 'none';
Settings.Scatter.MarkerSize = 12;
Settings.Scatter.Marker = '.';
Settings.Scatter.Color = [0 126 189]/255;
Settings.Line.LineStyle = '-';
Settings.Line.Color = [189 63 0]/255;
Settings.Line.LineWidth = 3;
Settings.RemoveOutliers = 0;
Settings.ShowEquation = true;
Settings.RobustOpts = 'on';
Style.eqnformat = 'y=%.0f [%.0f,%.0f] %+2.2fx [%2.2f,%2.2f]';

for av = 1:2
    pc = 1;
    for pi = 1:length(P(av).peakzero)
        fig(av);
        set(gcf,'position',[3 265 300 678],'color','white')
        subplot(3,1,pi)
        p = P(av).Peaks(:,pc,pi);
        win = P(av).peakzero(pi)+P(av).peakwin{pi};
        p(p == win(1)) = NaN;
        p(p == win(2)) = NaN;
        P(av).Peaks(:,pc,pi) = p;
        mdl = dp_plotFitRegression(age, p, Settings);
        ci = mdl.coefCI;
        
        fprintf('%d %d Intercept: CI: %2.2f %2.2f \n', av, pi, ci(1,:));
        fprintf('%d %d Slope: CI: %2.2f %2.2f \n', av, pi, ci(2,:));
        
        ylabel('Peak Latency (ms)');
        axis tight
        xlim([18 88])
    end
    export_fig(sprintf('/imaging/dp01/results/av_nature/Figures/Peak_Fitting_Real_Data%d',av),'-pdf')
end


%% Plot peak to peak latency
close all
Settings.Scatter.LineStyle = 'none';
Settings.Scatter.MarkerSize = 12;
Settings.Scatter.Marker = '.';
Settings.Scatter.Color = [0 126 189]/255;
Settings.Line.LineStyle = '-';
Settings.Line.Color = [189 63 0]/255;
Settings.Line.LineWidth = 3;
Settings.RemoveOutliers = 0;
Settings.ShowEquation = true;
Settings.RobustOpts = 'on';
Style.eqnformat = 'y=%.0f [%.0f,%.0f] %+2.2fx [%2.2f,%2.2f]';

fig(1, 'position',[560 330 670 613])
subplot(2,2,1)
n1p1 = squeeze(P(2).Peaks(:,:,2) - P(2).Peaks(:,:,1));
[mdl, tl] = dp_plotFitRegression(age, n1p1, Style);
ylabel('Vis. Peak2-Peak1 (ms)')
xlabel('Age (years)')
xlim([18 88])

subplot(2,2,2)
n1p1 = squeeze(P(1).Peaks(:,:,2) - P(1).Peaks(:,:,1));
n1p1(n1p1 < 0) = NaN;
[mdl, tl] = dp_plotFitRegression(age, n1p1, Style);
ylabel('Aud. Peak2-Peak1 (ms)')
xlabel('Age (years)')
xlim([18 88])

subplot(2,2,4)
n1p1 = squeeze(P(1).Peaks(:,:,3) - P(1).Peaks(:,:,2));
[mdl, tl] = dp_plotFitRegression(age, n1p1, Style);
ylabel('Aud. Peak3-Peak2 (ms)')
xlabel('Age (years)')
xlim([18 88])

export_fig('/imaging/dp01/results/av_nature/Figures/peak-to-peak_latencies', '-pdf')

%% Derive cumulative delay from peak data.
close all

Peaks = squeeze(P(av).Peaks(:,pc,:));
PeaksDiff = diff(Peaks,1,2);

% Do not need a reference here, since the slope reference is always const =
% 0 and slope = 1;

for av = 1
    disp(av)
    P(av).intercept = zeros(size(P(av).Peaks,1),1);
    P(av).slope = P(av).intercept;
    for ii = 1:size(P(av).Peaks,1)
        if all(~isnan(P(av).Peaks(ii,1,:)))
            l = LinearModel.fit(P(av).peakzero, squeeze(P(av).Peaks(ii,1,:)));
            P(av).intercept(ii) = l.Coefficients{'(Intercept)','Estimate'};
            P(av).slope(ii) = l.Coefficients{'x1','Estimate'};
        else
            P(av).intercept(ii) = NaN;
            P(av).slope(ii) = NaN;
        end
    end
end

pcumtv = diff(P(2).peakzero) ./ squeeze(P(2).Peaks(:,1,2)-P(2).Peaks(:,1,1));
pconst = P(2).Peaks(:,1,1) - P(2).peakzero(1).*pcumtv;


%%
close all

Settings.Scatter.LineStyle = 'none';
Settings.Scatter.MarkerSize = 12;
Settings.Scatter.Marker = '.';
Settings.Scatter.Color = [0 126 189]/255;
Settings.Line.LineStyle = '-';
Settings.Line.Color = [189 63 0]/255;
Settings.Line.LineWidth = 3;
Settings.RemoveOutliers = 0;
Settings.ShowEquation = true;
Settings.RobustOpts = 'on';
Style.eqnformat = 'y=%.0f [%.0f,%.0f] %+2.2fx [%2.2f,%2.2f]';

fig(1,'color','w')
subplot(2,2,1)
mdl = dp_plotFitRegression(age, dp_removeoutliers(P(1).intercept), Settings);
ci = mdl.coefCI;
fprintf('%d %d Intercept: CI: %2.2f %2.2f \n', av, pi, ci(1,:));
fprintf('%d %d Slope: CI: %2.2f %2.2f \n', av, pi, ci(2,:));

ylabel('Auditory Const. Delay (ms)')
xlabel('Age (years)')
xlim([18 88])

subplot(2,2,2)
mdl = dp_plotFitRegression(age, dp_removeoutliers(pconst), Settings);
ci = mdl.coefCI;
fprintf('%d %d Intercept: CI: %2.2f %2.2f \n', av, pi, ci(1,:));
fprintf('%d %d Slope: CI: %2.2f %2.2f \n', av, pi, ci(2,:));
ylabel('Visual Const. Delay (ms)')
xlim([18 88])
xlabel('Age (years)')

subplot(2,2,3)
mdl = dp_plotFitRegression(age, dp_removeoutliers(P(1).slope)*100, Settings);
ci = mdl.coefCI;
fprintf('%d %d Intercept: CI: %2.2f %2.2f \n', av, pi, ci(1,:));
fprintf('%d %d Slope: CI: %2.2f %2.2f \n', av, pi, ci(2,:));
ylabel('Auditory Cumtv. Delay (%)')
xlim([18 88])
xlabel('Age (years)')

subplot(2,2,4)
mdl = dp_plotFitRegression(age, dp_removeoutliers(pcumtv)*100, Settings);
ci = mdl.coefCI;
fprintf('%d %d Intercept: CI: %2.2f %2.2f \n', av, pi, ci(1,:));
fprintf('%d %d Slope: CI: %2.2f %2.2f \n', av, pi, ci(2,:));
ylabel('Visual Cumtv. Delay (%)')
xlim([18 88])
xlabel('Age (years)')

export_fig('/imaging/dp01/results/av_nature/Figures/Peak_Fitting_Real_Data_concum','-pdf')

%% Check Fit

close all

FitShift = squeeze(Fit(1,:,1))';
FitStretch = squeeze(Fit(1,:,2))';

teststretch = [FitStructOldTemplate(1,:,2).stretch]';
testshift = [FitStructOldTemplate(1,:,2).shift]';

dp_plotFit(testshift,dp_removeoutliers(FitShift),1 )

figure 
hist(testshift)

figure
hist(FitShift)