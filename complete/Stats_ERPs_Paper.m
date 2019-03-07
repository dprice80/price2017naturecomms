clear all
close all
clc

addpath /imaging/camcan/QueryFunction/QueryFun_v1/
addpath /home/dp01/matlab/PolyfitnTools/
addpath(genpath('/home/dp01/matlab/mediation_toolbox/'))
addpath /home/dp01/matlab/lib/fisherstats
addpath(genpath('/home/dp01/matlab/SCN_Core_Support/'))
addpath /home/dp01/matlab/lib/DelayProjAnalysis/
addpath /imaging/dp01/scripts/av_nature/
addpath /imaging/dp01/toolboxes/mni2fs/
addpath /home/dp01/matlab/lib/
addpath /home/dp01/matlab/export_fig/
addpath /home/dp01/matlab/lib/paruly/

DAT = [];
DAT.SessionList = {
    'AVtaskRTs'           '/imaging/camcan/cc700-scored/MEG/release001/data/<CCID>/MEG*scored.txt'
    };

DAT = CCQuery_CheckFiles(DAT);
DAT.dataflag = 1;
DAT = CCQuery_LoadData(DAT);
% 
Info = LoadSubIDs;
Gender = Info.GenderNum;
hand = Info.AddData.handedness;

RT = DAT.dataset.AVtaskRTs;

Age = Info.Age;

ERPvis  = load('/imaging/dp01/results/av_nature/Fits_AVpassiveVisual.mat');
ERPaud  = load('/imaging/dp01/results/av_nature/Fits_AVpassiveAudio.mat');

VisBL = ERPvis.BaseLine;
AudBL = ERPaud.BaseLine;

VisAdjTemp = ERPvis.AdjTempM;
AudAdjTemp = ERPaud.AdjTempM;

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

[AudStr_temp, AudLat_temp, ind] = complete(AudStr, AudLat);
AudStr(ind == 0) = NaN;
AudLat(ind == 0) = NaN;
AudAmp(ind == 0) = NaN;
AudBL(ind == 0) = NaN;
[VisStr_temp, VisLat_temp, ind] = complete(VisStr, VisLat);
VisStr(ind == 0) = NaN;
VisLat(ind == 0) = NaN;
VisAmp(ind == 0) = NaN;
VisBL(ind == 0) = NaN;

AudBL = dp_removeoutliers(AudBL);
VisBL = dp_removeoutliers(VisBL);

%%%%%%%%%%%%%%%%% AV TASK %%%%%%%%%%%%%%%%%%%%%
TaskERPvis  = load('/imaging/dp01/results/av_nature/Fits_AVtaskVis_ApplyPassive.mat');
TaskERPaud  = load('/imaging/dp01/results/av_nature/Fits_AVtaskAud_ApplyPassive.mat');

% Remove zeros (these are just subjects that had no data)
ind = TaskERPvis.FitOut.NstepsFit == 0;
TaskERPvis.FitOut(ind, :) = dataset(NaN);

% Remove zeros
ind = TaskERPaud.FitOut.NstepsFit == 0;
TaskERPaud.FitOut(ind, :) = dataset(NaN);

% Remove cases where amplitude was less than 0
ind = TaskERPaud.FitOut.Scaling < 0;
TaskERPaud.FitOut(ind, :) = dataset(NaN);

ind = TaskERPvis.FitOut.Scaling < 0;
TaskERPvis.FitOut(ind, :) = dataset(NaN);


TaskVisLat = dp_removeoutliers(TaskERPvis.FitOut.Latency);
TaskVisStr = dp_removeoutliers(TaskERPvis.FitOut.Stretch);
TaskAudLat = dp_removeoutliers(TaskERPaud.FitOut.Latency);
TaskAudStr = dp_removeoutliers(TaskERPaud.FitOut.Stretch);
TaskAudAmp = dp_removeoutliers(TaskERPaud.FitOut.Scaling);
TaskVisAmp = dp_removeoutliers(TaskERPvis.FitOut.Scaling);

[~, ~, ind] = complete(TaskAudStr, TaskAudLat);
TaskAudStr(ind == 0) = NaN;
TaskAudLat(ind == 0) = NaN;
TaskAudAmp(ind == 0) = NaN;
[~, ~, ind] = complete(TaskVisStr, TaskVisLat);
TaskVisStr(ind == 0) = NaN;
TaskVisLat(ind == 0) = NaN;
TaskVisAmp(ind == 0) = NaN;


% Auditory and visual acuity
Dem = dataset('File','/imaging/dp01/results/proj1/measures/Subset_Vars_From_Karen.txt');
audac = Dem.proportion_tones_heard;
audac(audac > 1) = 1;
visac = Dem.vision_both_eyes;


%% Check whether latency measures are correlated 

close all
figure
dp_plotFitRegression(AudStr, VisLat)
figure
dp_plotFitRegression(separateSignal(Age, AudStr), VisLat)


%% Check RMSE 

rmsevis = log(ERPvis.FitOut.RMSE);
rmseaud = log(ERPaud.FitOut.RMSE);

dp_plotFitRegression(Age, rmsevis);
dp_plotFitRegression(Age, rmseaud);


%% Plot Amplitude Scaling and Mean Power 

VisLat = dp_removeoutliers(ERPvis.FitOut.Latency);
VisStr = dp_removeoutliers(ERPvis.FitOut.Stretch);
AudLat = dp_removeoutliers(ERPaud.FitOut.Latency);
AudStr = dp_removeoutliers(ERPaud.FitOut.Stretch);
AudAmp = dp_removeoutliers(ERPaud.FitOut.Scaling);
VisAmp = dp_removeoutliers(ERPvis.FitOut.Scaling);

VisBL = dp_removeoutliers(ERPvis.BaseLine);
AudBL = dp_removeoutliers(ERPaud.BaseLine);
% VisBL = (ERPvis.BaseLine);
% AudBL = (ERPaud.BaseLine);

[AudStr_temp, AudLat_temp, ind] = complete(AudStr, AudLat);
AudStr(ind == 0) = NaN;
AudLat(ind == 0) = NaN;
AudAmp(ind == 0) = NaN;
AudBL(ind == 0) = NaN;
[VisStr_temp, VisLat_temp, ind] = complete(VisStr, VisLat);
VisStr(ind == 0) = NaN;
VisLat(ind == 0) = NaN;
VisAmp(ind == 0) = NaN;
VisBL(ind == 0) = NaN;

% Export data for external users
T = table();
T.AudCum = AudStr;
T.AudCon = AudLat;
T.AudAmp = AudAmp;
T.AudCum = AudStr;
T.AudCon = AudLat;
T.AudAmp = AudAmp;
T.Age = Age;
writetable(T, '/imaging/dp01/results/av_nature/export/sharing/latency_values_anon.csv');

Settings.Scatter.LineStyle = 'none';
Settings.Scatter.MarkerSize = 12;
Settings.Scatter.Marker = '.';
Settings.Scatter.Color = [0 126 189]/255;
Settings.Line.LineStyle = '-';
Settings.Line.Color = [189 63 0]/255;
Settings.Line.LineWidth = 3;
Settings.RemoveOutliers = 0;
Settings.ShowEquation = false;
Settings.RobustOpts = 'on';
    
close all
fig(1, 'position',[200 102 685 698])
subplot(2,2,1)
dp_plotFitRegression(Age, AudBL*100000, Settings)
xlabel('Age (yrs)')
axis square
ylabel('Auditory Baseline Shift (A.U)')
xlim([10 96])

X = Age;

subplot(2,2,2)
mdl = dp_plotFitRegression(X,AudAmp*100,Settings);
ci = mdl.coefCI;
fprintf('AudAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])

subplot(2,2,3)
dp_plotFitRegression(Age, VisBL*100000, Settings)
xlabel('Age (yrs)')
ylabel('Visual Baseline Shift (A.U)')
axis square
xlim([10 96])

subplot(2,2,4)

mdl = dp_plotFitRegression(X,VisAmp*100,Settings);
ci = mdl.coefCI;
fprintf('VisAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])

export_fig('/imaging/dp01/results/av_nature/Figures/Age_vs_Baseline','-pdf')


%% Compare delay between tasks
clc

TaskAudStr(isnan(AudStr)) = NaN;
TaskAudLat(isnan(AudLat)) = NaN;

TaskVisStr(isnan(VisStr)) = NaN;
TaskVisLat(isnan(VisLat)) = NaN;
TaskAudLat = separateSignal(TaskAudStr', TaskAudLat')';

test = TaskAudStr;
test(isnan(AudStr)) = NaN;
dp_plotFitRegression(Age, test);
dp_plotFitRegression(Age, AudStr);

[p,~,h] = signrank(AudStr, TaskAudStr);
disp(sprintf('AudStr: z=%2.3f, p=%2.3f',h.zval,p))

[p,~,h] = signrank(AudLat, TaskAudLat);
disp(sprintf('AudLat: z=%2.3f, p=%2.3f',h.zval,p))

[p,~,h] = signrank(VisStr, TaskVisStr);
disp(sprintf('VisStr: z=%2.3f, p=%2.3f',h.zval,p))

[p,~,h] = signrank(VisLat, TaskVisLat);
disp(sprintf('VisLat: z=%2.3f, p=%2.3f',h.zval,p))

disp('Tests for differences in slope between tasks')

[r,p] = corr(Age, AudLat - TaskAudLat, 'type', 'spearman', 'rows','complete');
disp(sprintf('AudLat: Rho=%2.3f, P=%2.3f, N=%d', r, p, sum(~isnan(AudLat - TaskAudLat))))

[r,p] = corr(Age, AudStr - TaskAudStr, 'type', 'spearman', 'rows','complete');
disp(sprintf('AudStr: Rho=%2.3f, P=%2.3f, N=%d',r,p, sum(~isnan(AudStr - TaskAudStr))))

[r,p] = corr(Age, VisLat - TaskVisLat, 'type', 'spearman', 'rows','complete');
disp(sprintf('VisLat: Rho=%2.3f, P=%2.3f, N=%d',r,p, sum(~isnan(VisLat - TaskVisLat))))

[r,p] = corr(Age, VisStr - TaskVisStr, 'type', 'spearman', 'rows','complete');
disp(sprintf('VisStr: Rho=%2.3f, P=%2.3f, N=%d',r,p, sum(~isnan(VisStr - TaskVisStr))))


%% Plot Task Correlations

clear *_temp

rmsevis = log(TaskERPvis.FitOut.RMSE);
rmseaud = log(TaskERPaud.FitOut.RMSE);
% rmsevis = ones(708,1);
% rmseaud = ones(708,1);

audac_sep = separateSignal(Age,[audac]);
visac_sep = separateSignal(Age,[visac]);

% [AudStr, AudAmp, AudLat, TaskVisStr, TaskVisLat, TaskVisAmp] = dp_removeoutliers(AudStr, AudAmp, AudLat, TaskVisStr, TaskVisLat, TaskVisAmp);

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
y = TaskAudLat;
mdl = dp_plotFitRegression(X,y,Style);
ci = mdl.coefCI;
fprintf('TaskAudLat Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
ylabel('Constant Delay (ms)')
xlabel('Age (years)')
axis square
line('XData',[10 96],'ydata',[0 0],'LineStyle','--','LineWidth',2)
ylim([-25 25])
xlim([10 96])

subplot(2,3,2)
y = TaskAudStr;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TaskAudStr Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
hold on
ylabel('Cumulative Delay (%)')
xlabel('Age (years)')
axis square
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
ylim([60 140])
xlim([10 96])

subplot(2,3,3)
y = TaskAudAmp;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TaskAudAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])

subplot(2,3,4)
y = TaskVisLat;
mdl = dp_plotFitRegression(X,y,Style);
ci = mdl.coefCI;
fprintf('TaskVisLat Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Constant Delay (ms)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[0 0],'LineStyle','--','LineWidth',2)
ylim([-65 +65])
xlim([10 96])

subplot(2,3,5)
y = TaskVisStr;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TaskVisStr Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Cumulative Delay (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
ylim([30 170])
xlim([10 96])

subplot(2,3,6)
y = TaskVisAmp;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TaskVisAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])


% Replace passive data with temporary variables, which are complete sets
% across both tasks (NANs removed)

TempAudStr = AudStr;
TempAudLat = AudLat;
TempVisStr = VisStr;
TempVisLat = VisLat;

TempAudStr(isnan(TaskAudStr)) = NaN;
TempAudLat(isnan(TaskAudLat)) = NaN;
TempAudAmp = AudAmp;
TempVisStr(isnan(TaskVisStr)) = NaN;
TempVisLat(isnan(TaskVisLat)) = NaN;
TempVisAmp = VisAmp;

fig(2,'color','w','position',[200 200 1050 640])

X = Age;

subplot(2,3,1)
y = TempAudLat;
mdl = dp_plotFitRegression(X,y,Style);
ci = mdl.coefCI;
fprintf('TempAudLat Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
ylabel('Constant Delay (ms)')
xlabel('Age (years)')
axis square
line('XData',[10 96],'ydata',[0 0],'LineStyle','--','LineWidth',2)
ylim([-25 25])
xlim([10 96])

subplot(2,3,2)
y = TempAudStr;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TempAudStr Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
hold on
ylabel('Cumulative Delay (%)')
xlabel('Age (years)')
axis square
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
ylim([60 140])
xlim([10 96])

subplot(2,3,3)
y = TempAudAmp;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TempAudAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])

subplot(2,3,4)
y = TempVisLat;
mdl = dp_plotFitRegression(X,y,Style);
ci = mdl.coefCI;
fprintf('TempVisLat Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Constant Delay (ms)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[0 0],'LineStyle','--','LineWidth',2)
ylim([-65 +65])
xlim([10 96])

subplot(2,3,5)
y = TempVisStr;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TempVisStr Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Cumulative Delay (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
ylim([30 170])
xlim([10 96])

subplot(2,3,6)
y = TempVisAmp;
mdl = dp_plotFitRegression(X,y*100,Style);
ci = mdl.coefCI;
fprintf('TempVisAmp Intercept: CI: %2.2f %2.2f \n', ci(1,:));
fprintf('       Slope: CI: %2.2f %2.2f \n', ci(2,:));
axis square
ylabel('Scaling Factor (%)')
xlabel('Age (years)')
line('XData',[10 96],'ydata',[100 100],'LineStyle','--','LineWidth',2)
xlim([10 96])


%% Plot Robust Regression of Delay Parameters

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
audind = ~isnan(y);


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
visind = ~isnan(y);

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

export_fig('/imaging/dp01/results/av_nature/Figures/Age_vs_Delay_Robust_Passive.pdf','-pdf')



%% Analyse Reaction Times

close all

ind = zeros(1,708);
ind(ERPvis.FindFiles) = 1;
Age = Info.Age;
Age(ind == 0) = NaN;
gender = Info.Genderc;
genderval = Info.GenderNum;

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
C = mat2dataset(zeros(10,3));
C.Properties.ObsNames = {'age' 'page' 'audcon' 'paudcon' 'audcum' 'paudcum' 'viscon' 'pviscon' 'viscum' 'pviscum'};
C.Properties.VarNames = {'RTmeg' 'RTmeg_conage' 'N'};

RTmeg = RT.mdnRT;

bc = 1; % bonferroni correlation

% MEG RTs
ctype = 'spearman';
[R,P,N] = create_RT_table_cell(ctype,Age, RTmeg);
C('age', 'RTmeg') = dataset(R);
C('page', 'RTmeg') = dataset(P*bc);
C('age', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,AudLat, RTmeg);
C('audcon', 'RTmeg') = dataset(R);
C('paudcon', 'RTmeg') = dataset(P*bc);
C('audcon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,AudStr, RTmeg);
C('audcum', 'RTmeg') = dataset(R);
C('paudcum', 'RTmeg') = dataset(P*bc);
C('audcum', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,VisLat, RTmeg);
C('viscon', 'RTmeg') = dataset(R);
C('pviscon', 'RTmeg') = dataset(P*bc);
C('viscon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,VisStr, RTmeg);
C('viscum', 'RTmeg') = dataset(R);
C('pviscum', 'RTmeg') = dataset(P*bc);
C('viscum', 'N') = dataset(N);

% control for age.
[R,P,N] = create_RT_table_cell(ctype, Age, RTmeg);
C('age', 'RTmeg_conage') = dataset(R);
C('page', 'RTmeg_conage') = dataset(P*bc);
C('age', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,AudLat, RTmeg, Age);
C('audcon', 'RTmeg_conage') = dataset(R);
C('paudcon', 'RTmeg_conage') = dataset(P*bc);
C('audcon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,AudStr, RTmeg, Age);
C('audcum', 'RTmeg_conage') = dataset(R);
C('paudcum', 'RTmeg_conage') = dataset(P*bc);
C('audcum', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,VisLat, RTmeg, Age);
C('viscon', 'RTmeg_conage') = dataset(R);
C('pviscon', 'RTmeg_conage') = dataset(P*bc);
C('viscon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,VisStr, RTmeg, Age);
C('viscum', 'RTmeg_conage') = dataset(R);
C('pviscum', 'RTmeg_conage') = dataset(P*bc);
C('viscum', 'N') = dataset(N);

addpath(genpath('/home/dp01/matlab/MediationToolbox/mediation_toolbox/'))
addpath(genpath('/home/dp01/matlab/RobustToolbox/robust_toolbox/'))
addpath(genpath('/home/dp01/matlab/SCN_Core_Support'))

% mediation model was not significant.
[paths, stats] = mediation(Age, RTmeg, AudStr);

dp_plotFitRegression(AudLat, RTmeg)

C

%% Reaction Times (TASK)
close all

ind = zeros(1,708);
ind(ERPvis.FindFiles) = 1;
Age = Info.Age;
Age(ind == 0) = NaN;
gender = Info.Genderc;
genderval = Info.GenderNum;

close all
C = mat2dataset(zeros(10,3));
C.Properties.ObsNames = {'age' 'page' 'audcon' 'paudcon' 'audcum' 'paudcum' 'viscon' 'pviscon' 'viscum' 'pviscum'};
C.Properties.VarNames = {'RTmeg' 'RTmeg_conage' 'N'};

RTmeg = RT.mdnRT;

bc = 1; % bonferroni correlation

% MEG RTs
ctype = 'spearman';
[R,P,N] = create_RT_table_cell(ctype,Age, RTmeg);
C('age', 'RTmeg') = dataset(R);
C('page', 'RTmeg') = dataset(P*bc);
C('age', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskAudLat, RTmeg);
C('audcon', 'RTmeg') = dataset(R);
C('paudcon', 'RTmeg') = dataset(P*bc);
C('audcon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskAudStr, RTmeg);
C('audcum', 'RTmeg') = dataset(R);
C('paudcum', 'RTmeg') = dataset(P*bc);
C('audcum', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskVisLat, RTmeg);
C('viscon', 'RTmeg') = dataset(R);
C('pviscon', 'RTmeg') = dataset(P*bc);
C('viscon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskVisStr, RTmeg);
C('viscum', 'RTmeg') = dataset(R);
C('pviscum', 'RTmeg') = dataset(P*bc);
C('viscum', 'N') = dataset(N);

% control for age.
[R,P,N] = create_RT_table_cell(ctype, Age, RTmeg);
C('age', 'RTmeg_conage') = dataset(R);
C('page', 'RTmeg_conage') = dataset(P*bc);
C('age', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskAudLat, RTmeg, Age);
C('audcon', 'RTmeg_conage') = dataset(R);
C('paudcon', 'RTmeg_conage') = dataset(P*bc);
C('audcon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskAudStr, RTmeg, Age);
C('audcum', 'RTmeg_conage') = dataset(R);
C('paudcum', 'RTmeg_conage') = dataset(P*bc);
C('audcum', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskVisLat, RTmeg, Age);
C('viscon', 'RTmeg_conage') = dataset(R);
C('pviscon', 'RTmeg_conage') = dataset(P*bc);
C('viscon', 'N') = dataset(N);

[R,P,N] = create_RT_table_cell(ctype,TaskVisStr, RTmeg, Age);
C('viscum', 'RTmeg_conage') = dataset(R);
C('pviscum', 'RTmeg_conage') = dataset(P*bc);
C('viscum', 'N') = dataset(N);

addpath(genpath('/home/dp01/matlab/MediationToolbox/mediation_toolbox/'))
addpath(genpath('/home/dp01/matlab/RobustToolbox/robust_toolbox/'))
addpath(genpath('/home/dp01/matlab/SCN_Core_Support'))

% mediation model was not significant.
[paths, stats] = mediation(Age, RTmeg, AudStr);


C


%% visual and auditory acuity fits using spearman's partial correlations 

AudStr = dp_removeoutliers(ERPaud.FitOut.Stretch);
AudAmp = dp_removeoutliers(ERPaud.FitOut.Scaling);
AudLat = dp_removeoutliers(ERPaud.FitOut.Latency);
VisStr = dp_removeoutliers(ERPvis.FitOut.Stretch);
VisLat = dp_removeoutliers(ERPvis.FitOut.Latency);
VisAmp = dp_removeoutliers(ERPvis.FitOut.Scaling);

AudStr(audind == 0) = NaN;
AudLat(audind == 0) = NaN;
AudAmp(audind == 0) = NaN;

VisStr(visind == 0) = NaN;
VisLat(visind == 0) = NaN;
VisAmp(visind == 0) = NaN;


[audac_nan,~,Naud] = dp_complete_nan(audac, AudStr);
Naud = sum(Naud);

[visac_nan,~,Nvis] = dp_complete_nan(visac, VisLat);
Nvis = sum(Nvis);

audconvar = audac;
visconvar = visac;

spearmans_Pcorr_acuity_RMSE_amp

%% Correcting for visual / auditory AMPLITUDE
% Col 1: Plot acuity against age


AudStr = dp_removeoutliers(ERPaud.FitOut.Stretch);
AudAmp = dp_removeoutliers(ERPaud.FitOut.Scaling);
AudLat = dp_removeoutliers(ERPaud.FitOut.Latency);
VisStr = dp_removeoutliers(ERPvis.FitOut.Stretch);
VisLat = dp_removeoutliers(ERPvis.FitOut.Latency);
VisAmp = dp_removeoutliers(ERPvis.FitOut.Scaling);

AudStr(audind == 0) = NaN;
AudLat(audind == 0) = NaN;
AudAmp(audind == 0) = NaN;

VisStr(visind == 0) = NaN;
VisLat(visind == 0) = NaN;
VisAmp(visind == 0) = NaN;

audconvar = AudAmp;
visconvar = VisAmp;
Naud = sum(audind);
Nvis = sum(visind);

spearmans_Pcorr_acuity_RMSE_amp

%% Correcting for Baseline 

AudStr = dp_removeoutliers(ERPaud.FitOut.Stretch);
AudAmp = dp_removeoutliers(ERPaud.FitOut.Scaling);
AudLat = dp_removeoutliers(ERPaud.FitOut.Latency);
VisStr = dp_removeoutliers(ERPvis.FitOut.Stretch);
VisLat = dp_removeoutliers(ERPvis.FitOut.Latency);
VisAmp = dp_removeoutliers(ERPvis.FitOut.Scaling);

AudStr(audind == 0) = NaN;
AudLat(audind == 0) = NaN;
AudAmp(audind == 0) = NaN;
AudBL(audind == 0) = NaN;

VisStr(visind == 0) = NaN;
VisLat(visind == 0) = NaN;
VisAmp(visind == 0) = NaN;
VisBL(visind == 0) = NaN;

audconvar = AudBL;
visconvar = VisBL;

spearmans_Pcorr_acuity_RMSE_amp

%% Correcting for Visual / Auditory RMSE
% Col 1: Plot acuity against age

AudRMSE = ERPaud.FitOut.RMSE;
AudRMSE(audind == 0) = NaN;
VisRMSE = ERPvis.FitOut.RMSE;
VisRMSE(audind == 0) = NaN;

audconvar = AudRMSE;
visconvar = VisRMSE;

spearmans_Pcorr_acuity_RMSE_amp



%% Plot MSP Results Auditory

close all
figure('Color','w','position',[101 87 1383 980])
dec = 0;
climsperc = 0.995;
clims = [0 2.5];
% Load and Render the FreeSurfer surface
S = [];
S.hem = 'rh'; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 4; % 1 no inflation, 6 fully inflated
S.decimation = dec;
S = mni2fs_brain(S);

% Plot an ROI, and make it semi transparent
% ROI = load_nii('/imaging/camcan/templates/HarvardOxford-combo-maxprob-thr25-2mm.nii');
% ROI.img(ROI.img ~= 93 & ROI.img ~= 45) = 0; % Visual Cortex = 48 and 96
% S.mnivol = ROI;
% S.roicolorspec = 'b'; % color. Can also be a three-element vector
% S.roialpha = 0.5; % transparency 0-1
% S.roismoothdata = 0;
% S = mni2fs_roi(S); 

% Add overlay, theshold
S.mnivol = ['/imaging/dp01/results/proj2trans/MSPinversion/AudMean_MSP.nii'];
% S.mnivol =;
S.clims_perc = climsperc; % overlay masking below 98th percentile
% S.clims = [];
S.colormap = paruly;
S.interpmethod = 'linear';
S = mni2fs_overlay(S); 

S = [];
S.hem = 'lh'; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 4; % 1 no inflation, 6 fully inflated
S.decimation = dec;
S = mni2fs_brain(S);

% Plot an ROI, and make it semi transparent
% ROI = load_nii('/imaging/camcan/templates/HarvardOxford-combo-maxprob-thr25-2mm.nii');
% ROI.img(ROI.img ~= 93 & ROI.img ~= 45) = 0; % Visual Cortex = 48 and 96
% S.mnivol = ROI;
% S.roicolorspec = 'b'; % color. Can also be a three-element vector
% S.roialpha = 0.5; % transparency 0-1
% S.roismoothdata = 0;
% S = mni2fs_roi(S); 

% Add overlay, theshold
S.mnivol = ['/imaging/dp01/results/proj2trans/MSPinversion/AudMean_MSP.nii'];
% S.mnivol =;
S.clims_perc = climsperc; % overlay masking below 98th percentile
% S.clims = clims;
S.colormap = paruly;
S.interpmethod = 'linear';
S = mni2fs_overlay(S); 

view([ -110 10]) % change camera angle
mni2fs_lights % Dont forget to turn on the lights!

%% Plot MSP Results Visual
% Obtained from run_MSP_Inversion.m
ROI = load_nii('/imaging/camcan/templates/HarvardOxford-combo-maxprob-thr25-2mm.nii');
ROI.img(ROI.img ~= 48 & ROI.img ~= 96) = 0; % Visual Cortex = 48 and 96

close all
figure('Color','w','position',[101 87 1383 980])

clims = [0.02 0.1];
sephem = 10;
dec = 0;

% Load and Render the FreeSurfer surface
S = [];
S.hem = 'lh'; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 4; % 1 no inflation, 6 fully inflated
S.separateHem = sephem;
S.decimation = dec;
S = mni2fs_brain(S);

% % Plot an ROI, and make it semi transparent
% S.mnivol = ROI;
% S.roicolorspec = 'b'; % color. Can also be a three-element vector
% S.roialpha = 0.3; % transparency 0-1
% S.roismoothdata = 3;
% S = mni2fs_roi(S); 

% Add overlay, theshold to 98th percentile
S.mnivol = ['/imaging/dp01/results/proj2trans/MSPinversion/VisMean_MSP.nii'];
S.interpmethod = 'linear';
% S.clims_perc = 0.995;
S.clims = clims;
S.colormap = paruly;
S = mni2fs_overlay(S);

%---------------------
% % Right side
S = [];
S.hem = 'rh'; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 4; % 1 no inflation, 6 fully inflated
S.separateHem = sephem;
S.decimation = dec;
S = mni2fs_brain(S);

% Plot an ROI, and make it semi transparent
% S.mnivol = ROI;
% S.roicolorspec = 'b'; % color. Can also be a three-element vector
% S.roialpha = 0.3; % transparency 0-1
% S.roismoothdata = 3;
% S = mni2fs_roi(S); 

% Add overlay, theshold to 98th percentile
S.mnivol = ['/imaging/dp01/results/proj2trans/MSPinversion/VisMean_MSP.nii'];
% S.clims_perc = 0.99;
S.clims = clims;
S.interpmethod = 'linear';
S.colormap = paruly;
S = mni2fs_overlay(S); 

view([60 30]) % change camera angle
mni2fs_lights % Dont forget to turn on the lights!



