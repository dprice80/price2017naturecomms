close all
clear all
clc

addpath /imaging/dp01/scripts/ERPfitting
addpath /imaging/dp01/scripts/av_nature/code_tests/
addpath /home/dp01/matlab/lib/paruly

timex = linspace(-1000,1500,2500)'; % peri-stimulus time column vector

noiseRange = 0;
spoint = 75; % stretch point / delay origin

showplots = false;

Nrand = 1;

constRange = linspace(-11,11,11);
cumtvRange = linspace(1/1.05, 1.05, 11);
[constGrid, cumtvGrid] = ndgrid(constRange, cumtvRange);
constGrid = reshape(constGrid, numel(constGrid), 1);
cumtvGrid = reshape(cumtvGrid, numel(cumtvGrid), 1);

timeWin = [-100 500];
peaksearchwin = [-50 50];

% CREATE TEMPLATE

for li = 1:length(noiseRange)
    noiseLevel = noiseRange(li);
    disp(str('Noise: ', noiseLevel))
    tstart = tic;
    for ii = 1:length(cumtvGrid)
        printProgress(ii,length(cumtvGrid))
        template_peaklat1 = 100;
        template_peaklat2 = 200;
        template_peaklat3 = 300;
        
        peakO = [template_peaklat1, template_peaklat2, template_peaklat3]; % original peak locations
        
        peakA = [-1, 0.5, -1];
        peakT = peakO;
        peakW = [25, 35, 100];
        template = createERPTemplate(timex, peakA, peakW, peakT);
        
        % CREATE DELAYED SIGNAL
        cumtvCentre = 75;
        cumtvVal = cumtvGrid(ii);
        constVal = constGrid(ii);
        peakA = [-1, 0.5, -1];
        peakT = (peakO-cumtvCentre)*cumtvVal+cumtvCentre+constVal;
        peakW = [25, 35, 100]; %*cumtvVal;
        tc = createERPTemplate(timex, peakA, peakW, peakT);
        
        
        % Create delayed signal from template - using interp
        timexint = (timex - cumtvCentre)*inv(cumtvVal)+cumtvCentre-constVal;
        tcint = interp1(timex, template, timexint);
        tcint(isnan(tcint)) = 0;
        tcint = zscore(tcint);
        
        %
        Fbp = [0.5 32];
        Fs = 1000;
        
        tc = tc ./ std(tc);
        
        
        for ni = 1:Nrand
            noise = randn(size(timex)) * noiseLevel;
            tcn = tc + noise;
            tcintn = tcint + noise;
            
            [~, t1] = closest(timex, timeWin(1));
            [~, t2] = closest(timex, timeWin(2));
            
            TempFit = ERPfit(tcn(t1:t2), template(t1:t2), timex(t1:t2), 75, [-100 500], showplots, 1e-10);
            peakconst.cumtv(ii, ni, li) = TempFit.stretch;
            peakconst.const(ii, ni, li) = TempFit.shift;
            
            TempFit = ERPfit(tcintn(t1:t2), template(t1:t2), timex(t1:t2), 75, [-100 500], showplots, 1e-10);
            interpconst.cumtv(ii, ni, li) = TempFit.stretch;
            interpconst.const(ii, ni, li) = TempFit.shift;
            
            tcintn_cs = cumsum(abs(tcintn(t1:t2)));
            t = timex(t1:t2);
            [val, ind] = closest(tcintn_cs, max(tcintn_cs)/2);
            interpconst.halfmax(ii, ni, li) = t(ind);
            
            tcn_cs = cumsum(abs(tcn(t1:t2)));
            t = timex(t1:t2);
            [val, ind] = closest(tcintn_cs, max(tcn_cs)/2);
            peakconst.halfmax(ii, ni, li) = t(ind);
            
            peakwin1 = peakO(1) + peaksearchwin;
            [~, t1] = closest(timex, peakwin1(1));
            [~, t2] = closest(timex, peakwin1(2));
            interpconst.peaklatency1(ii, ni, li) = timex(find(tcintn(t1:t2) == min(tcintn(t1:t2)),1,'first')+t1-1);
            
            peakwin1 = peakO(1) + peaksearchwin;
            [~, t1] = closest(timex, peakwin1(1));
            [~, t2] = closest(timex, peakwin1(2));
            peakconst.peaklatency1(ii, ni, li) = timex(find(tcn(t1:t2) == min(tcn(t1:t2)),1,'first')+t1-1);
            
            peakwin2 = peakO(2) + peaksearchwin;
            [~, t1] = closest(timex, peakwin2(1));
            [~, t2] = closest(timex, peakwin2(2));
            interpconst.peaklatency2(ii, ni, li) = timex(find(tcintn(t1:t2) == max(tcintn(t1:t2)),1,'first')+t1-1);
            
            peakwin2 = peakO(2) + peaksearchwin;
            [~, t1] = closest(timex, peakwin2(1));
            [~, t2] = closest(timex, peakwin2(2));
            peakconst.peaklatency2(ii, ni, li) = timex(find(tcn(t1:t2) == max(tcn(t1:t2)),1,'first')+t1-1);
            
            peakwin3 = peakO(3) + peaksearchwin;
            [~, t1] = closest(timex, peakwin3(1));
            [~, t2] = closest(timex, peakwin3(2));
            interpconst.peaklatency3(ii, ni, li) = timex(find(tcintn(t1:t2) == min(tcintn(t1:t2)),1,'last')+t1-1);

            peakwin3 = peakO(3) + peaksearchwin;
            [~, t1] = closest(timex, peakwin3(1));
            [~, t2] = closest(timex, peakwin3(2));
            peakconst.peaklatency3(ii, ni, li) = timex(find(tcn(t1:t2) == min(tcn(t1:t2)),1,'first')+t1-1);
            
            
            % Make 
        end
    end
    dispTimeLeft(1, 1,length(noiseRange), li, tstart)
end



%% Plot template
close all
fig(1,'color','w','Position',[560 727 414 216])
t = zscore(template);
t = t-mean(t(closestind(timex, -100):closestind(timex, 0)));
plot(timex, t, 'color', [0 160 255]/255, 'LineWidth', 2)
xlim([0 600])
xlabel('Time (ms)')
ylabel('Amplitude')
export_fig('/imaging/dp01/results/av_nature/Figures/ERP_template_example','-pdf')


%% Plot results with conf intervals
close all
clims_const = [-11 11];
clims_cumtv = [1./1.05 1.05];

P = interpconst;
for li = 1:length(noiseRange)
    
    figure('position',[31 416 1107 487],'color','white')
    
    colormap(paruly)
    subplot(2,3,1)
    imagesc(cumtvRange, constRange, reshape(mean(P.cumtv(:,:,li),2), 11,11))
    colorbar
    caxis(clims_cumtv)
    title('Varying Cumulative Delay')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    subplot(2,3,2)
    imagesc(cumtvRange, constRange, reshape(mean(P.const(:,:,li),2), 11,11))
    colorbar
    caxis(clims_const)
    title('Varying Constant Delay')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    subplot(2,3,3)
    imagesc(cumtvRange, constRange, reshape(mean(P.halfmax(:,:,li),2), 11,11))
    colorbar
    title('Halfmax')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    subplot(2,3,4)
    imagesc(cumtvRange, constRange, reshape(mean(P.peaklatency1(:,:,li),2), 11,11))
    colorbar
    title('Peak 1 Latency')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    subplot(2,3,5)
    imagesc(cumtvRange, constRange, reshape(mean(P.peaklatency2(:,:,li),2), 11,11))
    colorbar
    title('Peak 2 Latency')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    subplot(2,3,6)
    imagesc(cumtvRange, constRange, reshape(mean(P.peaklatency3(:,:,li),2), 11,11))
    colorbar
    title('Peak 3 Latency')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    export_fig(sprintf('/imaging/dp01/results/av_nature/Figures/ERP_validation_fits_%2.2d', li),'-pdf')
end

% figure('color','white')
% imagesc(cumtvRange, constRange, reshape(cumtvGrid./max(cumtvGrid),11,11) + reshape(constGrid./max(constGrid),11,11))
% title('Sum of cumtv + const')
% xlabel('Cumtv %'), ylabel('Const ms')



%% Calculate Constant and Cumulative Delay from the Peak Fitting Procedures

close all

P = interpconst; % CHOOSE PEAK LATENCY MODEL TYPE

sz = size(P.peaklatency1, 1);
Nnoise = length(noiseRange);

Dconst = zeros(sz,Nnoise);
Dcumtv = zeros(sz,Nnoise);

for ni = 1:Nnoise
    for ii = 1:sz
        x = [template_peaklat1, template_peaklat2, template_peaklat3]'-spoint;
        y = [P.peaklatency1(ii,1,ni) P.peaklatency2(ii,1,ni) P.peaklatency3(ii,1,ni)]';
        m = LinearModel.fit(x, y);
        Dconst(ii,ni) = m.Coefficients('(Intercept)',:).Estimate;
        Dcumtv(ii,ni) = m.Coefficients('x1',:).Estimate;
    end
end

p = dp_plotFitRegression(x,y);

% Plot results with conf intervals
close all

d = reshape(Dconst(:,1), 11,11);
BLfitB0 = d(6,6);
d = reshape(Dcumtv(:,1), 11,11);
BLfitB1 = d(6,6);

Dconst = Dconst-BLfitB0;
Dcumtv = Dcumtv./BLfitB1;

P = interpconst;
for li = 1:length(noiseRange)
    
    figure('position',[31 416 1107 487],'color','white')
    colormap(paruly)
    subplot(2,3,1)
    imagesc(cumtvRange, constRange, reshape(Dcumtv(:,li), 11,11))
    colorbar
    caxis(clims_cumtv)
    title('Cumulative Delay')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')
    
    colormap(paruly)
    subplot(2,3,2)
    imagesc(cumtvRange, constRange, reshape(Dconst(:,li), 11,11))
    colorbar
    caxis(clims_const)
    title('Constant Delay')
    xlabel('Cumtv. Delay (%)'), ylabel('Const. Delay (ms)')


    export_fig(sprintf('/imaging/dp01/results/av_nature/Figures/ERP_const_cumtv_derived_%2.2d', li),'-pdf')
end


%% Test noise susceptibility
close all
constRange = 10;
cumtvRange = 1.05;
noiseRange = (linspace(0,5,10));
spoint = 75; % stretch point / delay origin
showplots = false;
Nrand = 10000;

peakconst = [];
interpconst = [];
peakO1 = peakO(1);
peakO2 = peakO(2);
peakO3 = peakO(3);
[constGrid, cumtvGrid] = ndgrid(constRange, cumtvRange);
constGrid = reshape(constGrid, numel(constGrid), 1);
cumtvGrid = reshape(cumtvGrid, numel(cumtvGrid), 1);

timeWin = [-100 500];
timeWin1 = -100
timeWin2 = 500
peaksearchwin = [-50 50];


for li = 1:length(noiseRange)
    noiseLevel = noiseRange(li);
    disp(str('Noise: ', noiseLevel))
    tstart = tic;
    for ii = 1:length(cumtvGrid)
        printProgress(ii,length(cumtvGrid))
        template_peaklat1 = 100;
        template_peaklat2 = 200;
        template_peaklat3 = 300;
        
        peakO = [template_peaklat1, template_peaklat2, template_peaklat3]; % original peak locations
        
        peakA = [-1, 0.5, -1];
        peakT = peakO;
        peakW = [25, 35, 100];
        template = createERPTemplate(timex, peakA, peakW, peakT);
        
        % CREATE DELAYED SIGNAL
        cumtvCentre = 75;
        cumtvVal = cumtvGrid(ii);
        constVal = constGrid(ii);
        peakA = [-1, 0.5, -1];
        peakT = (peakO-cumtvCentre)*cumtvVal+cumtvCentre+constVal;
        peakW = [25, 35, 100]; %*cumtvVal;
        tc = createERPTemplate(timex, peakA, peakW, peakT);
        
        % Create delayed signal from template - using interp
        timexint = (timex - cumtvCentre)*inv(cumtvVal)+cumtvCentre-constVal;
        tcint = interp1(timex, template, timexint);
        tcint(isnan(tcint)) = 0;
        tcint = zscore(tcint);
        
        %
        Fbp = [0.5 32];
        Fs = 1000;
        
        tc = tc ./ std(tc);
        
        parfor ni = 1:Nrand
            noise = randn(size(timex)) * noiseLevel;
            tcn = tc + noise;
            tcintn = tcint + noise;
            
            [~, t1] = closest(timex, timeWin1);
            [~, t2] = closest(timex, timeWin2);
            
            
            TempFit = ERPfit(tcintn(t1:t2), template(t1:t2), timex(t1:t2), 75, [-100 500], showplots, 1e-10);
            interpconst_cumtv(ii, ni, li) = TempFit.stretch;
            interpconst_const(ii, ni, li) = TempFit.shift;
            
            tcintn_cs = cumsum(abs(tcintn(t1:t2)));
            t = timex(t1:t2);
            [val, ind] = closest(tcintn_cs, max(tcintn_cs)/2);
            interpconst_halfmax(ii, ni, li) = t(ind);
            
%             
            peakwin1 = peakO1 + peaksearchwin;
            [~, t1] = closest(timex, peakwin1(1));
            [~, t2] = closest(timex, peakwin1(2));
            interpconst_peaklatency1(ii, ni, li) = timex(find(tcintn(t1:t2) == min(tcintn(t1:t2)),1,'first')+t1-1);
            
            peakwin2 = peakO2 + peaksearchwin;
            [~, t1] = closest(timex, peakwin2(1));
            [~, t2] = closest(timex, peakwin2(2));
            interpconst_peaklatency2(ii, ni, li) = timex(find(tcintn(t1:t2) == max(tcintn(t1:t2)),1,'first')+t1-1);
            
            peakwin3 = peakO3 + peaksearchwin;
            [~, t1] = closest(timex, peakwin3(1));
            [~, t2] = closest(timex, peakwin3(2));
            interpconst_peaklatency3(ii, ni, li) = timex(find(tcintn(t1:t2) == min(tcintn(t1:t2)),1,'last')+t1-1);

        end
    end
    dispTimeLeft(1, 1,length(noiseRange), li, tstart)
end

interpconst.cumtv = interpconst_cumtv;
interpconst.const = interpconst_const;
interpconst.halfmax = interpconst_halfmax;
interpconst.peaklatency1 = interpconst_peaklatency1;
interpconst.peaklatency2 = interpconst_peaklatency2;
interpconst.peaklatency3 = interpconst_peaklatency3;

save('/imaging/dp01/results/av_nature/supp_fig_simulation_temp_vs_peak.mat')

%% Calculate constant and cumulative delay

load('/imaging/dp01/results/av_nature/supp_fig_simulation_temp_vs_peak.mat')

close all

P = interpconst; % CHOOSE PEAK LATENCY MODEL TYPE

sz = size(P.peaklatency1, 1);
Nnoise = length(noiseRange);

Dconst = zeros(sz,Nrand,Nnoise);
Dcumtv = zeros(sz,Nrand,Nnoise);

parfor ri = 1:Nrand
    ri
    for ni = 1:Nnoise
        for ii = 1:sz
            x = [template_peaklat1, template_peaklat2, template_peaklat3]'-spoint;
            y = [P.peaklatency1(ii,ri,ni) P.peaklatency2(ii,ri,ni) P.peaklatency3(ii,ri,ni)]';
            m = LinearModel.fit(x, y);
            Dconst(ii,ri,ni) = m.Coefficients('(Intercept)',:).Estimate;
            Dcumtv(ii,ri,ni) = m.Coefficients('x1',:).Estimate;
        end
    end
end

%% Plot Above 
% First need to make reference point
close all
constRange = 0;
cumtvRange = 1;
noiseRange = 0;
spoint = 75; % stretch point / delay origin
showplots = false;
Nrand = 1;


[constGrid, cumtvGrid] = ndgrid(constRange, cumtvRange);
constGrid = reshape(constGrid, numel(constGrid), 1);
cumtvGrid = reshape(cumtvGrid, numel(cumtvGrid), 1);

timeWin = [-100 500];
peaksearchwin = [-50 50];


for li = 1:length(noiseRange)
    noiseLevel = noiseRange(li);
    disp(str('Noise: ', noiseLevel))
    tstart = tic;
    for ii = 1:length(cumtvGrid)
        printProgress(ii,length(cumtvGrid))
        template_peaklat1 = 100;
        template_peaklat2 = 200;
        template_peaklat3 = 300;
        
        peakO = [template_peaklat1, template_peaklat2, template_peaklat3]; % original peak locations
        
        peakA = [-1, 0.5, -1];
        peakT = peakO;
        peakW = [25, 35, 100];
        template = createERPTemplate(timex, peakA, peakW, peakT);
        
        % CREATE DELAYED SIGNAL
        cumtvCentre = 75;
        cumtvVal = cumtvGrid(ii);
        constVal = constGrid(ii);
        peakA = [-1, 0.5, -1];
        peakT = (peakO-cumtvCentre)*cumtvVal+cumtvCentre+constVal;
        peakW = [25, 35, 100]; %*cumtvVal;
        tc = createERPTemplate(timex, peakA, peakW, peakT);
        
        
        % Create delayed signal from template - using interp
        timexint = (timex - cumtvCentre)*inv(cumtvVal)+cumtvCentre-constVal;
        tcint = interp1(timex, template, timexint);
        tcint(isnan(tcint)) = 0;
        tcint = zscore(tcint);
        
        %
        Fbp = [0.5 32];
        Fs = 1000;
        
        tc = tc ./ std(tc);
        
        
        for ni = 1:Nrand
            ni
            noise = randn(size(timex)) * noiseLevel;
            tcn = tc + noise;
            tcintn = tcint + noise;
            
            [~, t1] = closest(timex, timeWin(1));
            [~, t2] = closest(timex, timeWin(2));
            
            TempFit = ERPfit(tcintn(t1:t2), template(t1:t2), timex(t1:t2), 75, [-100 500], showplots, 1e-10);
            reference.cumtv(ii, ni, li) = TempFit.stretch;
            reference.const(ii, ni, li) = TempFit.shift;
            
            tcintn_cs = cumsum(abs(tcintn(t1:t2)));
            t = timex(t1:t2);
            [val, ind] = closest(tcintn_cs, max(tcintn_cs)/2);
            reference.halfmax(ii, ni, li) = t(ind);
            
            peakwin1 = peakO(1) + peaksearchwin;
            [~, t1] = closest(timex, peakwin1(1));
            [~, t2] = closest(timex, peakwin1(2));
            reference.peaklatency1(ii, ni, li) = timex(find(tcintn(t1:t2) == min(tcintn(t1:t2)),1,'first')+t1-1);
            
            peakwin2 = peakO(2) + peaksearchwin;
            [~, t1] = closest(timex, peakwin2(1));
            [~, t2] = closest(timex, peakwin2(2));
            reference.peaklatency2(ii, ni, li) = timex(find(tcintn(t1:t2) == max(tcintn(t1:t2)),1,'first')+t1-1);
            
            peakwin3 = peakO(3) + peaksearchwin;
            [~, t1] = closest(timex, peakwin3(1));
            [~, t2] = closest(timex, peakwin3(2));
            reference.peaklatency3(ii, ni, li) = timex(find(tcintn(t1:t2) == min(tcintn(t1:t2)),1,'last')+t1-1);

        end
    end
    dispTimeLeft(1, 1,length(noiseRange), li, tstart)
end


P = reference; % CHOOSE PEAK LATENCY MODEL TYPE

sz = size(P.peaklatency1, 1);
Nnoise = length(noiseRange);

Dconstref = zeros(sz,Nrand,Nnoise);
Dcumtvref = zeros(sz,Nrand,Nnoise);

for ri = 1:Nrand
    ri
    for ni = 1:Nnoise
        for ii = 1:sz
            x = [template_peaklat1, template_peaklat2, template_peaklat3]'-spoint;
            y = [P.peaklatency1(ii,ri,ni) P.peaklatency2(ii,ri,ni) P.peaklatency3(ii,ri,ni)]';
            m = LinearModel.fit(x, y);
            Dconstref = m.Coefficients('(Intercept)',:).Estimate;
            Dcumtvref = m.Coefficients('x1',:).Estimate;
        end
    end
end


% The above portion of this cell is only intended to generate a reference
% point for calculating constant and cumulative delay
close all
constRange = 10;
cumtvRange = 1.05;
noiseRange = (linspace(0,5,10));
spoint = 75; % stretch point / delay origin
showplots = false;
Nrand = 100;


[constGrid, cumtvGrid] = ndgrid(constRange, cumtvRange);
constGrid = reshape(constGrid, numel(constGrid), 1);
cumtvGrid = reshape(cumtvGrid, numel(cumtvGrid), 1);

timeWin = [-100 500];
peaksearchwin = [-50 50];


close all
Dconst = squeeze(Dconst);
Dcumtv = squeeze(Dcumtv);
meanDcon = mean(Dconst-Dconstref, 1);
stdDcon = std(Dconst-Dconstref, 0, 1);
meanDcumtv = mean(Dcumtv./Dcumtvref, 1);
stdDcumtv = std(Dcumtv./Dcumtvref, 0, 1);


%% Plot delay estimations for both const and cumtv
addpath /home/dp01/matlab/lib/boundedline/
col1 = [0 126 130]/255;
col2 = [80 189 0]/255;

% Constant Delay
Tconst = squeeze(interpconst.const);
meanTcon = mean(Tconst, 1);
stdTcon = std(Tconst, 0, 1);
close all
f1 = fig(1,'color','w','position',[560 671 420 500])
subplot(2,1,1)
[f,p] = boundedline(noiseRange, meanDcon, stdDcon','alpha');
set(f,'Color',col1,'linewidth',2)
set(p,'FaceColor',col1)
hold on
[f,p] = boundedline(noiseRange, meanTcon, stdTcon','alpha');
set(f,'Color',col2,'linewidth',2)
set(p,'FaceColor',col2)
ylabel('Estimated Delay Error (ms)')
xlabel('Simulated Noise Level (SNR{^{-1}})')
line([noiseRange(1) noiseRange(end)],[10 10], 'color','k','linestyle','--')
axis tight
xlim([noiseRange(1) noiseRange(end)])

% Cumulative Delay 
Tcumtv = squeeze(interpconst.cumtv);
meanTcumtv = mean(Tcumtv, 1);
stdTcumtv = std(Tcumtv, 0, 1);
subplot(2,1,2)
[f,p] = boundedline(noiseRange, meanDcumtv, stdDcumtv','alpha');
set(f,'Color',col1,'linewidth',2)
set(p,'FaceColor',col1)
hold on
[f,p] = boundedline(noiseRange, meanTcumtv, stdTcumtv','alpha');
set(f,'Color',col2,'linewidth',2)
set(p,'FaceColor',col2)
ylabel('Estimated Delay Error (ms)')
xlabel('Simulated Noise Level (SNR{^{-1}})')
line([noiseRange(1) noiseRange(end)],[1.05 1.05], 'color','k','linestyle','--')
% axis tight
xlim([noiseRange(1) noiseRange(end)])
ylim([0.8 1.1])

export_fig('/imaging/dp01/results/av_nature/Figures/ERP_template_vs_peak_error','-bmp','-m1')
