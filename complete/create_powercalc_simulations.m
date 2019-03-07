%% Power calculations

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
addpath /home/dp01/matlab/lib/boundedline/


Info = LoadSubIDs;

Age = Info.Age;


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


%%

[audstr, ~, vislat, ~,  age] = complete(AudStr, AudLat, VisLat, VisStr, Age);

N = length(audstr);

levels = 10:10:N;
Nlevels = length(levels);

Nrand = 10000;
R = zeros(2,Nlevels,Nrand);
P = R;
S = R;
CI = zeros(2,Nlevels,Nrand,2);

dp_matlabpool_start(60)

audvis = [audstr vislat];

%
parfor ii = 1:Nlevels
    ii
    for ri = 1:Nrand
        
        randind = randi(N, [1, levels(ii)]);
        
        for avi = 1:2
            mdl = LinearModel.fit(age(randind), audvis(randind,avi));
            R(avi,ii,ri) = mdl.Rsquared.Ordinary;
            P(avi,ii,ri) = mdl.coefTest;
            ci = mdl.coefCI;
            S(avi,ii,ri) = double(mdl.Coefficients('x1', 'Estimate'));
            SCI(avi,ii,ri,:) = ci(2,:); % 2 = confidence interval of the slope
        end        
    end
end

save('/imaging/dp01/results/av_nature/PowerCalculations.mat')

%%

close all
load('/imaging/dp01/results/av_nature/PowerCalculations.mat')
clear Sp
clear CIp

tpeak = 200;

[~, Sp(1,:,:)] = convertTemplateSlopeToPeakSlope(zeros(size(S(1,:,:))), S(1,:,:), zeros(size(S(1,:,:))), S(1,:,:), tpeak, 50);
Sp(2,:,:) = S(2,:,:);

[~, CIp(1,:,:,:)] = convertTemplateSlopeToPeakSlope(zeros(size(SCI(1,:,:,:))), SCI(1,:,:,:), zeros(size(SCI(1,:,:,:))), SCI(1,:,:,:), tpeak, 50);
CIp(2,:,:,:) = SCI(2,:,:,:);

for avi = 1:2

    Sav = squeeze(Sp(avi,:,:)); % slope
    CIav = squeeze(CIp(avi,:,:,:)); % confidence intervals for each bootstrap
    
    % get the slope value that is needed for a p<0.05..
    for pi = 1:Nlevels
        [s, ind] = sort(Sav(pi, :));
        p = P(avi, pi, ind);
        L(pi) = (s(find(p < 0.05, 1, 'first')));
    end
    
    fig(1,'position',[558 256 523 392])
    conflims = abs(quantile(Sav, [0.025 0.975], 2)-repmat(median(Sav,2), [1 2]));
    medline = median(Sav,2);
    
    % conflims = abs(quantile(CI(:,:,1), [0.025 0.975], 2)-repmat(median(CI(:,:,1),2), [1 2]));
    % medline = median(CI(:,1),2);
    
    [ax1, ax2] = boundedline(levels, medline, conflims);
    set(ax2,'FaceColor',[146, 199, 244]/255, 'FaceAlpha',0.5);
    set(ax1,'LineWidth',2)
    xlim([10 max(levels)])
    ln = line([10 N], [0 0],'color',[0 0 0],'LineStyle','--');
    
    hold on
    s = 1; % select 1 lower or 2 upper confidence limit
    conflims = quantile(CIav(:,:,s), [0.05 0.95], 2);%-repmat(median(CIav(:,:,s),2), [1 2]);
    for pi = 1:Nlevels
        pos = CIav(pi,:,s);
        power(pi) = sum(pos > 0) ./ length(pos);
    end
    
    medline = median(CIav(:,:,s),2);

    grid on
    xlabel('N Samples')
    ylabel('Slope')
    
    labs = {'Auditory', 'Visual'};
    
    STD = std(CIav(:,:,1),[],2);
    MEAN = mean(CIav(:,:,1),2);
    for pi = 1:Nlevels
        n = levels(pi);
        powerPara(pi) = sampsizepwr('t', [0 STD(pi)], MEAN(pi), [], n);
    end
    
    export_fig(sprintf('/imaging/dp01/results/av_nature/Figures/PowerCalcSlope_%s',labs{avi}),'-pdf')
    
    fig(100,'position',[558 256 523 392])
    c = {[0 126 189]/255, [244 126 0]/255};
    hold on
    plot(levels, smooth(power,5)','color', c{avi}, 'linewidth',2)
    xlabel('N samples')
    ylabel('Power = P(reject H_0 | H_1 is true)')
    xlim([10 N])
    grid on
    
    export_fig(sprintf('/imaging/dp01/results/av_nature/Figures/PowerCalc2_%s',labs{avi}),'-pdf')
end


