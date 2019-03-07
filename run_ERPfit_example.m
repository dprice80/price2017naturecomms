close all
clear all
clc

addpath /imaging/dp01/scripts/ERPfitting

timex = (-100:500)'; % peri-stimulus time column vector

template = fspecial('gaussian',size(timex),50); % average ERP

fitcentre = 50;
tc(:,1) = interp1(timex, template, (timex - fitcentre)/1.1 + fitcentre - 10, 'spline');
tc(:,2) = interp1(timex, template, (timex - fitcentre)/1.2 + fitcentre - 20, 'spline')*2;
tc(:,3) = interp1(timex, template, (timex - fitcentre)/1.3 + fitcentre - 30, 'spline')*3;
tc(:,4) = interp1(timex, template, (timex - fitcentre)/1.4 + fitcentre - 40, 'spline')*4;

% Here the fit centre is at t = 50, and the gaussians are stretched
% outwards from that point. So the parameter FitCentre in the function is
% defined as this point. Getting the fit centre wrong throws off the
% constant delay estimations (see last cell), so this must be a sensible estimate. 

plot(timex, template, 'r', timex, tc, 'b')

%%
clc
clear fit
for fiti = 1:4
    fit(fiti) = ERPfit(tc(:,fiti), template, timex, 50, [0 500], 1, 1e-6); % Fit = [delay stretch];
    disp(fit(fiti))
end

shift = [fit.shift]
stretch = [fit.stretch]
amp = [fit.amp]


%% Incorrect fit centre 
% Setting FitCentre to 150 results in incorrect constant delay, 
% because a constant shift is required in order to compensate for incorrect
% cumulative fitcentre. i.e. imagine a time series being stretched from
% a point at -1000 seconds - expanding the time series by a factor of 1.5
% would result in the window of interest shifting +1500 seconds, so a -1500
% seconds constant shift would be needed to correct for this. The answer is
% to iteratively run the algorithm, shifting the fitcentre variable until
% the constant delay is minimised (again a simple gradient decent)
% However, this has not yet been incorporated into the code.

clear fit
for fiti = 1:4
    fit(fiti) = ERPfit(tc(:,fiti), template, timex, 150, [0 500], 1, 1e-6); % Fit = [delay stretch];
end

shift_incorrect = [fit.shift]
stretch_incorrect = [fit.stretch]
amp = [fit.amp]
