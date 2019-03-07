function [fit,x1] = ERPfit(x,y,timex,FitCentre,trange,showplots,difflim)
% ERP Fitting Function.
% fit = ERPfit(x,y,timex,FitCentre,trange,showplots,difflim)
% y is modified to fit x
% y is the template
% timex is the reference timecourse
% FitCentre = centre of stretch. Points will be stretched outwards from this point (time in ms)
% Xrange = range within which to calculate the fit (2 element vector) [tmin tmax]
% showplots = 1:show plots, 0 no plots
% difflim = 1e6 gives good fits in a reasonable time. However, Smaller values = faster fitting,
% with less precision, which may be necessary for some cases (i.e. fitting
% individual trials).
x = shiftdim(x);
y = shiftdim(y);

tmin = min(timex);
tmax = max(timex);
tindex = 1:length(x);
Xrange = tindex(timex >= trange(1) & timex <= trange(2));

GridShift = linspace(-100,100,20);
GridStrch = linspace(.5,2,20);
cf = zeros(length(GridShift),length(GridStrch));
if showplots == 1
    for sti = 1:length(GridStrch)
        for shi = 1:length(GridShift)
            ShiftX = GridShift(shi);
            Stretch = GridStrch(sti);
            
            ti = ((timex-FitCentre)/Stretch) + FitCentre - ShiftX;
            yi = shiftdim(interp1(timex,y,ti));
            yi(isnan(yi)) = 0;
            %  Cost
            cf(shi,sti) = corr(x(Xrange),yi(Xrange));
        end
    end
    
    [shf, stf] = find(cf == max(cf(:)));
    
    if showplots == 1 || showplots == 2
        figure(1)
        clf
        subplot(2,1,1)
        pcolor(GridStrch,GridShift,cf);shading flat
        colorbar
        xlabel('Stretch')
        ylabel('Latency')
        caxis([0 1])
    end
end

%%

% Initial Values
% currval = [
%     GridShift(shf) GridStrch(stf)
%     ];

% Random Initial Values
currval = [
    0 1
    ]; % initial values

ind  = 0;
inda = 0;
cfgrad = 10;

% currcf = cf(shf,stf); % current cost function value
% bestcf = currcf;

currcf = 0; % current cost function value
bestcf = -2;

if showplots == 1 || showplots == 2
    figure(1)
    subplot(2,1,1)
    hold on
    p1 = plot(currval(2),currval(1),'wo');
    hold off
end

eInc = exp(1i*linspace(0,2*pi,10))';
Increment = [real(eInc)*10,imag(eInc)*0.1];

Nsteps = 0;

while cfgrad > difflim
    
    Nsteps = Nsteps+1;
    
    for inci = 1:size(Increment,1)
        % pos dir
        testval = currval+Increment(inci,:);
        ShiftX = testval(1);
        Stretch = testval(2);
        ti = ((timex-FitCentre)/Stretch) + FitCentre - ShiftX;
        yi = shiftdim(interp1(timex,y,ti));
        yi(isnan(yi)) = 0;
        %  Cost
        hc(inci,1) = corr(x(Xrange),yi(Xrange));
        dir = find(hc == max(hc),1,'first');
    end
    
    updatedcf = max(hc(:));
    cfgrad = abs(currcf-updatedcf);
    currcf = updatedcf;
    
    % check if caught in loop
    if currcf > bestcf
        bestcf = currcf;
        currval = currval+Increment(dir,:);
        ind = 1;
        inda = inda + 1; % accelleration index
    else
        ind = ind + 1;
        inda = 1;
    end
    
    if ind > 1
        Increment = Increment * 0.75; % decellerate current direction
    end
    
    if inda > 1
        Increment = Increment .* 1; % accelerate current direction (1 = no accelleration)
    end
    
    if showplots == 1 || (showplots == 2 && cfgrad <= difflim)
        % plot time series
        figure(1)
        subplot(2,1,2)
        p2 = plot(timex(Xrange)',x(Xrange)'./std(x(Xrange)),'b-',timex(Xrange)',yi(Xrange)'./std(yi(Xrange)),'r--');
        xlabel('real time');ylabel('Amplitude')
        ylim([-4 4])
        set(p1,'XData',currval(2),'YData',currval(1))
        title(str('r=',num2str(bestcf,'%2.7f')))
        pause(0.05)
    end
end

if showplots == 3
        % plot time series
        fig(1,'position',[500 500 800 400])
        clf
        p2 = plot(timex(Xrange)',x(Xrange)'./std(x(Xrange)),'b-',timex(Xrange)',yi(Xrange)'./std(yi(Xrange)),'r--');
        xlabel('real time');ylabel('Amplitude')
        ylim([-4 4])
end
    
x1 = x(Xrange)';
y1 = yi(Xrange)';
bl = mean(x)-mean(y);
x1 = x1-mean(x1);
y1 = y1-mean(y1);
fit.shift = currval(1);
fit.stretch = currval(2);
fit.amp = x1/y1;
fit.error = sqrt(sum((x1-y1).^2));
fit.R2 = corr(x1',y1')^2;
fit.Nsteps = Nsteps;
fit.baseline = bl;
fit.AdjustedTemplate = yi+bl;

