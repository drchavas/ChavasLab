%sinfit.m -- fit a sinusoid to one period of data
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Dan Chavas
% EAPS, Purdue University
% email: drchavas@gmail.com
% Website: -
% 2016-10-13; Last revision:

%------------- BEGIN CODE --------------

function [xdata_fit,ydata_fit,cycleinfo,gof,f,fit1] = sinfit(xdata,ydata)

%sprintf('NOTE: sinfit() assumes that input xdata is exactly one period with no repetition (i.e. [0,2*pi) )')

%% Make sure xdata and ydata are Nx1 vectors
if(size(xdata,2)>1)
    xdata = xdata';
end
if(size(ydata,2)>1)
    ydata = ydata';
end

%xdata = (1:12)';   %months
xdata_offset = xdata(1);
xdata_scalefac = (2*pi*((length(xdata)-1)/length(xdata)))/(xdata(end)-xdata_offset);
% ydata = 10^3*[1.0174 1.3080 2.1874 3.0870 3.6900 3.8400 3.8800 3.5200 ...
%    2.8000 2.1200 1.4600 1.1248]';   %aribtrary CAPE values from NARR
ydata_std = std(ydata);
ydata_mean = mean(ydata);

%%Map data to [0,2*pi) interval, standardized amplitude
xx = (xdata - xdata_offset);
xx = xx * xdata_scalefac;
%%testing
xx_test = linspace(0,2*pi,length(xdata)+1)';
xx_test = xx_test(1:end-1);  %data does not actually repeat
assert(max(abs(xx-xx_test))<10^5,'problem with mapping to [0,2pi) interval')
yy = (ydata-ydata_mean)/ydata_std;

%%TESTING
% data = sin(xx) + .1*randn(size(xx));
% yy = data;

%% Sinusoidal regression
%From https://www.mathworks.com/help/curvefit/least-squares-fitting.html
% f = fittype('a*sin(b*x)+c');
% clear fit   %fittype creates a function fit, which makes NO SENSE AT ALL
% [fit1,gof,fitinfo] = fit(xx,yy,f,'StartPoint',[1 1 0]);
% yy_fit = fit1.a*sin(fit1.b * xx)+fit1.c;

f = fittype('a*sin(x-b)+c');    %period is known; phasing is not
clear fit   %fittype creates a function fit, which makes NO SENSE AT ALL
[fit1,gof,fitinfo] = fit(xx,yy,f,'StartPoint',[1 0 0]);

xx_hires = min(xx):.01:max(xx);
yy_fit = fit1.a*sin(xx_hires-fit1.b)+fit1.c;
xdata_fit = xx_hires/xdata_scalefac + xdata_offset;   %this is the midpoint of the cycle

%% Relevant regression parameters
ydata_fit = yy_fit*ydata_std+ydata_mean;
rsquare_adj = gof.adjrsquare;
rmse = gof.rmse;

cycleinfo.ampnorm = fit1.a; %normalized amplitude of cycle
cycleinfo.amp = cycleinfo.ampnorm*ydata_std+ydata_mean;  %true amplitude of cycle
cycleinfo.x_amp0 = fit1.b/xdata_scalefac + xdata_offset;   %this is the midpoint of the cycle
cycleinfo.amp_offset = fit1.c*ydata_std+ydata_mean;   %
cycleinfo.maxval = max(ydata_fit);
cycleinfo.x_maxval = xdata_fit(ydata_fit==cycleinfo.maxval);
cycleinfo.minval = min(ydata_fit);


end
