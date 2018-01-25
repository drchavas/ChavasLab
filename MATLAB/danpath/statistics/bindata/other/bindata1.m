%bindata.m -- bin data and return a statistic within each bin
%source: https://xcorr.net/2013/12/23/fast-1d-and-2d-data-binning-in-matlab/
%with error fix from first comment (dealing with data beyond bin range)
%
% Syntax:
%
% Inputs:
%
% Outputs:
%
% Example:
%
% Other m-files required:
% Subfunctions: none
% MAT-files required: none
%
% Notes:

% Author: Dan Chavas
% EAPS, Purdue University
% email: drchavas@gmail.com
% Website: -
% 20 Jan 2017; Last revision:
% 2017-01-17: code modified to account for error in dealing with data
% beyond bin range (orig: lumped into first/last bin; now: ignored). added
% code to return mean values of x within bins too. added capacity to return
% median rather than mean.

%------------- BEGIN CODE --------------


function [ym,yb,xm] = bindata(y,x,xrg)
    %function [ym,yb] = bindata(y,x,xrg)
    %Computes ym(ii) = mean(y(x>=xrg(ii) & x < xrg(ii+1)) for every ii
    %using a fast algorithm which uses no looping
    %If a bin is empty it returns nan for that bin
    %Also returns yb, the approximation of y using binning (useful for r^2
    %calculations). Example:
    %
    %x = randn(100,1);
    %y = x.^2 + randn(100,1);
    %xrg = linspace(-3,3,10)';
    %[ym,yb] = bindata(y,x,xrg);
    %X = [xrg(1:end-1),xrg(2:end)]';
    %Y = [ym,ym]'
    %plot(x,y,'.',X(:),Y(:),'r-');
    %
    %By Patrick Mineault
    %Refs: https://xcorr.net/?p=3326
    %      http://www-pord.ucsd.edu/~matlab/bin.htm
    [~,whichedge] = histc(x,xrg(:)');   %returns bin number corresponding to x0 <= x < x1
 
    %bins = min(max(whichedge,1),length(xrg)-1);    %original code: this
    %line puts all data outside of bin range (whichedge=0) into first/last bin!)
    bins=whichedge(whichedge>0);    %new code: error fix from commenter -- this will simply ignore data beyond bin range.
    y=y(whichedge>0);   %only data within bin range
    x=x(whichedge>0);   %only data within bin range
    
    xpos = ones(size(bins,1),1);
    ns = sparse(bins,xpos,1);   %count within each bin
    Nbins_input = length(xrg);
    Nbins_sparse = length(ns);
    ns(Nbins_sparse+1:Nbins_input-1) = sparse(1,1); %need to append empty sparse matrices for bins at upper end (lower end ones automatically created)
    xsum = sparse(bins,xpos,x);
    xsum(Nbins_sparse+1:Nbins_input-1) = sparse(1,1); %need to append empty sparse matrices for bins at upper end (lower end ones automatically created)
    xm = full(xsum)./(full(ns));    %mean value of x within bins
    ysum = sparse(bins,xpos,y);
    ysum(Nbins_sparse+1:Nbins_input-1) = sparse(1,1); %need to append empty sparse matrices for bins at upper end (lower end ones automatically created)
    ym = full(ysum)./(full(ns));    %mean value of y within bins
    yb = ym(bins);  %bin
    
end