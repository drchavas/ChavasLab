%bindata2.m -- bin 2D data and return a statistic within each bin
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


function [ym,yb,xm] = bindata2(y,x1,x2,x1rg,x2rg)
    %function [ym,yb] = bindata2(y,x1,x2,x1rg,x2rg)
    %Computes:
    %ym(ii,jj) = mean(y(x1>=x1rg(ii) & x1 < x1rg(ii+1) & x2>=x2rg(jj) & x2 < x2rg(jj+1))
    %for every ii, jj
    %If a bin is empty it returns nan for that bin
    %using a fast algorithm which uses no looping
    %Also returns yb, the approximation of y using binning (useful for r^2
    %calculations). Example:
    %
    %x = randn(500,2);
    %y = sum(x.^2,2) + randn(500,1);
    %xrg = linspace(-3,3,10)';
    %[ym,yb] = bindata2(y,x(:,1),x(:,2),xrg,xrg);
    %subplot(1,2,1);plot3(x(:,1),x(:,2),y,'.');
    %subplot(1,2,2);h = imagesc(xrg,xrg,ym);
    %set(h,'AlphaData',~isnan(ym)); box off;
    %
    %By Patrick Mineault
    %Refs: https://xcorr.net/?p=3326
    %      http://www-pord.ucsd.edu/~matlab/bin.htm
    [~,whichedge1] = histc(x1,x1rg(:)');
    [~,whichedge2] = histc(x2,x2rg(:)');
 
%    bins1 = min(max(whichedge1,1),length(x1rg)-1);
%    bins2 = min(max(whichedge2,1),length(x2rg)-1);
 
    %bins = (bins2-1)*(length(x1rg)-1)+bins1; %original code: this
        %line puts all data outside of bin range (whichedge=0) into first/last bin!)
    bins = (bins2-1)*(length(x1rg)-1)+bins1;
    bins1=whichedge1(whichedge1>0);    %new code: error fix from commenter -- this will simply ignore data beyond bin range.
    bins2=whichedge2(whichedge2>0);    %new code: error fix from commenter -- this will simply ignore data beyond bin range.
    
    y=y(whichedge1>0 & whichedge2>0);   %only data within bin range
    x1=x1(whichedge1>0 & whichedge2>0);   %only data within bin range
    x2=x2(whichedge1>0 & whichedge2>0);   %only data within bin range
    
    
    xpos = ones(size(bins,1),1);
    ns = sparse(bins,xpos,1,(length(x1rg)-1)*(length(x2rg)-1),1);
    Nbins1_input = length(xrg1);
    Nbins2_input = length(xrg2);
    Nbins1_sparse = size(ns,1);
    Nbins2_sparse = size(ns,2);

    
    ysum = sparse(bins,xpos,y,(length(x1rg)-1)*(length(x2rg)-1),1);
    ym = full(ysum)./(full(ns));
    yb = ym(bins);
    ym = reshape(ym,length(x1rg)-1,length(x2rg)-1);
    
%     Nbins_input = length(xrg);
%     Nbins_sparse = length(ns);
%     ns(Nbins_sparse+1:Nbins_input-1) = sparse(1,1); %need to append empty sparse matrices for bins at upper end (lower end ones automatically created)
%     xsum = sparse(bins,xpos,x);
%     xsum(Nbins_sparse+1:Nbins_input-1) = sparse(1,1); %need to append empty sparse matrices for bins at upper end (lower end ones automatically created)
%     xm = full(xsum)./(full(ns));    %mean value of x within bins
%     ysum = sparse(bins,xpos,y);
%     ysum(Nbins_sparse+1:Nbins_input-1) = sparse(1,1); %need to append empty sparse matrices for bins at upper end (lower end ones automatically created)
%     ym = full(ysum)./(full(ns));    %mean value of y within bins
%     yb = ym(bins);  %bin
    
end