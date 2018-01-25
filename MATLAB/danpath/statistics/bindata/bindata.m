%bindata.m -- bin 2D data and return a statistic within each bin
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
% 2017-01-17: based in part on https://xcorr.net/2013/12/23/fast-1d-and-2d-data-binning-in-matlab/.

%------------- BEGIN CODE --------------


function [ym,x1b,x2b,x1m,x2m,nn] = bindata(y,x1,x2,x1rg,x2rg,statistic)
    %bindata.m -- bin 2D data and return a statistic within each bin
    %See bindata_example.m for usage
    %Inputs:
    %y: input array of data
    %x1: input array of x1 values for data
    %x2: input array of x2 values for data
    %x1rg: input vector for binning x1 values (must be monotonically increasing)
    %x2rg: input vector for binning x2 values (must be monotonically increasing)
    %statistic: 'mean' or 'median' (optional; default = 'mean')
    %
    %Outputs:
    %ym: array of summary statistic within bins
    %x1b: array of x1 binned values (bin center)
    %x2b: array of x2 binned values (bin center)
    %x1m: array of x1 bin-mean values (mean of x1 values in bin)
    %x2m: array of x2 bin-mean values (mean of x2 values in bin)
    %nn: array of counts within bins (i.e. a histogram)
    
    if(nargin<6)
        statistic = 'mean';
    end
    
    sprintf('statistic returned is %s',statistic)

    [~,whichedge1] = histc(x1,x1rg(:)');
    [~,whichedge2] = histc(x2,x2rg(:)');
        %note: N(k) will count the value X(i) if EDGES(k) <= X(i) < EDGES(k+1)
        %note: the last bin will count any values of X that match EDGES(end)
        %note: anything outside the range is assigned whichedge = 0
    %Remove the final bin (values of X that match EDGES(end))
    whichedge1(whichedge1==length(x1rg(:))) = 0;
    whichedge2(whichedge2==length(x2rg(:))) = 0;
     
    %% Set data outside of range to NaN (just for precaution)
    x1(whichedge1==0) = NaN;   %only data within bin range
    x2(whichedge2==0) = NaN;   %only data within bin range
    y(whichedge1==0 | whichedge2==0) = NaN;   %only data within bin range
    
    
    %% Initialize output matrices
    nbins1 = length(x1rg) - 1;
    nbins2 = length(x2rg) - 1;
    x1b = NaN(nbins1,nbins2);
    x2b = NaN(nbins1,nbins2);
    x1m = NaN(nbins1,nbins2);
    x2m = NaN(nbins1,nbins2);
    ym = NaN(nbins1,nbins2);
    nn = zeros(nbins1,nbins2);
    
    %% Loop through bins
    iedges1 = unique(whichedge1(whichedge1>0));
    iedges2 = unique(whichedge2(whichedge2>0));
    assert(isequal(size(whichedge1),size(y)),'problem with matrices')
    assert(isequal(size(whichedge2),size(y)),'problem with matrices')
    assert(isequal(size(x1),size(y)),'problem with matrices')
    assert(isequal(size(x2),size(y)),'problem with matrices')
    for ii=1:nbins1
        
        sprintf('x1 bin: [%3.2f,%3.2f)',x1rg(ii),x1rg(ii+1))
        
        for jj=1:nbins2
            
            sprintf('x2 bin: [%3.2f,%3.2f)',x2rg(jj),x2rg(jj+1))
            
            %% Pure bin center (no NaNs for plotting)
            x1b(ii,jj) = mean([x1rg(ii) x1rg(ii+1)]); %bin center
            x2b(ii,jj) = mean([x2rg(jj) x2rg(jj+1)]); %bin center
            
            indices_good = find(whichedge1==ii & whichedge2==jj);
            
            if(~isempty(indices_good))
                
                nn_temp = length(indices_good);
                sprintf('N=%i',nn_temp)
                nn(ii,jj) = nn_temp;
                
                x1m(ii,jj) = nanmean(x1(indices_good));
                x2m(ii,jj) = nanmean(x2(indices_good));
                switch statistic
                    case 'mean'
                        ym(ii,jj) = nanmean(y(indices_good));
                    case 'median'
                        ym(ii,jj) = nanmedian(y(indices_good));
                    otherwise
                        error('input statistic must be mean or median')
                end
            end
            
        end        
    end
    
    sprintf('statistic returned is %s',statistic)
    sprintf('x1 bin range: [%3.2f,%3.2f); x2 bin range: [%3.2f,%3.2f);',x1rg(1),x1rg(end),x2rg(1),x2rg(end))
    sprintf('Total number of datapoints within input bin ranges: N=%i',nansum(nn(:)))
    
end