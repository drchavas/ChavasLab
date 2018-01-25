%VI_calculate.m
%Purpose: Calculate VI from components
%
% Syntax:  [VI] = VI_calculate(Vp,shr,entropy_var)
%
% Inputs:
%   Vp [m/s] - local potential intensity
%   shr [m/s] - bulk environmental wind shear
%   entropy_var [-] - structure; field 'chim' is non-dimensional entropy
%       deficit; o.w. fields 'sdef' entropy deficit and 'asdeq' air-sea
%       disequilibrium used to calculate chim
%
% Outputs:
%   VI [-] - ventilation index
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 2015-10-16; Last revision:

%------------- BEGIN CODE --------------


function [VI] = VI_calculate(Vp,shr,entropy_var)

if(~isfield(entropy_var,'chim'))    %chim not given, must be calculated
    
    if(~isfield(entropy_var,'sdef') || ~isfield(entropy_var,'asdeq'))
        error('cannot calculate chim; requires input fields sdef and asdeq in entropy structure variable')
    end

    %%Calculate chim = entropy deficit / air-sea disequilibrium
    
    %do not allow supersaturation (negative numerator); also limit smallness of denominator
    chim = max(0,sstarminnerbar-smouterbar)./max(1,asdeq);
    %eliminate points where net enthalpy flux into ocean
    chim(asdeq<0) = NaN;

else
    
    chim = entropy_var.chim;
    
end

%%Calculate VENTILATION INDEX
VI = shr.*chim./max(Vp,1);
VI(isnan(Vp)) = NaN;
VI(Vp==0) = NaN;

end