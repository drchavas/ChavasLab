%SS_cat_color.m
%Purpose: return Saffir-Simpson category and color
%
% Syntax: [SS_cat,SS_color,SS_str] = SS_cat_color(Vmms)
%
% Inputs:
%   Vmms [m/s] - vector of wind speeds
%
% Outputs:
%   SS_cat [] - vector of Saffir-Simpson categories (-1=TD, 0=TS, 1-5)
%   SS_color [] - cell array of corresponding RGB color vectors by category
%   SS_str [] - cell array of corresponding category strings
%
% Example: [SS_cat,SS_color,SS_str] = SS_cat_color(Vmms)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 17 Apr 2014; Last revision:

% Revision history:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------- BEGIN CODE --------------


function [SS_cat,SS_color,SS_str] = SS_cat_color(Vmms)

SS_thresh = [17.49 33.8 42.7 49.4 57.9 70.0];   %[m/s]; wind speed thresholds
SS_colors = {[.85 .85 .85] [.6 .6 .6] [.5 .75 0] [.8 .8 0] [1 .5 0] [1 0 0] [.5 0 0]};
SS_strs = {'TD' 'TS (17.5 m/s)' 'Cat1 (33.8 m/s)' 'Cat2 (42.7 m/s)' 'Cat3 (49.4 m/s)' 'Cat4 (57.9 m/s)' 'Cat5 (70.0 m/s)'};

%% Color following Saffir-Simpson scale
for ii=1:length(Vmms)
    
    Vmms_temp = Vmms(ii);
    
    if(Vmms_temp<SS_thresh(1))  %TD
        SS_color{ii} = SS_colors{1};
        SS_cat(ii) = -1;
        SS_str{ii} = SS_strs{1};
    elseif(Vmms_temp>=SS_thresh(1) && Vmms_temp<SS_thresh(2))  %TS
        SS_color{ii} = SS_colors{2};
        SS_cat(ii) = 0;
        SS_str{ii} = SS_strs{2};
    elseif(Vmms_temp>=SS_thresh(2) && Vmms_temp<SS_thresh(3))  %Cat1
        SS_color{ii} = SS_colors{3};
        SS_cat(ii) = 1;
        SS_str{ii} = SS_strs{3};
    elseif(Vmms_temp>=SS_thresh(3) && Vmms_temp<SS_thresh(4))  %Cat2
        SS_color{ii} = SS_colors{4};
        SS_cat(ii) = 2;
        SS_str{ii} = SS_strs{4};
    elseif(Vmms_temp>=SS_thresh(4) && Vmms_temp<SS_thresh(5))  %Cat3
        SS_color{ii} = SS_colors{5};
        SS_cat(ii) = 3;
        SS_str{ii} = SS_strs{5};
    elseif(Vmms_temp>=SS_thresh(5) && Vmms_temp<SS_thresh(6))  %Cat4
        SS_color{ii} = SS_colors{6};
        SS_cat(ii) = 4;
        SS_str{ii} = SS_strs{6};
    elseif(Vmms_temp>=SS_thresh(6))  %Cat5
        SS_color{ii} = SS_colors{7};
        SS_cat(ii) = 5;
        SS_str{ii} = SS_strs{7};
    end

end

%------------- END OF CODE --------------

