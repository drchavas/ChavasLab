%vectorize_your_code.m -- example of how vectorization goes very fast
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%
% Outputs:
%
% Example:
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
% 2016-09-08; Last revision:

%------------- BEGIN CODE --------------


clear
clc
close all
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Parameters
N=10^3; %[-]; number of datapoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create two matrices of random data
data_in1 = randn(N); %NxN matrix of random normally-distributed data
data_in2 = randn(N); %NxN matrix of random normally-distributed data
nx = size(data_in1,1);
ny = size(data_in2,2);

%% Multiply each element together
%% Option 1: for loops
tic     %start a clock counter
for ii=1:nx
    
    for jj=1:ny
        
        data_out_slow(ii,jj) = data_in1(ii,jj)*data_in2(ii,jj);
        
    end
    
end
time_slow = toc     %start a clock counter

%% Option 2: vectorized matrix multiplication
tic     %start a clock counter
data_out_fast = data_in1.*data_in2;  %the dot means "apply operation to each element individually"
time_fast = toc     %start a clock counter

%% Test the code to make sure they're doing the same thing
assert(isequal(data_out_slow,data_out_fast),'Problem with code! These two matrices should be identical')

%% Write to the screen the information you care about
sprintf('For N^2 = %i datapoints, The vectorized code is %10.0f times faster than the looped code',N^2,time_slow/time_fast)


%------------- END OF CODE --------------