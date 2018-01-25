%num2numword.m
%Purpose: Input a number, output string word of number written out
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
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 6 Jun 2014; Last revision:

%------------- BEGIN CODE --------------

function [numword_out] = num2numword(num_in)

switch num_in
    case 1
        numword_out = 'ONE';
    case 2
        numword_out = 'TWO';
    case 3
        numword_out = 'THREE';
    case 4
        numword_out = 'FOUR';
    case 5
        numword_out = 'FIVE';
    case 6
        numword_out = 'SIX';
    case 7
        numword_out = 'SEVEN';
    case 8
        numword_out = 'EIGHT';
    case 9
        numword_out = 'NINE';
    case 10
        numword_out = 'TEN';
    case 11
        numword_out = 'ELEVEN';
    case 12
        numword_out = 'TWELVE';
    case 13
        numword_out = 'THIRTEEN';
    case 14
        numword_out = 'FOURTEEN';
    case 15
        numword_out = 'FIFTEEN';
    case 16
        numword_out = 'SIXTEEN';
    case 17
        numword_out = 'SEVENTEEN';
    case 18
        numword_out = 'EIGHTEEN';
    case 19
        numword_out = 'NINETEEN';
    case 20
        numword_out = 'TWENTY';
    case 21
        numword_out = 'TWENTY-ONE';
    case 22
        numword_out = 'TWENTY-TWO';
    case 23
        numword_out = 'TWENTY-THREE';
    case 24
        numword_out = 'TWENTY-FOUR';
    case 25
        numword_out = 'TWENTY-FIVE';
    case 26
        numword_out = 'TWENTY-SIX';
    case 27
        numword_out = 'TWENTY-SEVEN';
    case 28
        numword_out = 'TWENTY-EIGHT';
    case 29
        numword_out = 'TWENTY-NINE';
    otherwise
        assert(1==2,sprintf('String for number (%i) not coded for',num_in))
end

end


%------------- END OF CODE --------------
