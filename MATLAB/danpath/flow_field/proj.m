function c = proj(a,b)
%PROJ  Projection of a vector A onto vector B.
%   C = PROJ(A,B) returns the vector of the projection of A onto B.
%   A and B must be vectors of the same length.
%
% Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
%
% Copyright. June 12, 2009.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A. and R. Hernandez-Walls. (2006). proj:Projection of a 
%   vector A onto vector B. A MATLAB file. [WWW document].
%   URL http://www.mathworks.com/matlabcentral/fileexchange/24435
%
% Reference:
% Gerber, H. (1990), Elementary Linear Algebra. Brooks/Cole Pub. Co. Pacific
%     Grove, CA. 
%

if nargin < 2, 
    error('proj:Input', 'Requires exactly two input arguments.'); 
end

% Check dimensions
if any(size(a) ~= size(b)),
   error('MATLAB:proj:InputSizeMismatch', 'A and B must be same size.');
end

c = dot(a,b)/norm(b)^2*b;