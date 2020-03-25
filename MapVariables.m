function x = MapVariables(x)
%  MapVariables interatively assign values to the elements in vector x 
%
%  Inputs:
%         Integer vector x
%        
%  Outputs:
%         A vector contains health index values
%            
%  NU_Version1.0
%  by Anqi He, Northeastern University, Boston

global xupper xdelt xlower;
% map integer variables to a discrete set

% The possible values for x(1) and x(2)
allX1 = xlower:xdelt:xupper;

% Map x(1) and x(2) from the integer values used by SA to the discrete values required.
x(1) = allX1(x(1));
x(2) = allX1(x(2));

