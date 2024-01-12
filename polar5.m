function [x,y] = polar5(bs,s)
% Create a unit circle centered at (0,0) using four segments.
switch nargin
 case 0
 x = 5; % four edge segments
 return
 case 1
 A = [0,2*pi/5,4*pi/5,6*pi/5,8*pi/5; % start parameter values
 2*pi/5,4*pi/5,6*pi/5,8*pi/5,2*pi; % end parameter values
 1,1,1,1,1; % region label to left
 0,0,0,0,0]; % region label to right
 x = A(:,bs); % return requested columns
 return
 case 2
 x = cos(s).*(1+1/6*cos(5*s));
 y = sin(s).*(1+1/6*cos(5*s));
end
