function [x,y] = polar4(bs,s)
% Create a unit circle centered at (0,0) using four segments.
switch nargin
 case 0
 x = 4; % four edge segments
 return
 case 1
 A = [0,pi/2,pi,3*pi/2; % start parameter values
 pi/2,pi,3*pi/2,2*pi; % end parameter values
 1,1,1,1; % region label to left
 0,0,0,0]; % region label to right
 x = A(:,bs); % return requested columns
 return
 case 2
 x = cos(s).*(1+0.2*cos(4*s));
 y = sin(s).*(1+0.2*cos(4*s));
end
