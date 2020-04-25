function h = daxes3(x, y, z, c)
% daxes3 -- Draw axes
%
%  Usage
%    daxes3(x, y, z[, c])
%
%  Inputs
%    x, y, z  origin coordinates
%    c        color
%
%  Outputs
%    h       integer, current figure handle
%
%  Description
%    Draw axes with (x,y,z) origin.
%
%  See also
%    axis, daxes
%

if (nargin == 0)
   x = 0;
   y = 0;
   z = 0;
   c = 'b';
elseif (nargin == 2)
   c = 'b';
end;

ax = axis;
hold on;
plot3([ax(1) ax(2)],[y y],[z z],c);
plot3([x x],[ax(3) ax(4)],[z z],c);
plot3([x x],[y y],[ax(5) ax(6)],c);
hold off;

h = gcf;

% Written on Fri Jan 29 15:21:59 MET 1999
% by Jean-Baptiste Caillau - ENSEEIHT-IRIT (UMR CNRS 5505)
