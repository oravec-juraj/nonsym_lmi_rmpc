% SET_AXIS
%
%   Function SET_AXIS returnes the cell-array for the function SET_LABEL.
%
%   ax = set_axis(Gmin, Gstep, Gmax, Tmin, Tstep, Tmax)
%
%   where:
%
%   ax(class:cell) - is output cell array of axis description
%   Gmin(class:double) - is input minimal number of axis grid
%   Gstep(class:double) - is input step of axis grid
%   Gmax(class:double) - is input maximal number of axis grid
%   Tmin(class:double) - is input minimal number of axis text
%   Tstep(class:double) - is input step of axis text
%   Tmax(class:double) - is input maximal number of axis text
%
%
%   ax = set_axis(1, 0.1, 2, 1.4, 0.2, 1.8)
%
%   juraj.oravec@stuba.sk
%
%   est.2014.12.20.


% Copyright is with the following author(s):
%
% (c) 2014 Juraj Oravec, Slovak University of Technology in Bratislava,
% juraj.oravec@stuba.sk

% ------------------------------------------------------------------------------
% Legal note:
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330,
% Boston, MA 02111-1307 USA
%
% ------------------------------------------------------------------------------

function ax = set_axis(Gmin, Gstep, Gmax, Tmin, Tstep, Tmax)

%% G-grid
G = [Gmin : Gstep : Gmax];

%% Initialize outputs
for k = 1 : length(G)
    ax{k} = '';
end % for k

%% Text of axis
Gidx = find(G == Tmin);
dT = Tstep / Gstep;
T = [Tmin : Tstep : Tmax];
for k = 1 : length(T)
    ax{Gidx + (k-1)*dT} = num2str(T(k));
end % for k

end % function