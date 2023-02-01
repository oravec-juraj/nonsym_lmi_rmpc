% MUP_VERBOSE
%
%   Function MUP_VERBOSE manages the verbose messages.
%
%   juraj.oravec@stuba.sk
%
%   est.:2013.07.17.
%

% Copyright is with the following author(s):
%
% (c) 2016 Juraj Oravec, Slovak University of Technology in Bratislava,
% juraj.oravec@stuba.sk
% (c) 2016 Monika Bakosova, Slovak University of Technology in Bratislava,
% monika.bakosova@stuba.sk
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

function mup_verbose(vbs_type,vbs,mes,varargin)

if(vbs >= vbs_type)
    cmd = ['disp(sprintf(''',mes,''','];
    for cnt = 1 : length(varargin)
        cmd = [cmd,'varargin{',num2str(cnt),'},'];
    end % for cnt
    cmd = [cmd(1:end-1),'));'];
    eval(cmd);
end % if

end % function