% SET_LABELS
%
%   Function SET_LABELS set required label notation.
%
%   set_labels(h_ax,q,data_label,data_label_str)
%
%   where
%
%   h_ax            - is input axis handle (e.g.: h_ax = gca)
%   q               - is input char of axis to be set (e.g. for x-label: q = 'x')
%   data_label      - is input column vector of label values to be set
%   data_label_str  - is optional column vector of string label notation
%
%   Example:
% 
%   x = [1:10];
%   y = x.^3;
%   plot(x,y)
%   set_labels(gca,'x',[1;2;5;6],{'t_1';'t_2';'t_5';'t_6'})
%
%   juraj.oravec@stuba.sk
%
%   est.2012.06.18.
%   rev.2013.12.18.

% Copyright is with the following author(s):
%
% (c) 2012 Juraj Oravec, Slovak University of Technology in Bratislava,
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

function set_labels(h_ax,q,data_label,data_label_str)

% Ensure column vector
%
if (size(data_label,2) > 1)
    data_label = data_label';
end

% Check DATA_LABEL_STR setting
%
if ( nargin < 4 )
   data_label_str = num2str(data_label); % In case DATA_LABEL_STR not set, use DATA_LABEL as default
else
    % Ensure the same number of characters in strings
    %
    data_label_str = force_string_length(data_label_str);
    %
    % Ensure column vector
    %
    if (size(data_label_str,1) < 2)
        data_label_str = data_label_str';
    end
end

q_lim = [q,'Lim'];
q_tick = [q,'Tick'];
q_tick_label = [q,'TickLabel'];

ax = axis;
set(h_ax,q_lim,[min(data_label),max(data_label)],q_tick,data_label,q_tick_label,data_label_str)
axis(ax)

end % function

% FORCE_STRING_LENGTH
%
%   Functino FORCE_STRING_LENGTH ensure the same number of characters in each
%   string.
%
%   str_cell_out = force_string_length(str_cell_in)
%
%   where:
%
%   str_cell_out - is output cell-array of strings
%   str_cell_in  - is input cell-array of strings
%
%   juraj.oravec@stuba.sk
%
%   est.2013.12.18.
%
function str_cell_out = force_string_length(str_cell_in)

% Align type Right, Left
%
type = 'r';
% type = 'l';

n = length(str_cell_in);

if(length(str_cell_in) == 0)
    
    % Initial value of output
    %
    str_cell_out{1} = [];
    return
    
end % if

% Find maximal length of string
%
str_max = -Inf;
for k = 1 : n
    if(length(str_cell_in{k}) > str_max)
        str_max = length(str_cell_in{k});
    end % if
end % for k

% Ensure the same length of strings
%
for k = 1 : n
    str_cell_out{k,1} = force_string_length_item(str_cell_in{k},str_max,type);
end % for k

end % function

function str_out = force_string_length_item(str_in,n,type)

% Empty symbol
%
EMPTY = ' ';

if(length(str_in) < n)
    
    str_diff = n - length(str_in);
    str_empty = EMPTY;
    
    for k = 1 : str_diff-1
        str_empty = [str_empty,EMPTY];
    end % for k

    if(isequal(type,'r'))
        str_out = [str_empty,str_in];
	elseif(isequal(type,'l'))
        str_out = [str_in,str_empty];
    end % if
    
else
    
    str_out = str_in;
    
end % if

end % function