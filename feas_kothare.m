% FEAS_KOTHARE
%
%   FEAS_KOTHARE evaluates feasibility check of controller design
%   optimization problem.
%
%   juraj.oravec@stuba.sk
%
%   est.:2013.07.17.
%   rev.:2016.05.12.

% Copyright is with the following author(s):
%
% (c) 2016 Juraj Oravec, Slovak University of Technology in Bratislava,
% juraj.oravec@stuba.sk
% (c) 2016 Michal Kvasnica, Slovak University of Technology in Bratislava,
% michal.kvasnica@stuba.sk
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

%% Output initialization
chk = 0;
chk_eig = -Inf;

%% Timer initilaization
t0 = toc;

%% Counter initialization
cnt = 0;

%% Feasibility-Check tolerance
ZERO_tol = -1e-6;

%% Inverted Lyapunov matrix
entity = min(eig(X_opt));
mup_verbose(2,vbs,' FEAS_CHK: Inverse Lyapunov matrix MIN_EIG = %x.',entity)
if(entity < ZERO_tol);
    cnt = cnt + 1;
    chk(cnt) = 1;
    chk_eig(cnt) = entity;
    mup_verbose(2,vbs,' FEAS_CHK: Inverse Lyapunov matrix failed!')
end % if

%% Robust invariant ellipsoid
entity = min(eig([1, xk'; xk, X_opt]));
mup_verbose(2,vbs,' FEAS_CHK: Invariant ellipsoid matrix MIN_EIG = %x.',entity)
if(entity < ZERO_tol);
   cnt = cnt + 1;
   chk(cnt) = 1;
   chk_eig(cnt) = entity;
   mup_verbose(2,vbs,' FEAS_CHK: Invariant ellipsoid failed!')
end % if

%% Convergency
for v = 1 : nv
entity = min(eig([...
    X_opt,...
    X_opt*A{v}' + Y_opt'*B{v}',...
    X_opt*Wx^(1/2),...
    Y_opt'*Wu^(1/2);...
    A{v}*X_opt + B{v}*Y_opt,...
    X_opt,...
    ZEROx,...
    ZEROxu;...
    Wx^(1/2)*X_opt,...
    ZEROx,...
    g_opt*Ix,...
    ZEROxu;...
    Wu^(1/2)*Y_opt,...
    ZEROux,...
    ZEROux,...
    g_opt*Iu,...
   ]));
mup_verbose(2,vbs,' FEAS_CHK: VTX:%d: Stability condition matrix MIN_EIG = %x.',v,entity)
if(entity < ZERO_tol);
   cnt = cnt + 1;
   chk(cnt) = 1;
   chk_eig(cnt) = entity;
   mup_verbose(2,vbs,' FEAS_CHK: VTX:%d: Stability condition failed!',v)
end % if 
end % for v : nv

%% Input constraints
if(isempty(u_sat) == 0)
%% Part I
entity = min(eig([diag(u_sat.^2), Y_opt; Y_opt', X_opt ]));
mup_verbose(2,vbs,' FEAS_CHK: Input contraints matrix L2 MIN_EIG = %x.',entity)
if(entity < ZERO_tol);
    cnt = cnt + 1;
    chk(cnt) = 1;
    chk_eig(cnt) = entity;
    mup_verbose(2,vbs,' FEAS_CHK: Input contraints matrix failed!')
end
%% Part II
entity = min(eig([ U_opt, Y_opt; Y_opt', X_opt ]));
mup_verbose(2,vbs,' FEAS_CHK: Input contraints matrix L1 MIN_EIG = %x.',entity)
if(entity < ZERO_tol);
    cnt = cnt + 1;
    chk(cnt) = 1;
    chk_eig(cnt) = entity;
    mup_verbose(2,vbs,' FEAS_CHK: Input contraints matrix failed!')
end % if
%% Part III
entity = min([diag(U_opt) - u_sat.^2]);
mup_verbose(2,vbs,' FEAS_CHK: Input contraints matrix MIN_EIG = %x.',entity)
if(entity > -ZERO_tol);
    cnt = cnt + 1;
    chk(cnt) = 1;
    chk_eig(cnt) = entity;
    mup_verbose(2,vbs,' FEAS_CHK: Input contraints matrix failed!')
end % if
end % if

%% State constraints
if(isempty(y_sat) == 0)
for v = 1 : nv
entity = min(eig([
    X_opt,...
    (A{v}*X_opt + B{v}*Y_opt)'*C{v}';...
    C{v}*(A{v}*X_opt + B{v}*Y_opt),...
    diag(y_sat{v}.^2),...
    ]));
mup_verbose(2,vbs,' FEAS_CHK: VTX:%d: State contraints matrix MIN_EIG = %x.',v,entity)
if(entity < ZERO_tol);
   cnt = cnt + 1;
   chk(cnt) = 1;
   chk_eig(cnt) = entity;
   mup_verbose(2,vbs,' FEAS_CHK: VTX:%d: State contraints failed!',v)
end % if 
end % for v
end % if

%% Feasibility Check - Summary
t_load2 = toc;
t_load = t_load2 - t0;
if(sum(chk) > 0)
    chk_feas = 0;
    feas_info = 'Failed!';
    mup_verbose(2,vbs,'  %d. step: FEAS_CHK: Failed! (%d infeasible constraint(s) found) (%.2f)',k,sum(chk),t_load)
    mup_verbose(2,vbs,'  %d. step: FEAS_CHK: Maximal eigenvalue: %x ',k,max(chk_eig))
else
    chk_feas = 1;
    feas_info = 'Valid';
    mup_verbose(2,vbs,'  %d. step: FEAS_CHK: Valid. (%.2f)',k,t_load)
end % if
