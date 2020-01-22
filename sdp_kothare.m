%% SDP_KOTHARE
%
%   Script SDP_KOTHARE solves the SDP designed by Kothare et al. (1996).
%
%   Requires:
%   nv(double:1)  - real scalar of number of system vertics
%   nx(double:1)  - real scalar of number of system states
%   nu(double:1)  - real scalar of number of control inputs
%   x0(double:nx) - real-valued vector of system initial conditions
%
%   juraj.oravec@stuba.sk
%
%   est.2014.07.15.
%   rev.2016.05.12.


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

%% State Update
xk = x(:,k);

%% Optimized variables
X = sdpvar(nx);
Y = sdpvar(nu,nx);
U = sdpvar(nu,nu);
g = sdpvar(1);

%% Additional Parameters
ZERO = 1e-6; % ZERO-Tolerance of strict LMIs
ZEROx = zeros(nx,nx);
ZEROux = zeros(nu,nx);
ZEROxu = zeros(nx,nu);
Ix = eye(nx);
Iu = eye(nu);

%% Objective
obj = g + trace(X);

%% Constraints
constr = [];

%% Inverted Lyapunov matrix
lmi_lyap = [X >= ZERO];

%% Robust invariant ellipsoid
lmi_rie = [ [1, xk'; xk, X] >= ZERO ];

%% Convergency
lmi_conv = [];
for v = 1 : nv

lmi_conv_item = [ [X, (A{v}*X + B{v}*Y )', (sqrt(Wx)*X)', (sqrt(Wu)*Y)';...
    A{v}*X + B{v}*Y, X, ZEROx, ZEROxu;...
    sqrt(Wx)*X, ZEROx, g*Ix, ZEROxu;...
    sqrt(Wu)*Y, ZEROux, ZEROux, g*Iu] >= ZERO ];

lmi_conv = lmi_conv + lmi_conv_item;

end % for v

%% Input constraints
lmi_u_sat = [];
if(isempty(u_sat) == 0)
    %% L2-norm
    lmi_u_sat = [ [ diag(u_sat.^2), Y;...
        Y', X] >= ZERO ];
    %% L1-norm
    lmi_u_sat = lmi_u_sat + [ [ U, Y;...
        Y', X] >= ZERO ];
    for j = 1 : nu
        lmi_u_sat = lmi_u_sat + [ U(j,j) <= u_sat(j)^2 ];
    end % for j
end % if

%% Output constraints
lmi_y_sat = [];
if(isempty(y_sat) == 0)
    for v = 1 : nv

    lmi_y_sat_item = [ [diag(y_sat{v}.^2), C{v}*(A{v}*X + B{v}*Y);...
        (C{v}*(A{v}*X + B{v}*Y))', X ] >= ZERO ];

    lmi_y_sat = lmi_y_sat + lmi_y_sat_item;

    end % for v
end % if

%% All constraints
constr = lmi_lyap + lmi_rie + lmi_conv + lmi_u_sat + lmi_y_sat;

%% Optimization
sol = solvesdp(constr, g + trace(X), yalmip_opt);

%% Results
% For controller design:
X_opt = double(X);
Y_opt = double(Y);
% For Feasibility-check:
g_opt = double(g);
U_opt = double(U);

%% Controller design
Fk = Y_opt*X_opt^(-1);