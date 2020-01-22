%% RMPC_NONSYM
%
%   Function RMPC_NONSYM implements the LMI-based quasi-non-symmetric
%   constrained robust MPC.
%
%   data = rmpc_nonsym(model,design,setup)
%
%   where:
%
%   data(class:struct)   - is output structure with robust MPC results,
%   model(class:struct)  - is input structure with uncertain system,
%   design(class:struct) - is input structure with robust MPC design parameters,
%   setup(class:struct)  - is input structure with robust MPC setup.
%
%   To run DEMO type:
%   >>rmpc_nonsym demo
%
%   juraj.oravec@stuba.sk
%
%   est.:2016.05.12.

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

function [data] = rmpc_nonsym(model,design,setup)

%% ------------------------------------------------ %
%  PARSE INPUTS
%  ------------------------------------------------ %

%% DEMO
if(nargin == 1)
    if(isequal(model,'demo'))
        disp(sprintf(' RMPC_NONSYM: Run DEMO.'))
        %% Load DEMO model
        model = load_model_demo;
        %% Load DEMO design
        design = load_design_demo;
    end % if
end % if

%% SETUP
if(nargin < 3)
    % Set default values
    setup.vbs = 1; % Verbose-mode
    setup.flag_fig = 1; % Show figures
    setup.yalmip_verbose = 0; % YALMIP verbose-mode
    setup.sdp_solver = 'mosek'; % SDP solver
    % setup.sdp_solver = 'sedumi';
    % setup.sdp_solver = 'sdpt3';
    % setup.sdp_solver = 'sdpa';
end % if

%% Extract structure MODEL
A = model.A;
B = model.B;
C = model.C;

%% Extract structure DESIGN
simstep = design.simstep;
%
x0 = design.x0;
%
Wu = design.Wu;
Wx = design.Wx;
%
u_satL = design.u_satL;
u_satU = design.u_satU;
%
y_satL = design.y_satL;
y_satU = design.y_satU;

%% Extract structure DESIGN
vbs = setup.vbs;
flag_fig = setup.flag_fig;
yalmip_verbose = setup.yalmip_verbose;
sdp_solver = setup.sdp_solver;

%% YALMIP Setup
yalmip_opt = sdpsettings('solver',sdp_solver,'verbose',yalmip_verbose);

%% ------------------------------------------------ %
%  DATA PROCESSING
%  ------------------------------------------------ %

%% Problem size
% Number of states
nx = size(A{1},1);
% Number of inputs
nu = size(B{1},2);
% Number of outputs
ny = size(C{1},1);
% Number of vertex systems
nv = size(A,1);

%% Support Controller Design
% Tight input constraints
u_sat = min(u_satL,u_satU);
% Tight outputs constraints
for vtx = 1 : nv
    y_sat{vtx,1} = min(y_satL,y_satU);
end % for vtx
% SDP initialization
x(:,1) = x0; 
k = 1;
% solve SDP
sdp_support
% Gain matrix of support controller 
F0 = Fk;
% Feasibility check
tic
feas_support
% Verbose
mup_verbose(1,vbs,' Support Controller design:   FEAS_CHK: %s, YALMIP: %s.',feas_info,sol.info)

%% Closed-loop control
% For each vertex
for vtx = 1 : nv

%% Initialization
x = [x0, x0];
u = zeros(nu,1);
F{1} = zeros(nu,nx);

% For each control step
for k = 2 : simstep

%% Support control input
ut = F0*x(:,k-1);

%% Support System Output
if(isempty(y_satL) == 0)
    for vt = 1 : nv
        yt{vt,1} = C{vt}*(A{vt} + B{vt}*F0)*x(:,k-1);
    end % for vt
end % if

%% Update the time-varying vector of control input
% Check presence of non-symmetric input constraints
if(sum(u_satL ~= u_satU) > 0)
    %% Non-symmetric constraints found
    % Previous control input
    u_pre = u(:,k-1);
    % Update "u_sat"
    u_sat = update_u_sat(ut,u_pre,u_satL,u_satU);
elseif(sum(u_satL == u_satU) == length(u_satU))
    %% Symmetric constraints found
    u_sat = u_satL;
else
    %% No action
end % if

%% Update the Time-varying Vector of System Outputs
% Check presence of non-symmetric output constraints
if(sum(y_satL ~= y_satU) > 0)
    %% Non-symmetric constraints found
    % Previous system output
    y_pre = C{vtx}*x(:,k-1);
    % Update "y_sat"
    y_sat = update_y_sat(yt,y_pre,y_satL,y_satU);
elseif(sum(y_satL == y_satU) == length(y_satU))
    %% Symmetric constraints found
    for v = 1 : nv
        y_sat{v,1} = y_satL;
    end % for vtx
else
    %% No Action
end % if

%% Store time-varying constraints
U_sat{vtx,1}{k,1} = u_sat;
Y_sat{vtx,1}{k,1} = y_sat;

%% Solve SDP designed by Kothare et al., Automatica, 1996. (YALMIP required)
sdp_kothare
% Store optimization outputs
chk_info{k} = sol.info;
chk_state(k) = sol.problem;
% Check feasibility
feas_kothare
% Store "chk_feas"
chk_feas_all{vtx,1}(k,1) = chk_feas;

%% Handle results of feasibility check
if(chk_state(k) == 1)
    %% Infeasible problem found
    % Implement support control action
    u(:,k) = ut;
    F{k} = F0;
    % Verbose
    mup_verbose(1,vbs,' Support Controller implemented!')
    % Store "chk_support"
    chk_support{vtx,1}(k,1) = 1;
else
    %% SDP is feasible 
    % Implmennt designed controller
    F{k} = Fk;
    % Evaluate control law
    uk = Fk*xk;
    %% Check constraints on control inputs
    if(uk < -u_satL)|(uk > u_satU)
        %% Found violation of control inputs constraints
        % Implement support control action
        u(:,k) = ut;
        % Verbose
        mup_verbose(1,vbs,' Support Controller implemented!')
        % Store "chk_support"
        chk_support{vtx,1}(k,1) = 1;
    else
        %% Control input constraints hold
        % Implement support control action
        u(:,k) = uk;
        % Store "chk_support"
        chk_support{vtx,1}(k,1) = 0;
    end % if
end % if

%% Closed-loop control of "vtx"-th vertex of uncertain system
x(:,k+1) = A{vtx}*x(:,k) + B{vtx}*u(:,k);
% Verbose
mup_verbose(1,vbs,' Vertex: %d of %d, step %d of %d, FEAS_CHK: %s, YALMIP: %s.',vtx,nv,k,simstep,feas_info,chk_info{k})

end % for k

%% Quality criterion J
[J,Jt] = get_cost(x,u,Wx,Wu);

% Store vertex-dependent RMPC results into structure "data"
data.x{vtx,1} = x;
data.u{vtx,1} = u;
data.F{vtx,1} = F;
data.Jk{vtx,1} = Jt;
data.J{vtx,1} = J;

end % vtx

%% Store RMPC results into structure "data"
data.y_satL = y_satL;
data.y_satU = y_satU;
data.u_satL = u_satL;
data.u_satU = u_satU;
data.Wx = Wx;
data.Wu = Wu;
data.chk_feas = chk_feas_all;
data.chk_support = chk_support;
data.U_sat = U_sat;
data.Y_sat = Y_sat;
% Export RMPC results into MAT-file
chk = store_data('results_',data);

%% Show Results
if(flag_fig == 1)
% Initialize figure
figure(1), hold on
temp_pos = get(gcf,'position');
set(gcf,'position',[temp_pos(1:3), temp_pos(4)*2])
set(gcf,'color','white')
%% System state x_1
subplot(3,1,1), hold on, box on
plot([0,simstep],[y_satL(1),y_satL(1)],'d--',[0,simstep],[0,0],'--',[0,simstep],[-y_satL(1),-y_satL(1)],'d--')
plot([0,simstep],[y_satU(1),y_satU(1)],'d--',[0,simstep],[0,0],'--',[0,simstep],[-y_satU(1),-y_satU(1)],'d--')
axis([0, simstep, -max(y_satU(1),y_satL(1))*1.05, max(y_satU(1),y_satL(1))])
%% System state x_2
subplot(3,1,2), hold on, box on
plot([0,simstep],[y_satL(2),y_satL(2)],'s--',[0,simstep],[0,0],'--',[0,simstep],[-y_satL(2),-y_satL(2)],'d--')
plot([0,simstep],[y_satU(2),y_satU(2)],'s--',[0,simstep],[0,0],'--',[0,simstep],[-y_satU(2),-y_satU(2)],'d--')
axis([0, simstep, -max(y_satU(2),y_satL(2))*1.05, max(y_satU(2),y_satL(2))])
%% Control inputs u
subplot(3,1,3), hold on, box on
plot([0,simstep],[u_satL,u_satL],'--',[0,simstep],[0,0],'--',[0,simstep],[-u_satL,-u_satL],'--')
plot([0,simstep],[u_satU,u_satU],'--',[0,simstep],[0,0],'--',[0,simstep],[-u_satU,-u_satU],'--')
axis([0, simstep, -max(u_satL,u_satU)*1.05, max(u_satL,u_satU)*1.05])
%% Plot control trajectories
for vtx = 1 : nv
x = data.x{vtx};
u = data.u{vtx};
F = data.F{vtx};
%% System state x_1
subplot(3,1,1)
plot(0:simstep,x(1,:),'-')
xlabel('t')
ylabel('x_1')
title('Robust Constrained MPC using LMIs')
%% System state x_2
subplot(3,1,2)
hold on
plot(0:simstep,x(2,:),'-')
xlabel('t')
ylabel('x_2')
%% Control inputs u
subplot(3,1,3)
hold on
plot(1:k,u(1,:),'-')
xlabel('t')
ylabel('u')
end % for v

end % if flag_fig

end % function

function u_sat = update_u_sat(ut,u_pre,u_min,u_max)
%% UPDATE_U_SAT
%
%   Function UPDATE_U_SAT returnes the time-varying vector of control input
%   saturation. This function serves to design robust MPC with
%   quasi-non-symmetric constraints.
%
%   u_sat = update_u_sat(ut,u_pre,u_min,u_max)
%
%   where:
%
%   u_sat(double:nu) - is output real-valued vector of time-varying input constraints
%   ut(double:nu)    - is input real-valued vector of support control input
%   u_pre(double:nu) - is input real-valued vector of last control input
%   u_min(double:nu) - is input real-valued vector of lower-bound control input constraints
%   u_max(double:nu) - is input real-valued vector of upper-bound control input constraints
%   nu               - is number of control inputs
%
%   Type:
%
%   u_sat = update_u_sat([-5;0;2],[-6;-1;3],[-7;-8;-9],[9;8;7])
%
%   See also: UPDATE_Y_SAT.
%
%   juraj.oravec@stuba.sk
%
%   est.:2014.07.15
%   rev.:2016.05.12

%% Signum Function Evaluation
s = sign(ut);

%% Fix Crossing Zero
s0_idx = find(s == 0);
for k = 1 : length(s0_idx)
    s(s0_idx(k)) = -sign(u_pre(s0_idx(k)));
end % for k

%% Switching-Indicator Function
for k = 1 : length(ut)
    Iu(k,1) = max(s(k),0);
end % for k

%% Separation Matrices
U_max = diag(Iu);
U_min = diag(1-Iu);

%% Time-varying Vector of Saturation
u_sat = U_max*u_max + U_min*u_min;

end % function

function y_sat = update_y_sat(yt,y_pre,y_min,y_max)
%% UPDATE_Y_SAT
%
%   Function UPDATE_Y_SAT returnes the time-varying vector of system outputs
%   saturation. This function serves to design robust MPC with
%   quasi-non-symmetric constraints.
%
%   y_sat = update_y_sat(yt,y_pre,y_min,y_max)
%
%   where:
%
%   y_sat(cell:nv)   - is output real-valued vector of time-varying system outputs
%   yt(cell:nv)      - is input real-valued vector of support system outputs
%   y_pre(double:ny) - is input real-valued vector of last system outputs
%   y_min(double:ny) - is input real-valued vector of lower-bound system outputs constraints
%   y_max(double:ny) - is input real-valued vector of upper-bound system outputs constraints
%   nv               - is number of system vertices
%   ny               - is number of system outputs
%
%   See also: UPDATE_U_SAT.
%
%   juraj.oravec@stuba.sk
%
%   est.:2014.07.15
%   rev.:2016.05.12

%% Number of System Vertices
nv = length(yt);

%% Evaluate For Each Vertex System
for v = 1 : nv

%% Signum Function Evaluation
s{v,1} = sign(yt{v});

%% Fix Crossing Zero
s0_idx{v,1} = find(s{v} == 0);
for k = 1 : length(s0_idx{v})
    s{v}(s0_idx{v}(k)) = -sign(y_pre(s0_idx{v}(k)));
end % for k

%% Switching-Indicator Function
for k = 1 : length(yt{v})
    Iy{v}(k,1) = max(s{v}(k),0);
end % for k

%% Separation Matrices
Y_max{v} = diag(Iy{v});
Y_min{v} = diag(1-Iy{v});

%% Time-varying Vector of Saturation
y_sat{v,1} = Y_max{v}*y_max + Y_min{v}*y_min;

end % for v

end % function

function chk = store_data(filename,data)
%% STORE_DATA
%
%   Function STORE_DATA stores the data into external file.
%
%   chk = store_data(filename,data)
%
%   juraj.oravec@stuba.sk
%
%   est.:2016.05.12.

filename = [filename,datestr(now,30)];
save(filename,'data','-v6')

chk = (exist([filename,'.mat']) == 2);

end % function

function mup_verbose(vbs_type,vbs,mes,varargin)
%% MUP_VERBOSE
%
%   Function MUP_VERBOSE manages the verbose messages.
%
%   juraj.oravec@stuba.sk
%
%   est.:2013.07.17.
%   rev.:2016.05.12.

if(vbs >= vbs_type)
    cmd = ['disp(sprintf(''',mes,''','];
    for cnt = 1 : length(varargin)
        cmd = [cmd,'varargin{',num2str(cnt),'},'];
    end % for cnt
    cmd = [cmd(1:end-1),'));'];
    eval(cmd);
end % if

end % function

function [J,Jt] = get_cost(x,u,Wx,Wu)
%% GET_COST
%
%   Function GET_COST computes the quardratic quality criterion value.
%
%   [J,Jt] = get_cost(x,u,Wx,Wu)
%
%   juraj.oravec@stuba.sk
%
%   est.:2016.05.12.

tF = size(u,2);

J = 0;

for t = tF : -1 : 1

xt = x(:,t);
ut = u(:,t);

%% Stage Cost
Jt(t,1) = xt'*Wx*xt + ut'*Wu*ut;

end % for t

%% Total cost
J = sum(Jt);

end % function

function [design] = load_design_demo()
%% Load DEMO design

%% Number of simulation steps
design.simstep = 55;

%% Input constraints
% Upper bound
u_satU = 0.3; % u_satL = 1.3; % Non-symmetric ZZZ
% u_satU = 0.3; % % Symmetric ZZZ
% Lower bound (absolute value)
u_satL = 1.3; % Non-symmetric ZZZ
% u_satL = 0.3; % Symmetric ZZZ

%% Ouptut constraints (not handled)
% Upper bound
y_satU = [50; 5]; % Non-symmetric ZZZ
% y_satU = [45; 5]; % Symmetric ZZZ

% Lower bound (absolute value)
% y_satL = [45; 6]; % Non-symmetric ZZZ
y_satL = [45; 5]; % Symmetric ZZZ

%% Initial conditions
x0 = [-1; 4];
u0 = [0];

%% Parameters of quality criterion weighting matrices
Wx = diag([1; 1]);
Wu = diag([1]);

%% Generate function outputs
design.x0 = x0;
%
design.Wu = Wu;
design.Wx = Wx;
%
design.u_satL = u_satL;
design.u_satU = u_satU;
%
design.y_satL = y_satL;
design.y_satU = y_satU;

end % function

function [model] = load_model_demo()
%% Load DEMO system

%% Uncertain state matrix A
A{1,1} = [0.94, 1; 0, 1];
A{2,1} = [1.00, 1; 0, 1];
A{3,1} = [0.94, 0.94; 0, 1];
A{4,1} = [1.00, 0.94; 0, 1];
A{5,1} = [0.94, 1; 0, 1];
A{6,1} = [1.00, 1; 0, 0.94];
A{7,1} = [0.94, 0.94; 0, 1];
A{8,1} = [1.00, 0.94; 0, 0.94];

%% Uncertain input matrix B
B{1,1} = [0; 1];
B{2,1} = [0; 1];
B{3,1} = [0; 1];
B{4,1} = [0; 1];
B{5,1} = [0; 1];
B{6,1} = [0; 1];
B{7,1} = [0; 1];
B{8,1} = [0; 1];

%% Uncertain output matrix C
C{1,1} = eye(2);
C{2,1} = eye(2);
C{3,1} = eye(2);
C{4,1} = eye(2);
C{5,1} = eye(2);
C{6,1} = eye(2);
C{7,1} = eye(2);
C{8,1} = eye(2);

%% Generate function outputs
model.A = A;
model.B = B;
model.C = C;

end % function