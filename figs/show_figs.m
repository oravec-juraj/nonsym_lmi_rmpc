%% Non & Sym results

close all
nk = 52; % Number of sim-steps

%% Constraints
%% X1
figure(1)
axis([0,nk-2,-2,18])
aX1 = set_axis(-2, 1, 18, -2, 2, 18);
set_labels(gca,'y',[-2:1:18],aX1)
plot([0,nk],[0,0],'r:')
xlabel('t')
ylabel('y1')

%% X2
figure(2)
axis([0,nk-2,-2,3])
aX2 = set_axis(-2, 0.5, 3, -2, 1, 3);
set_labels(gca,'y',[-2:0.5:3],aX2)
plot([0,nk],[0,0],'r:')
xlabel('t')
ylabel('y2')

%% FIGURE MV
figure(3)
axis([0,nk-2,-1.5,0.5])
aU = set_axis(-1.5,0.5,0.5,-1.5,0.5,0.5);
set_labels(gca,'y',[-1.5:0.5:0.5],aU)
plot([0,nk],[-u_satL,-u_satL],'r-.')
plot([0,nk],[u_satU,u_satU],'r--')
plot([0,nk],[-u_satU,-u_satU],'r--')
xlabel('t')
ylabel('u')

%% Non-symmetric
load('results_nonsym.mat')

for v = 1 : nv
x = data.x{v};
u = data.u{v};
F = data.F{v};

%% FIGURE X1
figure(1)
hold on, box on
plot(0:nk-2,x(1,2:nk),'-','color',[0.2 0.7 0])
xlabel('t')
ylabel('x_1')

%% FIGURE X2
figure(2)
hold on, box on
plot(0:nk-2,x(2,2:nk),'-','color',[0.2 0.7 0])
xlabel('t')
ylabel('x_2')

%% FIGURE MV
figure(3)
hold on, box on
plot(0:nk-1,u(1,1:nk),'-','color',[0.2 0.7 0])
xlabel('t')
ylabel('U')

%% FIGURE NORM_F
figure(4), box on
for k = 1 : simstep
    norm_F(k,1) = norm(F{k});
end % for k
plot(1:k,norm_F,'s')
xlabel('t')
ylabel('|| F ||_2')
%
end % for v

%% Symmetric
load('results_symm.mat')

for v = 1 : nv
x = data.x{v};
u = data.u{v};
F = data.F{v};

%% FIGURE X1
figure(1)
plot(0:nk-2,x(1,2:nk),'b-')

%% FIGURE X2
figure(2)
plot(0:nk-2,x(2,2:nk),'b-')

%% FIGURE MV
figure(3)
plot(0:nk-1,u(1,1:nk),'b-')

%% FIGURE NORM_F
figure(4), box on
for k = 1 : simstep
    norm_F(k,1) = norm(F{k});
end % for k
plot(1:k,norm_F,'s')
xlabel('t')
ylabel('|| F ||_2')
%
end % for v

