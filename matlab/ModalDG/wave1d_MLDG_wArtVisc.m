%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Solving 1-D wave equation with Modal Local-DG discret.
% 					   with Artificial viscosity
%
%               du/dt + df/dx = e*u_xx,  for x \in [a,b]
%                  where f = f(u): linear & e = e(x)
%
%              coded by Manuel Diaz, NTU, 2013.10.13
%                 Copyright: all rights reserverd.
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: 
% 1.TVB Runge-Kutta Local Projection Discontinuous Galerkin Finite
% Element Method for conservation laws II: General Framework. (1989)
% 2. Tim Warburton, Numerical Partial Differential Equations, Lectures
% Notes. MA578-Section 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Basic Scheme Implementation with SSP-RK45 intergration scheeme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; %close all; clc;

%% Parameters
fluxfun = 'linear'; % select flux function
    cfl = 0.1;      % CFL condition
   tEnd = 1.0;      % final time
      P = 8;        % degree of accuaracy
     nE = 10;       % number of elements

% Define our Flux function
switch fluxfun
    case 'linear'
        a=-1; flux = @(w) a*w; 
        dflux = @(w) a*ones(size(w));
    otherwise
        error('case not listed')
end

% Build 1d mesh
xgrid = mesh1d([0,1],nE,'ChebyshevMod',P);
dx = xgrid.elementSize; J = xgrid.Jacobian; 
x = xgrid.nodeCoordinates;	xc = xgrid.elementCenter;

% Load DG tools
tool = DGtools(xgrid.solutionPoints);
V = tool.Vandermonde;
lR = tool.legRightEnd; lL = tool.legLeftEnd;
D = tool.MassMatrix*tool.CoefDiffMatrix;
invM = tool.invMassMatrix;

% Litf operators
F =  lL'*lL;  % F_nm =  L_n(-1)*L_m(-1)
G = -lL'*lR;  % G_nm = -L_n(-1)*L_m(+1)
H = -lR'*lR;  % H_nm = -L_n(+1)*L_m(+1)
I =  lR'*lL;  % J_nm =  L_n(+1)*L_m(-1)

% IC and k(x)
q0 = IC(x,6); %square Jump

% Set plot range
plotrange = [xgrid.range(1),xgrid.range(2),...
    min(min(min(0.9*q0)),min(min(1.1*q0))),1.1*max(max(q0))];

%% Solver Loop

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Set initial time & load IC
t = 0; q = q0; qt = V\q; it = 0;

% Using a 4rd Order 5-stage SSPRK time integration
res_ut = zeros(P+1,nE); % Runge-Kutta residual storage

% Using a 3-stage TVD Runge Kutta time integration
while t < tEnd
    qo=qt;
    
    % update time
    dt=cfl*dx^2; t=t+dt;
    
    % iteration counter
    it=it+1; 
    
    % Artificial viscosity
    epsilon = Smoothfunc(qt,V,dx,P);
    
    for RKs = 1:5
        t_local = t + rk4c(RKs)*dt;
        dF = RHS4(a,qt,D,J,F,G,H,I,invM,epsilon,P);
        res_ut = rk4a(RKs)*res_ut - dt*dF;
        qt = qt - rk4b(RKs)*res_ut;
    end
    
    % Transform legendre coefs into nodal values.
    q = V*qt;
    
    % build cell averages
    q_bar = qt(1,:);
    
    % Plot u 
    if rem(it,10) == 0
        subplot(1,2,1); 
        plot(x,q,x,q0,'-+'); axis(plotrange); grid on;
        xlabel('x'); ylabel('u'); title('DG-FEM')
        subplot(1,2,2); 
        plot(xc,q_bar,'ro'); axis(plotrange); grid on;
        drawnow;
    end
end

%% Final Plot for IC 2
subplot(1,2,1); plot(x,q,x,q0,'-+'); axis(plotrange);
title('Modal LDG','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);
subplot(1,2,2); plot(x,q0,'k-',xc,q_bar,'ro'); axis(plotrange);
title('Cell Averages','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);