function dF = RHS4(a,uto,D,J,F,G,H,I,invM,k,P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 1d heat equation using Nodal DG 
%
%         du/dt = -df/dx + epsilon*d^2u/dx^2,  for x \in [a,b]
%
%              coded by Manuel Diaz, NTU, 2013.12.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following ideas in: 
% Tim Warburton, Numerical Partial Differential Equations, Lectures
% Notes. MA578-Section 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lax-Friedrich switch parameters for advection
phiL = (1/a)*(a+abs(a))/2;
phiR = (1/a)*(a-abs(a))/2;

utm=circshift(uto,[0,+1]); % ut(i-1)
utp=circshift(uto,[0,-1]); % ut(i+1)

% Apply Periodic BCs
% just comment all BCs

% Apply Dirichlet BCs
%utm(:,1) = -utp(:,1);     % left BD
%utp(:,end) = -utm(:,end); % right BD

% Apply Neumann BCs1
% utm(:,1) = utp(:,1);     % left BD
% utp(:,end) = utm(:,end); % right BD

% Compute ut_x variable in LDG
ut_x = invM*(D*uto + phiL*(F*uto + G*utm) + phiR*(H*uto + I*utp))/J;

% Compute q variable in LDG
qto = invM*(D*uto + 0.5*(F*uto + G*utm) + 0.5*(H*uto + I*utp))/J;

qtm=circshift(qto,[0,+1]); % qt(i-1)
qtp=circshift(qto,[0,-1]); % qt(i+1)

% Apply Periodic BCs
% just comment all BCs

% Apply Dirichlet BCs
%qtm(:,1) = qtp(:,1);     % left BD
%qtp(:,end) = qtm(:,end); % right BD

% Apply Neumann BCs in this step
% qtm(:,1) = -qtp(:,1);     % left BD
% qtp(:,end) = -qtm(:,end); % right BD

% Compute the 2nd derivate
ut_xx = invM*(D*qto + 0.5*(F*qto + G*qtm) + 0.5*(H*qto + I*qtp))/J;

% residual
dF = -a*ut_x + ones(P+1,1)*k.*ut_xx;