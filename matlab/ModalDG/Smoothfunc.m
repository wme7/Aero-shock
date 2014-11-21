function [epsilon] = Smoothfunc(ut,V,h,P)
%SMOOTHFUNC: Discontinuity sensor
%   This function helpme to determine how smooth the the data inside every
%   element.

% u: nodal values
u=V*ut;

% u_hat: truncated ut
ut(P+1,:)=0; u_hat=V*ut;

% smooth indicator
S=sum((u-u_hat).*(u-u_hat)./(u.*u));

% s variable
s = log10(S); % s_e
s(s==-Inf)=0; % get rid of NaN or Inf values

% parameters
k = 5.0;    % $\kappa$: empirically choosen to obtain a sharp profile
so = 1/P^4; % $s_0$
epso = h/P; % $epsilon_0$

% ranges
range1 = s<so-k;
range2 = so-k<s & s<so+k;
range3 = s>so+k;

% epsilon
epsilon = 0*range1 + epso/2*(1+sin(pi*(s-so)/(2*k))).*range2 + epso.*range3;