close('all')
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x, rho, u, pres, T, tau, q]...
  = textread('ns_shock_structure.dat','%f %f %f %f %f %f %f','headerlines',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)

subplot(2,3,1)
hold on; axis('square'); hold on; grid; title('\rho');
xlabel('x'); ylabel('Density \rho');
plot(x,rho,'r*');

subplot(2,3,2)
hold on; axis('square'); hold on; grid; title('u')
xlabel('x'); ylabel('Velocity u');
plot(x,u,'r*');

subplot(2,3,3)
hold on; axis('square'); hold on; grid; title('p')
xlabel('x'); ylabel('Pressure p');
plot(x,pres,'r*');

subplot(2,3,4)
hold on; axis('square'); hold on; grid; title('T')
xlabel('x'); ylabel('Temperature T');
plot(x,T,'r*');

subplot(2,3,5)
hold on; axis('square'); hold on; grid; title('\tau')
xlabel('x'); ylabel('Viscous stress \tau');
plot(x,tau,'r*');

subplot(2,3,6)
hold on; axis('square'); hold on; grid; title('q')
xlabel('x'); ylabel('Heat flux q');
plot(x,q,'r*');
