%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE DATA OF THE SOLUTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('solution.dat','r');
format long e  
p1=fscanf(fid,'%le %le %le %le %le',[5 inf]);
x = p1(1,:); r = p1(2,:); u = p1(3,:); P = p1(4,:); s=p1(5,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
hold on; axis('square'); hold on; grid; title('Density');
xlabel('x'); ylabel('Density \rho'); plot(x,r,'ko');
%text(0.01,0.96,'Van Leer');text(0.01,0.80,'T=1.74423');
%text(0.01,3.9,'Osher');text(0.01,3.7,'29 time Steps');
%text(0.01,3.5,'CFL=0.95');  %S-W t=0.70471
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,2)
hold on; axis('square'); hold on; grid; title('Velocity')
xlabel('x'); ylabel('Velocity u'); plot(x,u,'ko');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3)
hold on; axis('square'); hold on; grid; title('Pressure')
xlabel('x'); ylabel('Pressure P'); plot(x,P,'ko');
%text(0.01,14.5,'Osher');text(0.01,13.3,'29 time Steps');
%text(0.01,12,'CFL=0.95');  %S-W t=0.7022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4)
hold on; axis('square'); hold on; grid; title('Entropy/R_{gas}')
xlabel('x'); ylabel('Entropy/R_{gas}'); plot(x,s,'ko');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
