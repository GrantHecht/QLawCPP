%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q law plotting script by Donald Ellison
% March 15th 2014
% Uses earth_sphere function Will Campbell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% Simulation parameters
DU = 6371.0;
TU = 8.068231200000000e+02;
mu = 1.0;

% Extract COE data from file
filestring = ['.\qlaw_output.txt'];
fileptr = fopen(filestring);
data = textscan(fileptr,'%f %f %f %f %f %f %f %f','headerlines',2);
%data = textscan(fileptr,'%f %f %f %f %f %f %f %f');

p = data{1};
e = data{2};
inc = data{3};
ape = data{4};
ran = data{5};
mass = data{6};
time = data{7};
L = data{8};



% Close data file
fclose(fileptr);

% Calculate true anomaly and semimajor axis
tru = L - ran - ape;
sma = p./(1.0 - e.*e);

% Convert time from TU's to seconds
time = time * TU;

% Convert from COE's to cartesian state
r = sma.*(1.0 - e.*e)./(1.0 + e.*cos(tru));
v = sqrt(2*mu./r - mu./sma);
h = sqrt(mu*sma.*(1.0 - e.*e));

x = r.*(cos(tru+ape).*cos(ran) - sin(tru+ape).*cos(inc).*sin(ran));
y = r.*(cos(tru+ape).*sin(ran) + sin(tru+ape).*cos(inc).*cos(ran));
z = r.*(sin(tru+ape).*sin(inc));

vx = -(mu./h).*(cos(ran).*(sin(tru+ape) + e.*sin(ape))+sin(ran).*(cos(tru+ape)+e.*cos(ape)).*cos(inc));
vy = -(mu./h).*(sin(ran).*(sin(tru+ape) + e.*sin(ape))-cos(ran).*(cos(tru+ape)+e.*cos(ape)).*cos(inc));
vz = (mu./h).*(cos(tru+ape)+e.*cos(ape)).*sin(inc);

% Plot some results
figure(1);
earth_sphere(100,'DU');
hold on;
traj = plot3(x,y,z,'.r','MarkerSize',2);
xlabel('x (DU)');
ylabel('y (DU)');
zlabel('z (DU)');
set(gca,'DataAspectRatio',[1 1 1]);

figure(2);
semimajor = plot(time./86400,sma*DU,'.b');
xlabel('time (days)');
ylabel('semimajor axis (km)');

figure(3);
eccentricity = plot(time./86400,e,'.b');
xlabel('time (days)');
ylabel('eccentricity');

figure(4);
inclination = plot(time./86400,inc*180/pi,'.b');
xlabel('time (days)');
ylabel('inclination (deg)');

ape = mod(ape,2*pi);
figure(5);
aperiapse = plot(time./86400,ape*180/pi,'.b');
xlabel('time (days)');
ylabel('argument of periapse (deg)');

ran = mod(ran,2*pi);
figure(6);
RAN = plot(time./86400,ran*180/pi,'.b');
xlabel('time (days)');
ylabel('RAN (deg)');

tru = mod(tru,2*pi);
figure(7);
trueanom = plot(time./86400,tru*180/pi,'.b');
xlabel('time (days)');
ylabel('true anomaly (deg)');

figure(8);
massplot = plot(time./86400,mass,'.b');
xlabel('time (days)');
ylabel('mass (kg)');

figure(9);
rper = plot(time./86400,sma.*(1 - e),'.b');
xlabel('time (days)');
ylabel('periapse radius (DU)');

