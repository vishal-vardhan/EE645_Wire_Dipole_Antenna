% author : Adepu Vishal Vardhan (170108003)

clear all;
close all;
clc;

Nz = 200;  % number of spacial segments
ks = 1;  % source index
Nt = 800;  % time steps
c0 = 299792458;
eps0 = 8.854e-12;
mu0 = 1/(c0*c0*eps0);

% excitation
to = 20;  % mean of the gaussian excitation
spread = 8;  % sigma of gaussian excitation
fs = 0;  % frequency of excitation
% fs = 1.5; 

% step size
dz = 0.001;  % space step size
r = 1;  % factor for changing magic time step condition
dt = r*dz/c0;  % time step size

% permittivity variation
er = ones(1, Nz);
sigmaE = zeros(1, Nz);
sigmaM = zeros(1, Nz);
mur = ones(1, Nz);

er(Nz/2:Nz) = 4;  % change in dielectric media
sigmaE(Nz/2:Nz) = 0.0;  % lossy dielectric media

Ex = zeros(1, Nz);  % electric field
Hy = zeros(1, Nz);  % magnetic field intensity

A = zeros(1, Nz);
B = zeros(1, Nz);
C = zeros(1, Nz);
D = zeros(1, Nz);
Ce = dt/(2*eps0);
Ch = dt/(2*mu0);

A = (er-Ce*sigmaE)./(er+Ce*sigmaE);
B = (2*Ce)./(er+Ce*sigmaE);
C = (mur-Ch*sigmaM)./(mur+Ch*sigmaM);
D = (2*Ch)./(mur+Ch*sigmaM);

figure;
for t = 1:Nt
	% mur absorbing boundary condition
    Ex(1) = Ex(2)+((c0*dt/sqrt(mur(1)*er(1))-dz)/(c0*dt/sqrt(mur(1)*er(1))+dz))*(Ex(2)-Ex(1));
    Ex(Nz) = Ex(Nz-1)+((c0*dt/sqrt(mur(Nz)*er(Nz))-dz)/(c0*dt/sqrt(mur(Nz)*er(Nz))+dz))*(Ex(Nz-1)-Ex(Nz));

	% field update equations
    Ex(2:Nz-1) = A(2:Nz-1).*Ex(2:Nz-1)-B(2:Nz-1).*((Hy(2:Nz-1)-Hy(1:Nz-2))/dz);
    
	% apply source
    Ex(ks) = exp(-0.5*((t-to)/spread)^2)*cos(2*pi*fs*t);
    Hy(1:Nz-1) = C(1:Nz-1).*Hy(1:Nz-1)-D(1:Nz-1).*((Ex(2:Nz)-Ex(1:Nz-1))/dz);
   
    % plot electric field
    subplot(2, 1, 1);
    plot(Ex);
    ylabel('Ex -->');
    xlabel('z -->');
    title(['time step = ', num2str(t)]);
    % plot magnetic field
    subplot(2, 1, 2);
    plot(Hy);
    ylabel('Hy -->');
    xlabel('z -->');
    pause(0);
end


