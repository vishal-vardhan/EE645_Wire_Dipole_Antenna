clc
clear all
close all

ep0 = 8.854e-12;
mu = 4*pi*1e-7;
er = 1;
ep =er*ep0;
c = 3.8*1e8;
dz = 0.001;
dx = 0.001;
dt = dz/c;

n=100;
ey = zeros(n,n);
hz = zeros(n,n);
hx = zeros(n,n);

sig = 8;
nsteps = 500;
cz = n/2;
cs = nsteps/5;
cx = n/2;

%intializing electric field with guassian funtion
% for nx=1:n
%     for nz = 1:n
%         ey(nx,nz) = exp(-1*((nx-cx)^2 + (nz-cz)^2)/(2*sig*sig));
% %         hz(nx,nz) = exp(-1*((nx-cx)^2 + (nz-cz)^2)/(2*sig*sig));
% %         hx(nx,nz) = exp(-1*((nx-cx)^2 + (nz-cz)^2)/(2*sig*sig));
%     end
% end


%this code is for TE mode

figure(1)
for i =1:nsteps
    for nx=2:n
        ey(nx,1) = ey(nx,2);
        ey(nx,n) = ey(nx,n-1);
    end
    for nz=2:n
        ey(1,nz) = ey(2,nz);
        ey(n,nz) = ey(n-1,nz);
    end
    for nx=2:n-1
        ey(nx,1) = exp(-1*((i-cs)^2 + (nx-cx)^2)/(2*sig*sig));
    end
    for nx =1:n-1
        for nz =1:n-1
            ey(nx,nz) = ey(nx,nz) + (dt/(ep*dz))*(hx(nx,nz+1)-hx(nx,nz)) - (dt/ep*dx)*(hz(nx+1,nz)-hz(nx,nz));
        end
    end
    
    for nx =2:n
        for nz =2:n
            hx(nx,nz) = hx(nx,nz) + (dt/(mu*dz))*(ey(nx,nz)-ey(nx,nz-1));
            hz(nx,nz) = hz(nx,nz) - (dt/(mu*dx))*(ey(nx,nz)-ey(nx-1,nz));
        end
    end
    
    mesh(ey);
%     mesh(ey);
%     mesh(hz);
    drawnow;
    pause(0.01);
    
end



