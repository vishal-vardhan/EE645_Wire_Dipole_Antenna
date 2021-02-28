%authors
% Bharath Thakkalapally
% Peela Jaswanth Aravind Kumar


%Tangential electric field on the surface is zero.So, Sum of incident an scattered Field is 0.

%assuming the incidents wave as TM signal wave (Hz = 0 and Ez varies with x and y)

% E_inc + E_scat = 0

% using the coeffs that were found in wireDipoleAntenna,We determine linear
% current density.

function [sigma, alpha] = rcs(lambda,a,N,printData)
    tic;
    eta0 = 377;
    g = 1.78107; %gaama
    phii = pi;
    
    k = 2*pi/lambda; % wavenumber
    w = 2*pi*a/N;  %segment size
    E_inc = zeros(N,1); %incident electric field
    Z = zeros(N,N);

    for m = 1:N
        for n = 1:N
            rmn = abs(2*a*sin(pi*(m-n)/N));
            t = k*rmn;
            H = (1- (t*t/4)) - j*( ((2/pi)*log(g*t/2)) + (t*t/(2*pi)*(1 - log(g*t/2))) );
            Z(m,n) = (2*pi/lambda)*(eta0/4)*(w)*H;
            
            if m==n
                Z(m,n) = (k*eta0*2*pi*a/(4*N))*(1-(j*(2/pi)*(log(g*k*2*pi*a/(4*N))-1)));
            end
        end
    end
    
    phi = linspace(0,N-1,N)*(2*pi/N);
    for m = 1:N
        E_inc(m) = exp(j*k*(a*cos(phi(m)-phii)));
    end
    alpha = pinv(Z)*E_inc;
    
    sigma = zeros(N,1);
    
    for m=1:N
        temp = 0;
        a = phi(m);
        for n = 1:N
            ang = k*a*cos(phi(n)-a);
            temp = temp + (alpha(n)*w*(cos(ang) + j*sin(ang)));
        end
        temp = temp* temp;
        
        sigma(m) = k*eta0*eta0*temp/4;
    end
%     polarplot(phi,abs(sigma));
    if(printData == "Yes")
        plot(phi*(180/pi),abs(alpha));
        title('current distribution vs phi');
        xlabel('phi in degrees');
        ylabel('magnitude of current distribution');
    end
    toc;
    

end