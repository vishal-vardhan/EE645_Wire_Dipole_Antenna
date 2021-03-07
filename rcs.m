%authors
% Bharath Thakkalapally
% Peela Jaswanth Aravind Kumar


%Tangential electric field on the surface is zero.So, Sum of incident an scattered Field is 0.

%assuming the incidents wave as TM signal wave (Hz = 0 and Ez varies with x and y)

% E_inc + E_scat = 0

% using the coeffs that were found in wireDipoleAntenna,We determine linear
% current density.

function [sigma_, alpha_] = rcs(lambda_,a_,N_,printData_)
    tic;
    eta0_ = 377;
    g_ = 1.78107; %gaama
    phii_ = pi;
    
    k_ = 2*pi/lambda_; % wavenumber
    w_ = 2*pi*a_/N_;  %segment size
    E_inc_ = zeros(N_,1); %incident electric field
    Z_ = zeros(N_,N_);

    for m_ = 1:N_
        for n_ = 1:N_
            rmn_ = abs(2*a_*sin(pi*(m_-n_)/N_));
            t_ = k_*rmn_;
            H_ = (1- (t_*t_/4)) - j*( ((2/pi)*log(g_*t_/2)) + (t_*t_/(2*pi)*(1 - log(g_*t_/2))) );
            Z_(m_,n_) = (2*pi/lambda_)*(eta0_/4)*(w_)*H_;
            
            if m_ == n_
                Z_(m_,n_) = (k_*eta0_*2*pi*a_/(4*N_))*(1-(j*(2/pi)*(log(g_*k_*2*pi*a_/(4*N_))-1)));
            end
        end
    end
    
    phi_ = linspace(0,N_-1,N_)*(2*pi/N_);
    for m_ = 1:N_
        E_inc_(m_) = exp(j*k_*(a_*cos(phi_(m_)-phii_)));
    end
    alpha_ = pinv(Z_)*E_inc_;
    
    sigma_ = zeros(N_,1);
    
    for m_=1:N_
        temp_ = 0;
        an_ = phi_(m_);
        for n_ = 1:N_
            ang_ = k_*a_*cos(phi_(n_)-an_);
            temp_ = temp_ + (alpha_(n_)*w_*(cos(ang_) + j*sin(ang_)));
        end
        temp_ = temp_* temp_;
        
        sigma_(m_) = k_*eta0_*eta0_*temp_/4;
    end
%     polarplot(phi_,abs(sigma_));
    if(printData_ == "Yes")
        plot(phi_*(180/pi),abs(alpha_));
        title('current distribution vs phi');
        xlabel('phi in degrees');
        ylabel('magnitude of current distribution');
    end
    toc;
    

end
