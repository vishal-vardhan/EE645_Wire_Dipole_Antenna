function [coeff_, Zin_] = wireDipoleAntenna(L_, a_, nSegments_, excitedSeg_, freq_, k_, V_, printData_)
    tic;
    %%%%%%%%%%%%%%%%%
    % input variables
    %%%%%%%%%%%%%%%%%
    % L = length of the wire
    % a = radius of the wire
    % nSegments_ = number of segments in the wire
    % freq = frequency
    % k = free space wave number
    % V = applied voltage
    % printData = flag to print output data

    %%%%%%%%%%%%%%%%%%
    % output variables
    %%%%%%%%%%%%%%%%%%
    % coeff_ = weights of the basis functions
    % Zin_ = input impedance of the wire
    
    % solving using gap generation model method
    e0_ = 8.85418*1e-12;  % permittivity of free space
    w_ = 2*pi*freq_;  % angular frequency
    delta_ = L_/nSegments_;  % step size
    Ezi_ = zeros(nSegments_, 1);  % incident em wave
    Ezi_(excitedSeg_) = V_/delta_;  % electric field at the excited segment
    z_ = zeros(nSegments_, 1);  % central points of each section

    z_(1) = (-L_/2)+(delta_/2);
    for m_ = 2:nSegments_
        z_(m_) = z_(m_-1)+delta_; 
    end

    A_ = zeros(nSegments_, nSegments_);
    syms x_;
    for m_ = 1:nSegments_
        for n_ = 1:nSegments_
            r_ = sqrt(a_^2+(z_(m_)-x_)^2);
            Ker_ = (exp(-1j*k_*r_)/(r_^5))*((1+1j*k_*r_)*(2*r_*r_-3*a_*a_)+((k_*a_*r_)^2));
            A_(m_, n_) = (1/(4*pi*w_*e0_*1j))*int(Ker_, x_, z_(n_)-delta_/2, z_(n_)+delta_/2);
        end
    end

    coeff_ = pinv(A_)*(-Ezi_);  % the weights of pulses
    Zin_ = V_/coeff_(excitedSeg_);
    
    % TODO add code to print data
%     if printData == 'y' || printData == 'Y'
%         % do the printing here
%     end
    toc;
end
