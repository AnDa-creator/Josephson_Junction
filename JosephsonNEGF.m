% d_wave junction
clear all; close all;
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;qh=q/hbar; B=0; kT = 0.025;                    % Lattice Constants
%inputs
a=2.5e-9; t=1; Ny=25; Nx=15;L=zeros(Nx);                                      % Lattice Parameters
zplus=1i*1e-12;                                                                                  % Infinitesimal eta

% Levels and s-wave junction parameters
mu_f = 1.5 * t; del_d = t * 0.05; del0 = del_d * mu_f/t;                       % Pairing parameters
theta_1 = pi/4; theta_2 = -pi/4; 
phi_1 = 0; phi_2 = 0;
kf = acos(1 - mu_f/(2*t))/a; Nky = 10; ky = linspace(-kf, kf, Nky);
Emax = 5 * del0; Emin = - 5 * del0; Ne = 20; E = linspace(Emin, Emax, Ne);
dE = E(2) - E(1);
T0 = 1;
U0 = 2 * t * sin(kf * a) * sqrt(1/T0 - 1);
vol = 0; 
mu1 = vol/2; mu2 = -vol/2;
f1=1./(1+exp((E-mu1)./kT));
f2=1./(1+exp((E-mu2)./kT));
I0 = (2 * q * t * dE)/(hbar * 2* pi * Nky);

N_p = 100;
phi_mat = linspace(0, 2*pi, N_p); 
Ij_mat = zeros(1, N_p);
ppm = ParforProgressbar(N_p);

parfor k = 1:N_p

    phi = phi_mat(k);k
    I = 0;
    for i = 1:Nky
        % alpha matrices for the leads
        alpha_1 = zeros(2,2);
        alpha_1(1,1) = 4*t - 2*t*cos(ky(i) * a) - mu_f - mu1; alpha_1(1,2) = del_d *exp(-1i * phi)* cos(2 * theta_1) * cos(ky(i)*a);
        alpha_1(2,1) = alpha_1(1,2)'; alpha_1(2,2) = - 4 * t + 2 * t * cos(ky(i) * a) + mu_f - mu1;
    
        alpha_2 = zeros(2,2);
        alpha_2(1,1) = 4*t - 2*t*cos(ky(i) * a) - mu_f - mu2; alpha_2(1,2) = del_d * cos(2 * theta_2) * cos(ky(i)*a);
        alpha_2(2,1) = alpha_2(1,2)'; alpha_2(2,2) = - 4 * t + 2 * t * cos(ky(i) * a) + mu_f - mu2;
    
        % Beta matrices for the leads
        beta_1 = zeros(2,2);
        beta_1(1,1) = -t; 
        beta_1(1,2) = 0.5 * (-del_d * exp(-1i * phi) * cos(2 * theta_1) - 1i * sin(2*theta_1) * sin(ky(i)*a));
        beta_1(2,1) = beta_1(1,2)';
        beta_1(2,2) = t;
    
        beta_2 = zeros(2,2);
        beta_2(1,1) = -t;
        beta_2(1,2) = 0.5 * (-del_d * cos(2 * theta_2) - 1i * sin(2*theta_2) * sin(ky(i)*a));
        beta_2(2,1) = beta_2(1,2)';
        beta_2(2,2) = t;
        
    
        for j = 1:Ne
            gs1=inv((E(j) + zplus)*eye(2,2));gs2=inv((E(j) + zplus)*eye(2,2));
            g1 = gs1; g2 = gs2;
            % Green's function for Lead 1
            change=1;
            while change >1e-6
                g1=inv((E(j) + zplus)*eye(2,2) - alpha_1  -beta_1'*gs1*beta_1);
                change=sum(sum(abs(g1-gs1)))/(sum(sum(abs(g1)+abs(gs1))));
                gs1=0.5*g1+0.5*gs1;
            end
            % Green's function for Lead 2
            change=1;
            while change >1e-6
                g2=inv((E(j) + zplus)*eye(2,2) - alpha_2  -beta_2'*gs2*beta_2);
                change=sum(sum(abs(g2-gs2)))/(sum(sum(abs(g2)+abs(gs2))));
                gs2=0.5*g2+0.5*gs2;
            end
    
            % Coupling Matrices
    
            tau_1 = [-t 0; 0 t];
            tau_2 = [-t 0; 0 t];
    
            
            % Calculate the device green function
            g_d = zeros(2,2);
            g_d(1,1) = E(j) - 4 * t + 2* t * cos(ky(i) * a) - U0 + mu_f;
            g_d(2,2) = E(j) + 4 * t - 2* t * cos(ky(i) * a) + U0 - mu_f;
            g_d = inv(g_d);
    
            % Overall Green's function
            GD = zeros(2+2+2,2+2+2);
            GD(1:2,1:2) = inv(g1); GD(1:2,3:4) = -tau_1;
            GD(3:4,1:2) = -tau_1'; GD(3:4,3:4) = inv(g_d);
            GD(5:6,5:6) = inv(g2); GD(5:6,3:4) = -tau_2;
            GD(3:4,5:6) = -tau_2';
            
            GD = inv(GD);
    
            % Calculate the Current
            
            I = I + I0 * (real(GD(3,1) - GD(1,3) + GD(4,2) - GD(2,4))) * (f1(j));
        
        end

    end
    
    Ij_mat(k) = I;
    pause(100/N_p);
    ppm.increment(); 
end
delete(ppm);
%% Plots
figure(6)                                                                                                   % Subpart c                                             
h = plot(phi_mat/(2*pi),Ij_mat);
hold on
g = plot(phi_mat/(2*pi), max(Ij_mat) * -sin(phi_mat));
set(h,'linewidth',2.0)
set(gca,'Fontsize',24)
ylabel("current(arb units)");
xlabel("\Phi (rad)");
xlim padded
grid on
legend('Current', 'sinusoid')