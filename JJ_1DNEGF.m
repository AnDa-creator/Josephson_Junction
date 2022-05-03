clear all
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;qh=q/hbar; B=0; kT = 0.025;                    % Lattice Constants
%inputs
a=2.5e-9; t=1; NxL=15; NxD = 15; NxR = 15;N = NxL + NxD + NxR;  % Lattice Parameters
L=zeros(N);                                      
zplus=1i*1e-12;                                                        % Infinitesimal eta

% Levels and s-wave junction parameters
mu_f = 1.5 * t; del_d = t * 0.05; del0 = del_d * mu_f/t;                       % Pairing parameters
Emax = 5 * del0; Emin = - 5 * del0; Ne = 30; E = linspace(Emin, Emax, Ne);
vol = 0; 
mu1 = mu_f + vol/2; mu2 = mu_f + -vol/2;
mu = [mu1 * ones(1, NxL), linspace(mu1, mu2, NxD), mu2 * ones(1, NxR)];     % Chemical potential profile
f1=1./(1+exp((E-mu1)./kT));
f2=1./(1+exp((E-mu2)./kT));

ep = 0.2;                                                                   % onsite energy
alpha = zeros(2, 2*N);                                                        % onsite interaction
for i = 1:N
    alpha(:,2*i - 1:2*i) = [ep - mu(i), del_d; del_d, -(ep - mu(i))];
end
alpha_1 = alpha(:,1:2);                                                         % onsite interaction for lead 1
alpha_2 = alpha(:,2*N-1:2*N);                                                       % onsite interaction for lead 2

beta = [-t, 0; 0, t];                                                       % hopping parameters

H_bdg = kron(diag(ones(1,N-1),1), beta) + kron(diag(ones(1,N-1), -1), beta);  % TB Hamiltonian 
for i = 1:N
    mat = zeros(1,N);
    mat(i) = 1;
    H_bdg = H_bdg + kron(diag(mat) , alpha(:,2*i - 1:2*i));
end

I = 0;

for j = 1:Ne
    gs1=inv((E(j) + zplus)*eye(2,2));gs2=inv((E(j) + zplus)*eye(2,2));
    % Green's function for Lead 1
    change=1;
    while change >1e-6
        g1=inv((E(j) + zplus)*eye(2,2) - alpha_1  -beta'*gs1*beta);
        change=sum(sum(abs(g1-gs1)))/(sum(sum(abs(g1)+abs(gs1))));
        gs1=0.5*g1+0.5*gs1;
    end
    % Green's function for Lead 2
    change=1;
    while change >1e-6
        g2=inv((E(j) + zplus)*eye(2,2) - alpha_2  -beta*gs2*beta');
        change=sum(sum(abs(g2-gs2)))/(sum(sum(abs(g2)+abs(gs2))));
        gs2=0.5*g2+0.5*gs2;
    end
    sig_1 = zeros(2*N,2*N); sig_2 = zeros(2*N,2*N);
    sig_1(1:2,1:2) = beta' * gs1 * beta; sig_2((2*N - 1):2*N,(2*N - 1):2*N) = beta' * gs2 * beta;
    
    G_ret = inv((E(j) + zplus)*eye(2*N,2*N) - sig_1 - sig_2);

    sigLess_1 = 1i * (sig_1 - sig_1') * f1(j);
    sigLess_2 = 1i * (sig_2 - sig_2') * f2(j);

    G_less = G_ret * (sigLess_1 + sigLess_2) * G_ret';   

    I = I + (q/(hbar * 2*pi));

end


