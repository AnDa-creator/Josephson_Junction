%Fig.21.1
clear all
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;qh=q/hbar; B=0;
%inputs
a=2.5e-9;t0=1;
NW=25;Nx=15;L=zeros(Nx);R=L;L(1,1)=1;R(Nx,Nx)=1;zplus=1i*1e-12;
%Hamiltonian
al=4*t0; by=-t0;bx=-t0;
alpha=kron(eye(NW),al)+kron(diag(ones(1,NW-1),+1),by)+kron(diag(ones(1,NW-1),-1),by');
alpha=alpha+diag([1:1:NW].*0);
alpha=alpha+diag([zeros(1,8) 0*ones(1,9) zeros(1,8)]);
beta=kron(diag(exp(i*qh*B*a*a*[1:1:NW])),bx);
H=kron(eye(Nx),alpha);
if Nx>1
H=H+kron(diag(ones(1,Np-1),+1),beta)+kron(diag(ones(1,Np-1),-1),beta');end
ii=0;for EE=[-0.05:1e-2:1.05]*t0
ii=ii+1;ig0=(EE+zplus)*eye(NW)-alpha; if ii==1
gs1=inv(ig0);gs2=inv(ig0);end
change=1;
while change >1e-6
Gs=inv(ig0-beta'*gs1*beta);
change=sum(sum(abs(Gs-gs1)))/(sum(sum(abs(gs1)+abs(Gs))));
gs1=0.5*Gs+0.5*gs1;
end
sig1=beta'*gs1*beta;sig1=kron(L,sig1);gam1=i*(sig1-sig1');
change=1;
while change >1e-6
Gs=inv(ig0-beta*gs2*beta');
change=sum(sum(abs(Gs-gs2)))/(sum(sum(abs(gs2)+abs(Gs))));
gs2=0.5*Gs+0.5*gs2;
end
sig2=beta*gs2*beta';sig2=kron(R,sig2);gam2=i*(sig2-sig2');
G=inv((EE*eye(Nx*NW))-H-sig1-sig2);
DD=real(diag(i*(G-G')))./2/pi;
Tcoh(ii)=real(trace(gam1*G*gam2*G'));E(ii)=EE/t0;ii
end
ii=1;for kk=pi*[-1:0.01:1]
H=alpha+beta*exp(i*kk)+beta'*exp(-i*kk);
[V,D]=eig(H);EK(:,ii)=sort(abs(diag(D)))./t0;K(ii)=kk/pi;ii=ii+1;
end
X=linspace(0,9,101);Ean= 2*(1-cos(pi*X./(NW+1)));
hold on
figure(1)
h=plot(Tcoh,E,'k');
set(h,'linewidth',[3.0])
%h=plot(X,Ean,'k--');
%set(h,'linewidth',[1.2])
set(gca,'Fontsize',[36])
axis([0 10 -.1 1])
% Fig.21.1a, Transmission versus width, at E=t0
clear all
%Constants (all MKS, except energy which is in eV)
hbar=1.06e-34;q=1.6e-19;qh=q/hbar;B=0;
%inputs
a=2.5e-9;t0=1;
NW=25;Nx=1;L=zeros(Nx);R=L;L(1,1)=1;R(Nx,Nx)=1;zplus=i*1e-12;
%Hamiltonian
al=4*t0;by=-t0;bx=-t0;
alpha1=kron(eye(NW),al)+kron(diag(ones(1,NW-1),+1),by)+kron(diag(ones(1,NW-1),-1),by');
ii=0;EE=t0*1;
for NN=[0:1:NW-1]
    ii=ii+1;
    alpha=alpha1+diag([zeros(1,NN) 100*ones(1,NW-NN)]);
    beta=kron(diag(exp(i*qh*B*a*a*[1:1:NW])),bx);
    H=kron(eye(Nx),alpha); 
if Nx>1
H=H+kron(diag(ones(1,Np-1),+1),beta)+kron(diag(ones(1,Np-1),-1),beta');
end
ig0=(EE+zplus)*eye(NW)-alpha; if ii==1
gs1=inv(ig0);gs2=inv(ig0);
end 
change=1;
while change >1e-6
Gs=inv(ig0-beta'*gs1*beta);
change=sum(sum(abs(Gs-gs1)))/(sum(sum(abs(gs1)+abs(Gs))));
gs1=0.5*Gs+0.5*gs1;
end
sigl=beta'*gsl*beta;sigl=kron(L,sigl);gaml=i*(sigl-sigl');
change=l;
while change >le-6
Gs=inv(ig0-beta*gs2*beta');
change=sum(sum(abs(Gs-gs2)))/(sum(sum(abs(gs2)+abs(Gs))));
gs2=0.5*Gs+0.5*gs2;
end
sig2=beta*gs2*beta';sig2=kron(R,sig2);gam2=i*(sig2-sig2');

G=inv((EE*eye(Np*NW))-H-sig1 -sig2);

DD=real(diag(i*(G-G')))./2/pi; Tcoh(ii)=real(trace(gaml*G*gam2*G')) ;E(ii)=NN; X(ii)=(NN+l)*(acos(l-(EE/2/tO)))/pi;
Xl(ii)=(NN+l)*sqrt(EE/2/tO)/pi;ii
end
hold on
figure(l)
h=plot(E,Tcoh,'k');
set(h,'linewidth',[3.0])
%h=plot(E,X,'k');
%h=plot(E,Xl,'k-');
%set(h,'linewidth',[1.2])
set(gca,'Fontsize',[36])
axis([0 NW -l 10])
grid on