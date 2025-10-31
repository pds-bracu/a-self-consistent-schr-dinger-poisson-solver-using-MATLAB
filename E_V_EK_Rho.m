clc 
clear all 
close all 
  
%% 
q = 1.6e-19;    % charge of a proton 
NA = 1e24;      % doping concentration (m^-3) 
e0 = 8.854e-12; % permittivity of free space 
esio2 = 3.9*e0; % permittivity of SiO2 
esi = 11.7*e0;  % permittivity of Si 
VG = 1.2875;    % Gate voltage 
NV = 1.04e25;   % effective DOS in the valence band for Si 
(m^-3) 
  
nm = 1e-9;      % nanometer dimension 
tm = 6.25*nm;   % width of metal (nm) 
tsi = 25*nm;    % width of Si (nm) 
tox = 1.5*nm;   % width of SiO2 (nm) 
dz = 0.05 * nm; % differential unit (nm) 
z = [-tox:dz:0 dz:dz:(tsi-dz) tsi:dz:(tsi+tox)];   % without 
metal parts % z=0 is the 1st SiO2-Si junction 
len = length(z); 
nz = zeros(len,1);    % charge density (m^-3) 
c = 0;                % variable to count the # of iterations 
(not necessary) 
V = 0;                % potential profile 
  
for m=1:10000000 
    %% 
    c = c+1;          % counting the # of iterations 
    p = -(nz+NA);     % volume charge density, rho (m^-3) 
     
    I3 = (find(z==0));     % index of the left SiO2_Si 
junction 
    I4 = (find(z==tsi));   % index of the right SiO2_Si 
junction 
  
    A = zeros(len,len);    % left side matrix of the Poisson's 
eqn 
    A(1,1) = 1; 
    A(len,len) = 1; 
  
    B = zeros(len,1);      % right side matrix of the 
Poisson's eqn
B(1,1) = VG;           % Drichilet BC 
    B(len,1) = VG;         % Drichilet BC 
     
    for i = 2:len-1 
        if i==I3 || i==I4 
            B(i,1) = 0; 
        else 
            B(i,1) = -(q*p(i,1))/esi; 
        end 
    end 
  
    for i = 2:len-1 
        if i == I3      % Neumann BC at left SiO2-Si junction 
            A(i,i-1) = -esio2/dz; 
            A(i,i) = (esio2 + esi)/dz; 
            A(i,i+1) = -esi/dz; 
        elseif i == I4  % Neumann BC at the right SiO2-Si 
junction 
            A(i,i-1) = -esi/dz; 
            A(i,i)=(esio2 + esi)/dz; 
            A(i,i+1)= -esio2/dz; 
        else 
            A(i,i-1) = 1/(dz^2); 
            A(i,i) = -2/(dz^2); 
            A(i,i+1) = 1/(dz^2); 
        end 
    end 
  
    Vold = V;      % potential from the previous iteration 
    V = inv(A)*B;  % potential without the metal parts 
  
    % adding the constant potential on the both the sides for 
metal 
    znew = [-tm-tox:dz:-tox-dz z 
(tsi+tox+dz):dz:(tsi+tox+tm)]; 
  
    I1 = (find(znew==-tm-tox));   % index of leftmost point 
(metal) 
    I2 = (find(znew==-tox));      % index of the left Metal
SiO2 junction 
    I3 = (find(znew==0));         % new index of the left 
SiO2-Si junction 
    I4 = (find(znew==tsi));       % new index of the right 
SiO2-Si junction 
    I5 = (find(znew==tsi+tox));   % index of the right Metal
SiO2 junction 
    I6 = (find(znew==tsi+tox+tm));% index of rightmost point 
(metal) 
  
    lennew = length(znew);
     phi = zeros(lennew,1); % final potential profile with 
metal parts on the both sides 
    phi(I1:I2-1,1) = VG;   % adjusting the potential for metal 
on the left side 
    phi(I2:I5,1) = V;      % potential profile from Poisson's 
eqn found previously 
    phi(I5+1:I6,1) = VG;   % adjusting the potential for metal 
on the right side 
     
    if (m==1)              % plotting the Trial potential and 
initial E field 
        plot(znew/nm,phi,'Linewidth',2) 
        xlim([-10 35]); 
        xlabel('z (nm)'); 
        ylabel('Trial potential (Volts)'); 
        E_field = -gradient(phi);   % calculating E field from 
potential 
        set(gca,'FontSize',16); 
        figure() 
        plot(znew/nm,E_field,'Linewidth',2) 
        ylabel('Electric field (Vm^{-1})'); 
        xlim([znew(1)/nm znew(end)/nm]); 
        set(gca,'FontSize',16); 
    end 
     
    % adjusting the band diagram 
    E = -phi+VG;  %eV 
    E(I2:I5,1) = E(I2:I5,1) + 3.25; 
    E(I3:I4,1) = E(I3:I4,1) - 3.1; 
     
    if (m==1)    % plotting the initial band profile 
        figure() 
        plot(znew/nm,E,'Linewidth',1) 
        xlim([znew(1)/nm znew(end)/nm]); 
        ylabel('Energy (eV)'); 
        xlabel('z (nm)'); 
        set(gca,'FontSize',16); 
        hold on 
    end 
  
    %% 
    if abs(sum(Vold)-sum(V)) < 0.0001  % setting threshold 
        break 
    end 
    a = 0.95; 
    V = a*Vold + (1-a)*V;   % Gauss-Seidal convergence
    %% Schrodinger part 
     
    h = 6.626e-34; 
    h_ = h/(2*pi); 
    m0 = 9.1e-31; 
    m = 1.08*m0; 
  
    % adjusting the band diagram, this time without the metal 
portions 
    E_wo_metal = -V+VG; 
    E_wo_metal = E_wo_metal + 3.25; 
    E_wo_metal(find(z==0):find(z==tsi),1) = 
E_wo_metal(find(z==0):find(z==tsi),1) - 3.1; 
  
    % forming the Hamiltonial matrix 
    H0 = diag(((h_^2)/(m*dz^2))+E_wo_metal*q); 
    H1 = diag((ones(1,length(E_wo_metal)-1))*(
(h_^2)/(2*m*dz^2)),1); 
    H2 = diag((ones(1,length(E_wo_metal)-1))*(
(h_^2)/(2*m*dz^2)),-1); 
     
    H = H0+H1+H2; 
    [psi, En] = eig(H);   % eigenvalues 
      
    % normalizing psi 
    area=(dz*trapz(psi.^2)).^-1; 
    psi=psi.^2*diag(area); 
    psi=sqrt(psi); 
  
    %% charge density 
    nv = 6;     % degeneracy for Si 
    md = 0.358*m0;  % effective mass 
    K = 1.38e-23; 
    T = 300; 
    EF = q*(-1.04+VG+((K*T/q)*log(NV/NA)));   % Fermi energy 
  
    En1 = En(1:5,1:5);   % working with the first five 
wavefunctions 
    En1 = diag(En1); 
    psi1 = psi(:,1:5); 
     
    Nij = (nv*md*K*T)/(pi*h_*h_); 
    log_part = log(1+exp((EF-En)./(K*T))); 
    temp = (psi.^2)*diag(log_part); 
    nz = Nij*sum(temp,2); 
end 
  
plot(znew/nm,E,'r-.','Linewidth',2) 
xlim([znew(1)/nm znew(end)/nm]); 
legend('n(z)=0','n(z{\neq}0'); 
set(gca,'FontSize',16); 
figure() 
plot(znew/nm,phi,'Linewidth',2) 
xlim([-10 35]); 
xlabel('z (nm)'); 
ylabel('Actual potential (Volts)') 
set(gca,'FontSize',16); 
figure() 
plot(z/nm, nz,'Linewidth',2) 
hold on 
n_temp = zeros(length(nz),1); 
n_temp(find(z==0)+1:find(z==tsi),1) =  NA; 
plot(z/nm, nz+n_temp,'r--','Linewidth',2) 
xlim([z(1)/nm z(end)/nm]) 
xlabel('z (nm)'); 
ylabel('Inversion charge (m^{-2})'); 
legend('Injection charge','Total charge'); 
set(gca,'FontSize',16); 
