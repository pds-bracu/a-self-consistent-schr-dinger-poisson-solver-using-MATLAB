clc 
clear all 
close all 
c = 0; 
V = 0; 
VG_array = 0:0.64375/16:2.1;    
% Gate voltage array 
Index_VG_1_2875 = find(VG_array == 1.2875); 
Q = zeros(1,length(VG_array));   
for n=1:length(VG_array) 
for m=1:10000 
c = c+1; 
q = 1.6e-19;    
% Charge density 
% charge of a proton 
NA = 1e24;      
% doping concentration (m^-3) 
e0 = 8.854e-12; % permittivity of free space 
esio2 = 3.9*e0; % permittivity of SiO2 
esi = 11.7*e0;  % permittivity of Si 
%     VG = 1.2875; 
VG = VG_array(n); 
NV = 1.04e25;   
% effective DOS in the valence band 
for Si (m^-3) 
% 0 is the 1st SiO2-Si junction, so 
nm = 1e-9;      
% nanometer dimension
tm = 6.25*nm;   % width of metal (nm) 
        tsi = 25*nm;    % width of Si (nm) 
        tox = 1.5*nm;   % width of SiO2 (nm) 
        dz = 0.05 * nm; % differential unit (nm) 
        z = [-tox:dz:0 dz:dz:(tsi-dz) tsi:dz:(tsi+tox)];   % 
without metal parts % z=0 is the 1st SiO2-Si junction 
        len = length(z); 
         
        if (m == 1) 
            nz = zeros(len,length(VG_array)); 
        end 
        p = -(nz(:,n)+NA); % volume charge density, rho (m^-3) 
  
        I3 = (find(z==0)); % index of the left SiO2_Si 
junction 
        I4 = (find(z==tsi)); % index of the right SiO2_Si 
junction 
  
        A = zeros(len,len); % left side matrix of the 
Poisson's eqn 
        A(1,1) = 1; 
        A(len,len) = 1; 
  
        B = zeros(len,1); % right side matrix of the Poisson's 
eqn 
        B(1,1) = VG;      % Drichilet BC 
        B(len,1) = VG;    % Drichilet BC 
  
        for i = 2:len-1 
            if i==I3 || i==I4 
                B(i,1) = 0; 
            else 
                B(i,1) = -(q*p(i,1))/esi; 
            end 
        end 
  
        for i = 2:len-1 
            if i == I3   % Neumann BC at left SiO2-Si junction 
                A(i,i-1) = -esio2/dz; 
                A(i,i) = (esio2 + esi)/dz; 
                A(i,i+1)=-esi/dz; 
            elseif i == I4 % Neumann BC at the right SiO2-Si 
junction 
                A(i,i-1) = -esi/dz; 
                A(i,i) = (esio2 + esi)/dz; 
                A(i,i+1) = -esio2/dz; 
            else 
                A(i,i-1) = 1/(dz^2); 
                A(i,i) = -2/(dz^2); 
                A(i,i+1) = 1/(dz^2); 
 end 
        end 
  
                Vold = V;      % potential from the previous 
iteration 
                V = inv(A)*B;  % potential without the metal 
parts 
  
                % adding the constant potential on the both 
the sides for metal 
                znew = [-tm-tox:dz:-tox-dz z 
(tsi+tox+dz):dz:(tsi+tox+tm)]; 
  
                I1 = (find(znew==-tm-tox));   % index of 
leftmost point (metal) 
                I2 = (find(znew==-tox));      % index of the 
left Metal-SiO2 junction 
                I3 = (find(znew==0));         % new index of 
the left SiO2-Si junction 
                I4 = (find(znew==tsi));       % new index of 
the right SiO2-Si junction 
                I5 = (find(znew==tsi+tox));   % index of the 
right Metal-SiO2 junction 
                I6 = (find(znew==tsi+tox+tm));% index of 
rightmost point (metal) 
  
                lennew = length(znew); 
                phi = zeros(lennew,1); % final potential 
profile with metal parts on the both sides 
                phi(I1:I2-1,1) = VG;   % adjusting the 
potential for metal on the left side 
                phi(I2:I5,1) = V;      % potential profile 
from Poisson's eqn found previously 
                phi(I5+1:I6,1) = VG;   % adjusting the 
potential for metal on the right side 
         
        % adjusting the band diagram 
        E = -phi+VG;  %eV 
        E(I2:I5,1) = E(I2:I5,1) + 3.25; 
        E(I3:I4,1) = E(I3:I4,1) - 3.1; 
  
        %% 
        if abs(sum(Vold)-sum(V)) < 0.1 % setting threshold 
            break 
        end 
        a = 0.95; 
        V = a*Vold + (1-a)*V; % Gauss-Seidal convergence 
        %% Schrodinger part 
  
        h = 6.626e-34;
        h_ = h/(2*pi); 
        m0 = 9.1e-31; 
        m = 1.08*m0; 
         
        % adjusting the band diagram, this time without the 
metal portions 
        E_wo_metal = -V+VG; 
        % plot(z/nm, E_wo_metal) 
        E_wo_metal = E_wo_metal + 3.25; 
        E_wo_metal(find(z==0):find(z==tsi),1) = 
E_wo_metal(find(z==0):find(z==tsi),1) - 3.1; 
        % plot(z/nm, E_wo_metal) 
         
        % forming the Hamiltonial matrix 
        H0 = diag(((h_^2)/(m*dz^2))+E_wo_metal*q); 
        H1 = diag((ones(1,length(E_wo_metal)-1))*(
(h_^2)/(2*m*dz^2)),1); 
        H2 = diag((ones(1,length(E_wo_metal)-1))*(
(h_^2)/(2*m*dz^2)),-1); 
        %  
        H = H0+H1+H2; 
        [psi, En] = eig(H);  % eigenvalues 
         
        % normalizing psi 
        area=(dz*trapz(psi.^2)).^-1; 
        psi=psi.^2*diag(area); 
        psi=sqrt(psi); 
        % plot(z/nm,abs(psi(:,1))) 
  
        %% C-V 
        nv = 6;    % degeneracy for Si 
        md = 0.358*m0;  % effective mass 
        K = 1.38e-23; 
        T = 300; 
        EF = q*(-1.04+VG);   % Fermi energy 
  
        En1 = En(1:5,1:5);   % working with the first five 
wavefunctions 
        En1 = diag(En1); 
        psi1 = psi(:,1:5); 
        Nij = (nv*md*K*T)/(pi*h_*h_); 
        log_part = log(1+exp((EF-En)./(K*T))); 
        temp = (psi.^2)*diag(log_part); 
        nz(:,n) = Nij*sum(temp,2); 
    end 
    Q(n) = q*abs(trapz(dz,nz(:,n)));    % Charge density 
end 
 
C = diff(Q)./diff(VG_array);   % Capacitance per unit area
VG_array1 = (VG_array(2:end)+VG_array(1:(end-1)))/2; % 
adjusting the length due to use of the diff() 
figure() 
plot(VG_array1,C*1e2,'-r','LineWidth',2);
ylim([0 4.2]);
xlabel('Gate Voltage, V_{G} (Volts)');
ylabel('Gate Capacitance, C_{G} (\mu F/cm^{2})');
set(gca,'FontSize',18);
