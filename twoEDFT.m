function Etot1 = twoEDFT(Z,S,Nr)

% Z is nuclear charge
% S is an appropriate screening constant
% Nr is the number of steps in our linear grid

% This function is a simple density functional theory method of finding the
% ground state energy of a two electron atom/ion, using the local density
% approximation and electron gas density model to account for electron
% exchange energy.

bohr = 0.48*Nr; % Maximal probability of electron density

% To prevent computers from crashing
if Nr>10000
    error('Step size too large')
end
% To prevent non-physical answers
if S>2
    error('Electron shielding cannot be larger than total electron charge')
end
% To prevent poor answers
if rem(bohr,1) ~= 0 ;
    error('Step size poorly set')
end
% Grid and R definition 
Nr3 = Nr^3;
dist = 40;
 
% Define 3D box dimensions for use with radial distance from centre
m=linspace(-dist,dist,Nr);
h = m(2)-m(1);
x = m;
y = m;
z = m;
R=sqrt(x.^2+y.^2+z.^2);
% define x and y and z space
 
% Laplacian
one=ones(Nr,1);
Lap=spdiags([one -2*one one], -1:1, Nr, Nr)/h^2;
%Id = speye(Nr);
%L = kron(kron(L,Id),Id) + kron(kron(Id,L), Id) + kron(kron(Id,Id),L);
 
% Potentials
% DFT states that energy of a system can be represented as a functional of 
% electron density; which is in turn a function of R (distance from
% 'nucleus')
% Energy(p(r)) = T(p(r)) + Hh(p(r)) + Eext(p(r)) + Exc(p(r))
% We know all apart from Exc(p(r))
 
% Constants
hbar = 4.135667662 * 10^(-15); % eV.s
m = 0.511 * 10^6; % eV/c2
 
% Only 'known' part of energy, external (central) potential
Vext = -(Z-S)./R;
% Eext = n.*Vext;
% Coulombic interactions of electrons with nuclei, assuming 
% Born-Oppenheimer approximation (nuclei are 'fixed')
 
% Initial total potential definition
Vtot = Vext;
 
b = (-hbar^2/2*m)*Lap+diag(Vtot);
 
% Defining psi (KS orbital) (found using Hartree-Fock methods)
[psi, E] = eigs((-hbar^2/2*m)*Lap+diag(Vtot), 1, 'sa');
psi = Nr^0.5 * psi;

% Loop for generating K-S wavefunctions
j = 1;
while j<=5
% p(r) = sum(psi(i))^2 , where psi is Kohn-Sham orbital
n = Z*psi.^2; % Electron Density
 
% We use the local density approximation, and approximate the system to a
% homogeneous electron gas; the equation for this is the form of Vxc below
 
% Iterative loop until self-fulfilling equations
 
T = Z*psi'*(-0.5*(hbar^2/2*m)*Lap)*psi*h;  % Kinetic Energy Term
 
% External potential at each electron density
Eext = sum(n.*Vext*h);
 
% Solving Poisson equation to get electronic repulsion
Vh = zeros(length(x));
for x_1 = 1:Nr 
    for x_2 = 1:Nr
        if x_2 ~= x_1 % avoiding infinite hartree energies
            Vh(x_1, x_2) = n(x_1)*n(x_2)/(abs((x(x_1) - x(x_2))));
        end
    end
end
Eh = 0.5 * sum(n.*Vh*h^2); % Hartree Energy
% (electron-electron Coulombic interaction, neglects both exchange and
% correlation) (0.5 prevents double counting)
          
Vxc = -(3/pi)^(1/3)*n.^(1/3); % Electron Exchange and Correlation 
                              % Contributions
Exc = (3/4)*Vxc.*n*h;
 
% Total potential
Vtot = Vxc + Vh + Vext;
% Total Energy Therefore
Etot = T + Eext + Eh + Exc;
j = j+1;
end
Etot1 = Etot(1,bohr);
disp('Total Potential Energy is a combination of Kinetic, External, Hartree and Exchange Correlation Energies')
disp('in eV')

end