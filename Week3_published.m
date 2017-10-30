%% Phase 3 - Final steps with Schrodinger Equation, and then moving on to Kohn Sham Equations

%% Particle in a spherical box and central potential
% Before moving on to other methods for solving the Schrodinger Equation we
% decided to try a similar method as for particle in a box but using
% spherical coordinates.  Doing this we could attempt to get reasonable
% eigenfunctions for the radial and angular components.  The idea was that
% once this was working it should be straightforward to introduce a
% central potential and get the spherical harmonics.  A major issue we had
% when solving these eigenstates is that the MATLAB eig() function places
% the calculated eigenfunctions in seemingly random order.  This means that
% for this example we had to search through all 200 eigenfunctions to be
% able to identify candidates for the ground state and excited states.
%
% Below is the code to find the radial functions for the particle in a hard
% sphere and followed by that for a particle in a hard sphere witha
% central potential.  
 
N = 200; R = 1; r = linspace(0,R,N); dr = r(3)-r(2);
Dr = (1/(2*dr))*(diag(ones((N-1),1),1) - diag(ones((N-1),1),-1));
 
Vr = zeros(N,N);
Vr(N,N) = 100000;
% potential defined by zero everywhere except at the limit of the sphere
% where it is large.
 
r(1) = 0.000001;   % to avoid dividing by zero in the laplacian
lapr = diag(1./r.^2)*Dr*(diag(r.^2)*Dr);
H = ones(N,N,3); H(:,:,1) = (-1/2)*lapr + Vr;
[Efns_r, ~] = eig(H(:,:,1));
 
Efns_r_GS = Efns_r(1:2:end,95);
Efns_r_GS = Efns_r_GS./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_GS .* conj(Efns_r_GS)));
% possible ground state given by the 95th eigenfunction, with
% normalisation. 
 
Efns_r_1st = Efns_r(1:2:end,97);
Efns_r_1st = Efns_r_1st./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_1st .* conj(Efns_r_1st)));
% possible 1st excited state given by the 97th eigenfunction, with
% normalisation. 
 
Efns_r_2nd = Efns_r(1:2:end,116);
Efns_r_2nd = Efns_r_2nd./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_2nd .* conj(Efns_r_2nd)));
% possible 2nd excited state given by the 116th eigenfunction, with
% normalisation. 
 
hold on 
plot(r(1:2:end), Efns_r_GS)
plot(r(1:2:end), Efns_r_1st)
plot(r(1:2:end), Efns_r_2nd)
title('Particle in a sphere - radial eigenfunctions')
legend('Possible 1st excited state', 'Possible 2nd excited state', 'Possible 3rd excited state','location','southeast')
xlabel('r')
ylabel('\psi (r)')
grid on
hold off
 
%%
% Obtaining the angular wavefunctions:
N = 500; R = 1;
r = linspace(0,R,N); t = linspace(0,pi,N); a = linspace(0,2*pi,N);
dr = r(3)-r(2); dt = t(3)-t(2); da = a(3)-a(2);
Dt = (1/(2*dt))*(diag(ones((N-1),1),1) - diag(ones((N-1),1),-1));
 
r(1) = 0.000001; t(1) = 0.000001; a(1) = 0.000001;
lapt = diag(1./r.^2)*diag(1./sin(a).^2)*Dt*Dt;
H = ones(N,N,3);
H(:,:,2) = (-1/2)*lapt;
 
[Efns_t, ~] = eig(H(:,:,2));
Efns_t = Efns_t./sqrt(trapz(t,(Efns_t.*conj(Efns_t))));
 
for i=1:4
    subplot(2,2,i)
    plot(t,Efns_t(:,i+223).*conj(Efns_t(:,i+223)))
    xlabel('\theta')
    ylabel('\psi (\theta)')
end
suptitle('A selection of Angular Eigenfunctions')
 
%%
% As you can see these angular wavefunctions look well behaved.
% These are eigenfunctions 220-224.  An issue arises in that they are not
% symmetric, as one would expect for a periodic wavefunction.  In order to
% impose periodicity one could take this further by defining the values of
% the laplacian at the start and end point to be the same. The azimuthal
% eigenfunctions as of yet have not given us any well-behaved functions.
% Further work could go into obtaining these.  
% Now introducting a central potential,
 
N = 200; R = 10; r = linspace(0,R,N); dr = r(3)-r(2); 
Dr = (1/(2*dr))*(diag(ones((N-1),1),1) - diag(ones((N-1),1),-1));
r(1) = 0.0000000001; 
Vr = diag(-1./r);
% Potential defined here
 
lapr = diag(1./r.^2)*Dr*(diag(r.^2)*Dr);
H = ones(N,N,3);
H(:,:,1) = (-1/2)*lapr + Vr;
 
[Efns_r, E_r] = eig(H(:,:,1));
 
Efns_r_GS = Efns_r(2:2:end,119);
Efns_r_GS = Efns_r_GS./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_GS .* conj(Efns_r_GS)));
% possible ground state given by the 119th eigenfunction, with
% normalisation.  
 
Efns_r_1st = Efns_r(1:2:end,120);
Efns_r_1st = Efns_r_1st./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_1st .* conj(Efns_r_1st)));
% possible 1st state given by the 120th eigenfunction, with
% normalisation.  
 
Efns_r_2nd = Efns_r(2:2:end,133);
Efns_r_2nd = Efns_r_2nd./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_2nd .* conj(Efns_r_2nd)));
% possible 2nd state given by the 133th eigenfunction, with
% normalisation. 
 
Efns_r_3rd = Efns_r(2:2:end,134);
Efns_r_3rd = Efns_r_2nd./sqrt(trapz(r(1:2:end), 2*pi* diag(r(1:2:end))*...
    Efns_r_3rd .* conj(Efns_r_3rd)));
% possible 3rd state given by the 134th eigenfunction, with
% normalisation. 
 
figure
hold on
plot(r(1:2:end),Efns_r_GS)
plot(r(1:2:end),Efns_r_1st)
plot(r(2:2:end),Efns_r_2nd)
plot(r(2:2:end),Efns_r_3rd)
grid on 
title('Particle with a central potential - radial eigenfunctions')
legend('possible ground state','possible 1st state', 'possible 2nd state', 'possible 3rd state')
xlabel('r')
ylabel('\psi (x)')
hold off
 
%%
% It is important to note that due to the odd behaviour of the eig()
% function in that it gives vectors containing both the positive and
% negative eigenfunctions aswell as lots of zero values and also that the
% well-behaved eigenfunctions are hidden among them all, makes it difficult
% to write user-defined functions to do these calculations.  For example
% sometimes one has to extract every fourth element from the eigenfunction,
% or every 4th.  
 
%% Two interacting electrons in a 1D infinite square well - an attempt at DFT
% We then decided to attempt a simple application of Density functional
% theory - finding the eigenfunctions of two interacting electrons in an
% infinite square well. The first step in a DFT problem is to solve the
% single electron Schrodinger equation for each electron, assuming that
% they do not interact at all. This part of the code is simply solving the
% 1D infite potential well Schrodinger equation:
  
% Set parameters, constants + 'infinite' potential well of length L.
clc
clear
N = 200; L = 5;                                      
x = linspace(0,L,N)'; dx = x(2)-x(1);
  
% Effectively infinite potential well:
V = zeros(N); V(1) = 100000; V(N) = 100000; 
   
% Define Laplacian and Hamiltonian
Lap = (1/(dx^2))*(-2*diag(ones(N,1),0) + diag(ones((N-1),1),1) + diag(ones((N-1),1),-1));
H = (-1/2)*Lap + diag(V);
   
%Solve one-electron Schrodinger equation:
[psi_1, E1] = eigs(H);

%Removing the first, trivial eigenfunction extracting every other value,
%so that the functions are single-valued:

psi_1 = psi_1(:,2:size(psi_1,2));
psi_1 = psi_1(1:2:end, :);
  
% Re-sizing x, dx, Lap and potential:
x = x(1:2:end); dx = x(2)-x(1);
Lap = (1/(dx^2))*(-2*diag(ones(length(x),1),0) + diag(ones((length(x)-1),1),1) + diag(ones((length(x)-1),1),-1));
V = zeros(length(x)); V(1) = 1000000 ; V(length(x)) = 1000000;

  
%% Iteratively Solving the Kohn-Sham Equations for two Electrons in an Infinite Square Well
% Up to this point, all we have done is standard quantum mechanics, with no
% DFT aspects. However, if we wish to introduce a second particle, which
% interacts with the first, this results in a problem which cannot be
% solved by analytic quantum mechanics. The previous section gave us a 
% single electron wavefunction which can be used as a starting point for
% multiple-electron DFT.
  
%%
% Density funtional theory deals with this problem by considering single
% electron wavefunctions and how they are affected by one another. The
% steps for solving a problem using DFT are as follows:
%
% # Find a trial electron density (from single-electron wavefunctions).
% # Find the total energy for that density.
% # Solve the Kohn-Sham equation to find new single electron wavefunctions
% # Use these wavefunctions to find a a new electron density.
% # Repeat steps 2-4 until the total energy stops decreasing. The energy
% will keep decreasing with every iteration, but it will tend towards the
% 'true' energy, so the process only needs to be repeated until the energy
% stops changing up to a desired accuracy.
%
% Finding the total energy requires the following components:
%
% * The kinetic energy of the electrons. Found by taking the expectation
% value of the kinetic energy operator.
% * The Hartree energy. This is the electron-electron Coulomb interaction.
% It is found by double integrating the product of both densities.
% * The exchange-correlation energy. In this example, we use the local
% density approximation (LDA) for the exchange energy. The LDA uses an
% electron gas model find analytic expressions for this term. In this
% example, we ignore the correlation energy as there is no analytic
% expression for it - it must be approximated using quantum Monte-Carlo 
% methods It is not the biggest contribution to energy, so we hope it will 
% not affect our results too much.
% * The interaction energy between the electrons and the external
% potential. In our case, this is zero inside the well.
%%
% The Kohn-Sham equations requires the kinetic energy operator, the
% exchange correlation potential (found from the LDA), the Coulomb potential
% due to the change density distribution and the external potential (which
% in this case is zero inside the well).
%%
% In this example we will just consider the ground state:
  
% Extracting the ground state:
PSI = psi_1(:, 1);

% Normalising:
psiarea = sum((PSI.^2*dx));
PSI = (psiarea^(-1/2))*PSI; 
  
% LOOP THIS PART UNTIL ENERGY DOES NOT CHANGE (requires about 10 iterations)

for i = 1:10
    % Defining electron charge and density function, since we are working
    % in atomic units, the charge of the electron is taken to be 1.
    den = 2*PSI.^2 ; 
    
    % Plotting the electron density
    subplot(2,5,i)
    plot(x, den)
    set(gca,'fontsize',7)
      
    % Using the LDA to get the exchange potential functional (ignoring the
    % correlation functional):
    V_ex = -(3/pi)^(1/3)*sum((den.^(1/3))*dx);
   
    % Kinetic energy (expectation value of KE operator (x2))
    T = 2*PSI'*(-1/2)*Lap*PSI*dx ; 
   
    % Hartree energy
    % Creating a matrix (roughly equivalent to hartree potential) to sum the 
    % elements of in order to perform the required double integration:
   
    Vh = zeros(length(x));
    for x_1 = 1:length(den) % 
        for x_2 = 1:length(den)
            if x_2 ~= x_1 % avoiding infinite hartree energies
                Vh(x_1, x_2) = den(x_1)*den(x_2)/(abs((x(x_1) - x(x_2))));
            end
        end
    end
   
    Eh = -(1/2)*sum(sum(Vh*dx^2)); % Hartree Energy
   
    % Exchange energy (just integral of exchange potential):
    E_ex = -(3/4)*(3/pi)^(1/3)*sum(den.^(4/3).*dx);
  
    % Finding the total energy and adding to graph title:
    Etot = Eh + T + E_ex;
    title({['Iteration ', num2str(i),] ,['Energy = ', num2str(round(Etot, 1)), 'eV']}, 'Fontsize', 7)
      
    % Create Coulomb potential function:
    Vc = zeros(length(x));
    for x_1 = 1:length(x)
        for x_2 = 1:length(x)
            if x_1 ~= x_2
                Vc(x_1, x_2) = den(x_2)/(abs((x(x_1) - x(x_2))));
            end
        end
    end
   
    % Integrate Vc over x_2:
    Vc2 = sum(Vc*dx,2);
   
    % Finding the eigenfunctions from the Kohn-Sham equations:
    [PSI, E] = eig((-1/2)*Lap + diag(Vc2) + V_ex*eye(length(x)) );
  
    % Extracting the ground state only:
    PSI = PSI(:,1);
    
    % Normalising:
    psiarea = sum((PSI.^2*dx));
    PSI = (psiarea^(-1/2))*PSI; 
end
suptitle('')
  
%%
% Here, the y-axis represents electron density in all graphs.
% This method produces an energy which reaches a steady value after several
% iterations, as expected. However, it does not always decrease with each
% iteration, as should happen - the energy fluctuates up and down 
% occasionally before settling on a final energy. 
% Also, the electron density curve does not reach an unchanging
% distribution. It converges to a curve, but then alternates between having
% a maximum on the left side and the right side of the box, depending on
% the iteration.
% This is not an ideal 
% outcome and this result has probably arisen due to the fact that we only
% used one pseudo-wavefunction to describe both electrons, rather than two
% separate single electron wavefunctions. However, we are happy that the
% main DFT components of the code are functioning well enough to move on to 
% more complicated problems.
 
%% Helium Atom, and Two Electron Ions
% Now that we have managed to successfully code in a central potential, as 
% well as the electron-electron interaction energies associated with a 1D
% two electron infinite potential well, we may combine the two to create
% some simple density functional code for a two electron atom (such as
% Helium).
%
% This provides a good test of our code for us, as this is the first
% physically interesting system that we cannot solve analytically within
% the Schrodinger equation, requring us to use the approximations we have
% developed above.
% 
% Also, for ease of use and to demonstrate technical skill, we have placed
% the code inside a function: twoEDFT.m , which takes three inputs. Z, the
% nucleus charge; S, a static screening constant for one electron relative
% to the other; and Nr, the number of steps we wish to take in our
% calculations. 
% It has been necessary to introduce some error messages in the function
% both to prevent the average computer from crashing, and to prevent a user
% attempting to obtain non-physical values (incorrect use of the function).
%
%%
% As previously done with the two electron in a box model, we are using the
% local density approximation to model electron exchange, with the electron
% uniform gas model.

%%
% We may input the standard parameters for Helium at 1000 steps, obtaining 
% (in eV) a total energy for both electrons in the Helium system.
% Ground state ionisation energy:

A = twoEDFT(2,0.3,1000);
disp(A)


%%
% This is relatively close to an experimental value of -108.8 eV 
% (The PubChem Project) that we wish to obtain.
% Given that we are neglecting electron correlation, not using 
% error-correction methods such as compensation charges
% and artificially setting boundary conditions at
% the 'edge' of the atom to zero, this seems like quite a reasonable
% output. We consider the endeavour relatively successful.
%
% Of course, if we had more time, we would choose better (but more complex)
% approximations, perform a Monte Carlo simulation in order to better 
% account for electron correlation, and compare our results in more complex
% scenarios. 
%
% However, it is worth noting that DFT code is still a currently active
% area of research, and many scientists spend much time and effort on
% finding the best way to programme their quantum mechanics. As a student
% project, we believe that this is a very large step onto the first
% stepping stone of understanding and using one of atomic and molecular
% physics/ chemistry's most useful tools for research.
%
%% References
% The PubChem Project. USA: National Center for Biotechnology
% Information _"Helium - PubChem Public Chemical Database"._ Retrieved from
% https://pubchem.ncbi.nlm.nih.gov/compound/23987
% 