%% Phase 2 - Taking Schrodinger Equation Further
% Once we were confortable using the Schrodinger Equation for relatively
% simple potentials such as particle in a 2D box, we took it further by
% first looking at systems which need new graphical methods, and then
% systems that need polar coordinates. This allows us to develope
% understanding 2D and 3D polar coordinates within MATLAB which may be
% useful when considering separable wavefunctions for atomic density
% functional theory.
%
% We also look at preparing our code for extensions and more complex work
% by establishing a more adaptable non-uniform spacial grid density for
% saving computational time.
%% Particle in a 3D box
% The particle in a 3D box requires an extension to the 2D particle
% in a box.  The complications arise when creating the 3-dimensional
% wavefunction from the 3 1-dimensional wavefunctions.  Whereas for a 2D
% system you can just multiply them together to create an n-n array, for 3D
% it is not that easy.  To do it we created a triple loop which goes
% through every element in an n-n-n array with indices i,j and k.  It then
% generates that element by multiplying the ith element of the
% x-wavefunction with the jth elements of the y-wavefunction and the kth
% element of the z-wavefunction.   This satisfies the seperable assumption
% for this system,  where 
%
% $$  \Psi (x,y,z) = \psi(x) \, \psi(y) \, \psi(z)$$
% 
% The second complication was presenting the results.  There are various
% options, such as surface plots or contour slice plots.  We decided on
% balance the first was more appropriate and easy to read.
%% 
% As with the 2D particle in a box, with an added dimension.  To do this
% efficiently I wrote a user-defined-function called cart_eig() which 
% calculates the eigenfunctions for each individual dimension.  

clc
clear
N = 100;
V = zeros(N,N);  V(1,1) = 1000000;  V(N,N) = 1000000;

[ x,Efns_x ] = cart_eig( 'x', N, 0.5, V);
[ y,Efns_y ] = cart_eig( 'y', N, 0.5, V);
[ z,Efns_z ] = cart_eig( 'z', N, 0.5, V);

hold on
subplot(1,2,1)
[ ~ ] = isosurf_eig_plot( x,y,z,[1,1,1], Efns_x, Efns_y, Efns_z, 2);
title('Ground state')
subplot(1,2,2)
[ ~ ] = isosurf_eig_plot( x,y,z,[2,1,1], Efns_x, Efns_y, Efns_z, 2);
title('1 excitation in the x dimension')
hold off
%% 
% these user-defined-functions create the surface data to plot the 
% isoplanes of the wavefunction. The blue surfaces represent where the wavefunction is negative, 
% and the red surface is where it is positive.  
  
%% Defining a Periodically Varying Grid Density
% Thus far in developing our code for 1D, 2D and 3D quantum wavefunctions,
% we have been spacing our pre-defined x, y, and z spacial coordinates as a
% linearly spaced list of increasing coordinates between two defined
% 'start' and 'end' values, with 'N' numbers of linearly spaced steps.
%
% Whilst this makes sense for small systems and works well, if we wish to compute
% accurate wavefunctions and obtain well-defined potential energy
% curves this would require a very large density of points around
% areas of interest such as peaks.
%
% For a linear system this requires us to input a very high step count,
% which can cause our programme to slow down unnecessarily. We hypothesise
% that it might be less computationally stressful to develop a non-uniform
% grid density of points. This can be used when a linear system point 
% system is taking to long to calculate.

%%
% We start by trying to develop an exponentially decaying grid density
% using step size given by 
%
% $$ h(i) = x(i) - x(i-1) $$
%
% This would be useful for spherical systems in the ground
% state as it will give high resolution around the origin.  
%
% Basic code for exponential decay from origin, step size
%

d1 = 0.05; rmax1 = 5;rmin1 = 1; nmax1 = 1.5;
% Parameters define factor of exponential decay, and range of points

i = 1;  x1 = 1:1:100;
while i<100  
i = i+1;
x1(i) = rmax1 + (rmax1 - rmin1)*(exp(i*d1)-1)/(exp(nmax1*d1)-1);
end
% This loop alters the point-to-point distance in our non-uniform grid.
%%
% The next step is to patch together point distributions as above but
% switch the exponential round from 'decay' to 'gain', and then keep
% reversing it periodically as many times as is defined. This creates areas
% of periodically varying grid density. These may be tuned to matchup with
% suspected areas of interest. This point distribution would be a natural
% choice if doing calculations with periodic boundary conditions,  such as
% with lattices.
clf
d2 = 0.05; rmax2 = 5; rmin2 = 1; nmax2 = 1.5;

k2 = 2; j2 = 2; y2 = zeros(1,1000); i = 1; 
% Initialising; k determines number of periods, j is the element number for
% x, y tells us the difference between x values, i varies 1 to 100. 

x2 = 1:1:100; % x is step sizes

while k2<=10 % k/2 is the number of peaks/troughs
    while 1<=i && i<=100 % capping i 
x2(j2) = (rmax2 + (rmax2 - rmin2)*(exp(i*d2)-1)/(exp(nmax2*d2)-1));
y2(j2+1) = y2(j2)+x2(j2);
i = i+(-1)^k2; % i varies 1 to 100, then k changes, and we go 100 to 1
j2 = j2+1; 
    end
    k2 = k2+1;
    i = i+(-1)^k2;
end

sz = size(y2,2);

% A simple plot to illustrate the periodic and exponential varying grid
% densities.  Some scaling factors have been added to allow them to be on
% the same axis.  
test1 = ones(1,100);
test2 = ones(1,sz);
hold on
plot(y2,test2,'.', 'markers', 3);
plot(x1.*188,test1+1,'x', 'markers', 3);
ylim([0.4 2.6])
title('Graph showing two methods of spacing intervals')
hold off
%%
% As we can see, this can allow for much more enhanced precision at certain
% points of interest without needing the entire spacial grid to be very
% dense.
%
% In order to implement this it would be done in an analytical fashion: to
% begin with we run the programme with a low density linear grid,
% establish points of interest, then run with a higher density non-uniform
% grid in order to better investigate these points of interest. However
% this is not straight forward and so we suspect that this only becomes
% useful when dealing with more complex structures.
%
% Therefore we may not use this in our future applications within the
% scope of this MATLAB project. However it would certainly be
% useful for more complex extensions of this project.

%% Particle in a Corral - Solving the Schrodinger Equation
% We decided that if we wanted to progress to solving problems involving
% atoms we would need to work in polar coordinates. As an intermediate step
% we decided to solve the Schrodinger equation in a 2D circular infinite
% potential (also known as a 'corral') as well using the polar Laplacian.
% Since the polar Laplacian is not defined at r=0, we started the
% coordinate from r=0.01. We then solved the Schodinger equation using
% separation of variables in the same way as before.
 
clc
clear
  
N = 200; R = 2; theta = linspace(0,2*pi, N); r = linspace(0.01,R, N);            
dr = r(2) - r(1);  dtheta = theta(2) - theta(1);                                
  
% first derivative with respect to r and theta:
Dr = (1/(2*dr))*(diag(ones((N-1),1),1) - diag(ones((N-1),1),-1));
Dtheta = (1/(2*dtheta))*(diag(ones((N-1),1),1) - diag(ones((N-1),1),-1));
  
% Diagonalise r vector and (1/r) vector so it can multiply matrices:
rmat = diag(r);         % r
rmat_inv = diag(r.^-1); %(1/r)
  
% r component of the 2D polar laplacian:
lapr = rmat_inv*(Dr*(rmat*Dr));
  
% theta component of Laplacian
laptheta = (rmat_inv.^2)*(Dtheta*Dtheta);
 
% Setting boundary conditions:
laptheta(1,1) = 0; laptheta(1,2) = 0; laptheta(2,1) = 0;
laptheta(N,N-1) = 0; laptheta(N-1,N) = 0; laptheta(N,N) = 0;
  
% Hamiltonian inside the well (V = 0):
% Define the Hamiltonian as an NxNx2 multidimensional array
H = ones(N, N, 2);
H(:,:,1) = (-1/2)*laptheta;
H(:,:,2) = (-1/2)*lapr;
  
[Efns_theta, E_theta] = eig(H(:,:,1));
[Efns_r, E_r] = eig(H(:,:,2));
  
%%
% We then produced radial plots of the eigenfunctions.
% From these radial plots we could see that the eigenvector contains two
% eigenfunctions, different by a sign, as well as a large number of points
% with a zero value. It is evident that the 'eig' Matlab function is
% finding merged eigenfunctions rather than 'pure' ones. After
% trial and error we found that extracting every fourth value in the
% eigenfunction give a 'pure' wavefunction. Extracting every 4th element+2
% would give the same function but negative.

%% Particle in a Corral - the angular wavefunction
% Having obtained the radial wavefunctions we moved on to the angular
% eigenfunctions. Plotting them gave some odd-looking graphs.

Efns_r2 = Efns_r(1:4:end, :);
% extracting every 4th element to give pure eigenfunctions. 

for i = 1:9
    subplot(3,3,i)
    plot(theta, Efns_theta(:,20*i), '.')
    xlim([0,2*pi])
    xlabel('\theta')
    ylabel('\psi ( \theta)')
end
 suptitle('Angular wavefunctions - Test Plot')
 
%%
% These eigenfunctions do not look like superpositions of multiple
% functions, unlike the radial ones we found. We tried extracting  in the 
% same way as for the radial functions and though some of the functions 
% looked a little bit like sine waves we did not find anything that was
% regular and continuous, like an eigenfunction should.
 
%% 
% The lack of good angular eigenfunctions makes plotting the wavefunction
% difficult. Solving the equation analytically gives an angular 
% wavefunction of 
% 
% $$ \frac{1}{\sqrt{2 \pi}} \cdot e^{i \cdot m \cdot \theta} $$
%
% where m = 1,2,3 etc. The
% imaginary componenents may explain why matlab had trouble extracting
% eigenfunctions.
 
%% Particle in a Corral - the Probability Density
% To overcome this problem, we decided to plot the probability density
% instead, as it will not be affected by the lack of angular wavefunctions,
% since |exp(i*m*theta)|^2 = 1. Therefore the angular component will only
% contribute a factor of $1/\sqrt{2pi}$ to the total probability density. 

Efns_r2 = Efns_r(1:4:end, :);  
% extracting every 4th element to give pure eigenfunctions. 

pdensity = (1/2*pi)*Efns_r2.^2;
% Probability density function matrix of eigenfunctions squared
%% 
% Then, we plotted the probability density as a function of position. To
% start with, we only plotted the probability density for the first
% eigenfunction. Since the probability density does not change with theta, 
% we were able to define it as a matrix which changed in the r direction,
% but not in the theta direction.
 
% Define new theta and r column vectors consistent dimensions with pdensity 
theta = linspace(0,2*pi, size(pdensity, 1))';
r = linspace(0.01, R, size(pdensity, 1))';
 
% Using meshgrid to combine them in order to produce a 3D plot.
[tgrid, rgrid] = meshgrid(theta,  r);
 
% Extracting the first eigenfunction
firstefn = pdensity(:,1);
 
% Defining the matrix for z-coordinate (the probability density)
zpolar = repmat(firstefn,1,size(pdensity, 1));
 
%% 
% Matlab does not allow 3D plots in polar coordinates, we had to make use
% of the inbuilt function pol2cart before plotting the resulting surface:
% We then plotted the probability densities for some of the excited states.
  
for i = 1:2
    subplot(1,2,i)
    [X,Y,Z] = pol2cart(tgrid, rgrid, repmat(pdensity(:,i),1,size(pdensity, 1)));
    surf(X,Y,Z)
    xlabel('x')
    ylabel('y')
    zlabel('\psi (x,y)')
    title(['Excited State ', num2str(i-1)] )
end
 
suptitle('Probability Densities for Particle in a Corral')