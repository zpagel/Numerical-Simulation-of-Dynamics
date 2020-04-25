%% Simulation Parameters and Things.
%E_R is 160kHz for lithium; 8kHz for 133Cs
E_D=4; % Drive energy in recoils, chicago drives at 7kHz on Cs; 
moddepth=0.035;
delta0=2*pi*moddepth; % 5% modulation; the chicage is 65nm/1064nm 
numSteps=101; % Number of time steps
gammaa=1;
c=2E-1;
numBands=5;

hbar=1.0545718E-34; %hbar
amu=1.660468E-27;
m=7*amu;
% m=6*amu;
k_L=2*pi/(1064E-9);
Er=hbar^2*k_L^2/(2*m);

vr=hbar*k_L/m;                 % recoil velocity

fr=Er/(2*pi*hbar);

% Fourier Transform Stuff
numStates=51;%Number of states

stateHigh=10;

kpoints=501; %number of points in FBZ to calculate
kvec=linspace(-1,1,kpoints);

phiStates=2000;
phiBounds=[-pi/2 pi/2];
xData=transpose(linspace(phiBounds(1),phiBounds(2),phiStates));





