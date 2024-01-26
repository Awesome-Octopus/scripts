%%  This simulation depicts a typical 1H NMR experiment on ethanol

clear all

%% Creates a matrix of spins along the z axis

Spins.degrees = [0,0,0.7;0,0,0.25;0,0,0.5;0,0,0.5;0,0,0.25;0,0,0.5;0,0,1;0,0,0.5;];
Spins.offset = [4.7 3.59 3.57 3.55 3.53 1.12 1.11 1.10].*1e-6*200000000;
Spins.Tone = 1;
Spins.Ttwo = 0.05;

%% Sets the parameters and function for first pi/2 pulse

Sim.direction = 'y';
Sim.angle = 90;
Sim.framerate = 20;
Spins.degrees = pulsespin(Spins,Sim);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.1;
Sim.framerate = 20;
Sim.nframes = 1200;

Spins.degrees = evolution(Spins,Sim);