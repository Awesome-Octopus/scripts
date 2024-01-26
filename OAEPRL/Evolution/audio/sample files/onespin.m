%%  This simulation depicts a typical 1H NMR experiment on ethanol

clear all

%% Creates a matrix of spins along the z axis

Spins.degrees = [0,90,1];
Spins.offset = [100];
Spins.Tone = 0.01;
Spins.Ttwo = 0.05;

%% Sets the evolution parameters and function between pulses

Sim.length = 0.1;
Sim.framerate = 20;
Sim.nframes = 1200;

Spins.degrees = evolutionaudio(Spins,Sim);