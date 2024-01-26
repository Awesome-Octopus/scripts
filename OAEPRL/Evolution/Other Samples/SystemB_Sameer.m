%%  This simulation depicts a typical 1H NMR experiment on ethanol

clear all

Sim.video = 'avi'; % change this to 'gif' to get an animated GIF file
Sim.filename = 'systemB.avi'; % Change this to '900500.gif' if 'gif' is specified above
Sim.figuresize = [900 500];
Sim.numloops = 0;
Sim.framerate = 30;


%% Creates a matrix of spins along the z axis

Spins.degrees = [0,0,0.5;0,0,1];
Spins.offset = [3;0.75].*1e-6*200000000;
Spins.Tone = 1;
Spins.Ttwo = 0.02;

%% Sets the parameters and function for first pi/2 pulse

Sim.direction = 'y';
Sim.angle = 90;
Sim.nframes = 90;
[Spins.degrees Sim] = pulse(Spins,Sim);

%% Sets the evolution parameters and function between pulses

Sim.length = 0.05;
Sim.nframes = 600;

[Spins.degrees Sim] = evolution(Spins,Sim);