%% Spin Evolution Function
% Ohio Advanced EPR Laboratory
% Robert McCarrick
% Vector Simulation Package, Version 2.0
% Last Modified 1/23/2012

function [Spinsout Sim]=evolutionaudio(Spins,Sim)

%% Checks for variables and falls back on a default value if missing

if (isfield(Spins,'degrees') == 0);
    error('No spin array (.degrees) specified in structure.');
end

if (isfield(Spins,'offset') == 0);
    error('No spin offset (.offset) specified in structure.');
end

if (isfield(Spins,'Tone') == 0);
    Spins.Tone=100;
end

if (isfield(Spins,'Ttwo') == 0);
    Spins.Ttwo=100;
end

if (isfield(Sim,'length') == 0);
    error('No length (.length) specified in structure.');
end

if (isfield(Sim,'framerate') == 0);
    Sim.framerate=20;
end

if (isfield(Sim,'video') == 0);
    Sim.video='n';
end


if (isfield(Sim,'nframes') == 0);
    Sim.nframes=100;
end

if (isfield(Sim,'figuresize') == 0);
    Sim.figuresize= [900 500];
end


%% Gets the number of spins in the array

mdegrees = size(Spins.degrees);
starts = zeros(mdegrees(1),3);

%% Converts spins to rectangular coordinates

coordinates = zeros(mdegrees(1),3);

for i=1:mdegrees(1)
    x = sin(Spins.degrees(i,2)/360*2*pi)*cos(Spins.degrees(i,1)/360*2*pi)*Spins.degrees(i,3);
    y = sin(Spins.degrees(i,2)/360*2*pi)*sin(Spins.degrees(i,1)/360*2*pi)*Spins.degrees(i,3);
    z = cos(Spins.degrees(i,2)/360*2*pi)*Spins.degrees(i,3);
    coordinates(i,:)=[x,y,z];
end

%% Creates time axes for the plots

Sim.inc = Sim.length/Sim.nframes;
time = Sim.inc:Sim.inc:Sim.length;
xsum = zeros(1,numel(time));
ysum = zeros(1,numel(time));
zsum = zeros(1,numel(time));

%% Calculates the xy and z projections for T1 and T2 decays

xyprojection = zeros(1,mdegrees(1));
zprojection = zeros(1,mdegrees(1));

for m=1:mdegrees(1)
    xyprojection(m) = (coordinates(m,1)^2+coordinates(m,2)^2)^0.5;
    zprojection(m) = coordinates(m,3);
end


%% Performs the incremental rotation and plots each frame with no animation output

relaxation = zeros(numel(mdegrees(1)),3);
ends  = zeros(numel(mdegrees(1)),3);

for i=1:numel(time)
    for j=1:mdegrees(1)
        theta=2*pi*Spins.offset(j)*Sim.inc;
        rotmat=[cos(-theta) -sin(-theta) 0;sin(-theta) cos(-theta) 0;0 0 1];
        scale=xyprojection(j).*(2.71828183^(-time(i)./Spins.Ttwo));
        coordinates(j,:) = coordinates(j,:)*rotmat;
        relaxation(j,1)=coordinates(j,1).*scale/xyprojection(j);
        relaxation(j,2)=coordinates(j,2).*scale/xyprojection(j);
        relaxation(j,3)=1+(zprojection(j)-1).*(2.71828183^(-time(i)./Spins.Tone));
        ends(j,:)=relaxation(j,:);
    end
    xsum(1,i)=round(sum(ends(:,1))/mdegrees(1)*100)/100;
	ysum(1,i)=round(sum(ends(:,2))/mdegrees(1)*100)/100;
    zsum(1,i)=round(sum(ends(:,3))/mdegrees(1)*100)/100;

end

%% Audio Generation
max(time)
timeaudio = linspace(0,max(time)*100,max(time)*100*8000);
yinterp = interp1(time*100,xsum,timeaudio);
yaudio = sin(timeaudio*2*pi*500).*yinterp;

subplot(3,1,1)
plot(time,xsum)
subplot(3,1,2)
plot(time,ysum)
subplot(3,1,3)
plot(timeaudio,yinterp)

sound(yaudio,8000)

% Converts the rectangular coordinates back to spherical coordinates

r = zeros(1,mdegrees(1));

for k=1:mdegrees(1)
    r(k)=((ends(k,1))^2+(ends(k,2))^2+(ends(k,3))^2)^0.5;
end

angles = zeros(mdegrees(1),2);

for l=1:mdegrees(1)
    angles(l,:)= [atan2(ends(l,2),ends(l,1))*180/pi acos(ends(l,3)/r(l))*180/pi];
end

Spinsout = [angles r'];