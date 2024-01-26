%% Pulse Function
% Ohio Advanced EPR Laboratory
% Robert McCarrick
% Vector Simulation Package, Version 2.01
% Last Modified 10/31/2017

function [Spinsout Sim]=pulsespin(Spins,Sim)

%% Checks for variables and falls back on a default value if missing

if (isfield(Spins,'degrees') == 0);
    error('No spin array (.degrees) specified in structure.');
end

if (isfield(Sim,'angle') == 0);
    error('No angle (.angle) specified in structure.');
end

if (isfield(Sim,'direction') == 0);
    Sim.direction='-y';
end

if (isfield(Sim,'framerate') == 0);
    Sim.framerate=50;
end

if (isfield(Sim,'nframes') == 0);
    Sim.nframes=round(Sim.angle/2);
end

if (isfield(Sim,'video') == 0);
    Sim.video='n';
end

if strcmp(Sim.video,'gif') | strcmp(Sim.video,'avi');
    if (isfield(Sim,'filename') == 0);
        error('Please specify filename for animation.');
    end
end

if (isfield(Sim,'numloops') == 0);
    Sim.numloops=inf;
end

if (isfield(Sim,'figuresize') == 0);
    Sim.figuresize= [900 500];
end

%%  Sets the figure size and position

figurehandle = get(0,'CurrentFigure');
screensize = get(0,'ScreenSize');
currentposition = get(figurehandle,'Position');
desiredposition = [50, screensize(4)-Sim.figuresize(2)-150, Sim.figuresize(1), Sim.figuresize(2)];

    if isempty(figurehandle) == 1
        figure('Position',[50, screensize(4)-Sim.figuresize(2)-150, Sim.figuresize(1), Sim.figuresize(2)]);
    elseif isequal(desiredposition,currentposition) == 0
        set(gcf,'Position',[50, screensize(4)-Sim.figuresize(2)-150, Sim.figuresize(1), Sim.figuresize(2)]);
    else
    end

%% Sets the rotation axis based on the input

angle = Sim.angle/360*2*pi;
theta=angle/(Sim.nframes);

if (Sim.direction == 'x');
    rotmat = [1 0 0;0 cos(-theta) -sin(-theta);0 sin(-theta) cos(-theta)];
elseif (Sim.direction == 'y');
    rotmat = [cos(-theta) 0 sin(-theta);0 1 0;-sin(-theta) 0 cos(-theta)];
elseif (Sim.direction == '-x');
    rotmat = [1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
elseif (Sim.direction == '-y');
    rotmat = [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
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

time = 0:Sim.angle/Sim.nframes:Sim.angle;
xsum = zeros(1,numel(time));
ysum = zeros(1,numel(time));
zsum = zeros(1,numel(time));

if strcmp(Sim.video,'n');
%% Performs the incremental rotation and plots each frame with no animation output

for j=1:Sim.nframes+1
    if j ~= 1
        coordinates = coordinates*rotmat;
    end
    ends=coordinates;
	xsum(1,j)=round(sum(ends(:,1))/mdegrees(1)*100)/100;
	ysum(1,j)=round(sum(ends(:,2))/mdegrees(1)*100)/100;
    zsum(1,j)=round(sum(ends(:,3))/mdegrees(1)*100)/100;
    subplot('Position',[0.03 0.03 0.4 0.9]);
    quiver3(starts(:,1),starts(:,2),starts(:,3),coordinates(:,1),coordinates(:,2),coordinates(:,3));
    hold on
    plot3([1 0],[0 0],[0 0],'k')
    plot3([0 0],[1 0],[0 0],'k')
    plot3([0 0],[0 0],[1 0],'k')
    text(1.1,0,0,'X_R')
    text(0,1.1,0,'Y_R')
    text(0,0,1.1,'Z_R')
    hold off
    set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
    title('Vector Diagram Simulation')
    axis ([-1 1 -1 1 -1 1]);
    subplot('Position',[0.5 0.76 0.47 0.2]);
    plot(time,xsum)
    axis([min(time) max(time) -1 1])
    ylabel('M_X');
    xlabel('degree rotation');
    subplot('Position',[0.5 0.43 0.47 0.2]);
    plot(time,ysum)
    axis([min(time) max(time) -1 1])
    ylabel('M_Y');
    xlabel('degree rotation');
    subplot('Position',[0.5 0.1 0.47 0.2]);
    plot(time,zsum)
    axis([min(time) max(time) -1 1])
    ylabel('M_Z');
    xlabel('degree rotation');
    pause(1/Sim.framerate);
end

% Converts the rectangular coordinates back to spherical coordinates

r = zeros(1,mdegrees(1));

for k=1:mdegrees(1)
    r(k)=((ends(k,1)-0)^2+(ends(k,2)-0)^2+(ends(k,3)-0)^2)^0.5;
end

angles = zeros(mdegrees(1),2);

for l=1:mdegrees(1)
    angles(l,:)= [atan2(ends(l,2),ends(l,1))*180/pi acos(ends(l,3)/r(l))*180/pi];
end

Spinsout = [angles r'];

elseif strcmp(Sim.video,'gif');
%% Performs the incremental rotation and plots each frame with animated GIF output

for j=1:Sim.nframes+1
    if j ~= 1
        coordinates = coordinates*rotmat;
    end
    ends=coordinates;
	xsum(1,j)=round(sum(ends(:,1))/mdegrees(1)*100)/100;
	ysum(1,j)=round(sum(ends(:,2))/mdegrees(1)*100)/100;
    zsum(1,j)=round(sum(ends(:,3))/mdegrees(1)*100)/100;
    subplot('Position',[0.03 0.03 0.4 0.9]);
    quiver3(starts(:,1),starts(:,2),starts(:,3),coordinates(:,1),coordinates(:,2),coordinates(:,3));
    hold on
    plot3([1 0],[0 0],[0 0],'k')
    plot3([0 0],[1 0],[0 0],'k')
    plot3([0 0],[0 0],[1 0],'k')
    text(1.1,0,0,'X_R')
    text(0,1.1,0,'Y_R')
    text(0,0,1.1,'Z_R')
    hold off
    set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
    title('Vector Diagram Simulation')
    axis ([-1 1 -1 1 -1 1]);
    subplot('Position',[0.5 0.76 0.47 0.2]);
    plot(time,xsum)
    axis([min(time) max(time) -1 1])
    ylabel('M_X');
    xlabel('degree rotation');
    subplot('Position',[0.5 0.43 0.47 0.2]);
    plot(time,ysum)
    axis([min(time) max(time) -1 1])
    ylabel('M_Y');
    xlabel('degree rotation');
    subplot('Position',[0.5 0.1 0.47 0.2]);
    plot(time,zsum)
    axis([min(time) max(time) -1 1])
    ylabel('M_Z');
    xlabel('degree rotation');
    frame = getframe(gcf);
    image = frame2im(frame);
    [imind,cm] = rgb2ind(image,256);
    if exist(Sim.filename,'file') == 0;
        imwrite(imind,cm,Sim.filename,'gif','Loopcount',Sim.numloops,'DelayTime',1/Sim.framerate);
        Sim.gifout = 1;
    else
        if (isfield(Sim,'gifout') == 0);
            error('File already exists!');
        else
            imwrite(imind,cm,Sim.filename,'gif','WriteMode','append','DelayTime',1/Sim.framerate);
        end
    end
end

% Converts the rectangular coordinates back to spherical coordinates

r = zeros(1,mdegrees(1));

for k=1:mdegrees(1)
    r(k)=((ends(k,1)-0)^2+(ends(k,2)-0)^2+(ends(k,3)-0)^2)^0.5;
end

angles = zeros(mdegrees(1),2);

for l=1:mdegrees(1)
    angles(l,:)= [atan2(ends(l,2),ends(l,1))*180/pi acos(ends(l,3)/r(l))*180/pi];
end

Spinsout = [angles r'];

elseif strcmp(Sim.video,'avi');
%% Performs the incremental rotation and plots each frame with AVI output

    if exist(Sim.filename,'file') == 0;

        aviobj = VideoWriter(Sim.filename);
        open(aviobj)

        for j=1:Sim.nframes+1
            if j ~= 1
                coordinates = coordinates*rotmat;
            end
            ends=coordinates;
        	xsum(1,j)=round(sum(ends(:,1))/mdegrees(1)*100)/100;
        	ysum(1,j)=round(sum(ends(:,2))/mdegrees(1)*100)/100;
            zsum(1,j)=round(sum(ends(:,3))/mdegrees(1)*100)/100;
            subplot('Position',[0.03 0.03 0.4 0.9]);
            quiver3(starts(:,1),starts(:,2),starts(:,3),coordinates(:,1),coordinates(:,2),coordinates(:,3));
            hold on
            plot3([1 0],[0 0],[0 0],'k')
            plot3([0 0],[1 0],[0 0],'k')
            plot3([0 0],[0 0],[1 0],'k')
            text(1.1,0,0,'X_R')
            text(0,1.1,0,'Y_R')
            text(0,0,1.1,'Z_R')
            hold off
            set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
            title('Vector Diagram Simulation')
            axis ([-1 1 -1 1 -1 1]);
            subplot('Position',[0.5 0.76 0.47 0.2]);
            plot(time,xsum)
            axis([min(time) max(time) -1 1])
            ylabel('M_X');
            xlabel('degree rotation');
            subplot('Position',[0.5 0.43 0.47 0.2]);
            plot(time,ysum)
            axis([min(time) max(time) -1 1])
            ylabel('M_Y');
            xlabel('degree rotation');
            subplot('Position',[0.5 0.1 0.47 0.2]);
            plot(time,zsum)
            axis([min(time) max(time) -1 1])
            ylabel('M_Z');
            xlabel('degree rotation');
            aviframes(j) = getframe(gcf);
        end

        writeVideo(aviobj,aviframes)
        close(aviobj)
        close(gcf)
 
        % Converts the rectangular coordinates back to spherical coordinates
        
        r = zeros(1,mdegrees(1));
        
        for k=1:mdegrees(1)
            r(k)=((ends(k,1)-0)^2+(ends(k,2)-0)^2+(ends(k,3)-0)^2)^0.5;
        end

        angles = zeros(mdegrees(1),2);
        
        for l=1:mdegrees(1)
            angles(l,:)= [atan2(ends(l,2),ends(l,1))*180/pi acos(ends(l,3)/r(l))*180/pi];
        end

        Sim.aviout = aviframes;
        Spinsout = [angles r'];

    else

        
        if (isfield(Sim,'aviout') == 0);
            error('File already exists!');
        end
        
        aviobj = VideoWriter(Sim.filename);
        open(aviobj)
        aviframes = Sim.aviout;
        framenumber = size(Sim.aviout);

        for j=1:Sim.nframes+1
            if j ~= 1
                coordinates = coordinates*rotmat;
            end
            ends=coordinates;
            xsum(1,j)=round(sum(ends(:,1))/mdegrees(1)*100)/100;
            ysum(1,j)=round(sum(ends(:,2))/mdegrees(1)*100)/100;
            zsum(1,j)=round(sum(ends(:,3))/mdegrees(1)*100)/100;
            subplot('Position',[0.03 0.03 0.4 0.9]);
            quiver3(starts(:,1),starts(:,2),starts(:,3),coordinates(:,1),coordinates(:,2),coordinates(:,3));
            hold on
            plot3([1 0],[0 0],[0 0],'k')
            plot3([0 0],[1 0],[0 0],'k')
            plot3([0 0],[0 0],[1 0],'k')
            text(1.1,0,0,'X_R')
            text(0,1.1,0,'Y_R')
            text(0,0,1.1,'Z_R')
            hold off
            set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
            title('Vector Diagram Simulation')
            axis ([-1 1 -1 1 -1 1]);
            subplot('Position',[0.5 0.76 0.47 0.2]);
            plot(time,xsum)
            axis([min(time) max(time) -1 1])
            ylabel('M_X');
            xlabel('degree rotation');
            subplot('Position',[0.5 0.43 0.47 0.2]);
            plot(time,ysum)
            axis([min(time) max(time) -1 1])
            ylabel('M_Y');
            xlabel('degree rotation');
            subplot('Position',[0.5 0.1 0.47 0.2]);
            plot(time,zsum)
            axis([min(time) max(time) -1 1])
            ylabel('M_Z');
            xlabel('degree rotation');
            aviframes(j+framenumber(2)) = getframe(gcf);
        end
    
        writeVideo(aviobj,aviframes)
        close(aviobj)
        close(gcf)
        
        % Converts the rectangular coordinates back to spherical coordinates
        
        r = zeros(1,mdegrees(1));
        
        for k=1:mdegrees(1)
            r(k)=((ends(k,1)-0)^2+(ends(k,2)-0)^2+(ends(k,3)-0)^2)^0.5;
        end

        angles = zeros(mdegrees(1),2);
        
        for l=1:mdegrees(1)
            angles(l,:)= [atan2(ends(l,2),ends(l,1))*180/pi acos(ends(l,3)/r(l))*180/pi];
        end

        Sim.aviout = aviframes;
        Spinsout = [angles r'];

    end
end
