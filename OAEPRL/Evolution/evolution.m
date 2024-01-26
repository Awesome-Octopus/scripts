%% Spin Evolution Function
% Ohio Advanced EPR Laboratory
% Robert McCarrick
% Vector Simulation Package, Version 2.01
% Last Modified 10/31/2017

function [Spinsout Sim]=evolution(Spins,Sim)


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

if strcmp(Sim.video,'gif') | strcmp(Sim.video,'avi');
    if (isfield(Sim,'filename') == 0);
        error('Please specify filename for animation.');
    end
end

if (isfield(Sim,'nframes') == 0);
    Sim.nframes=100;
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


if strcmp(Sim.video,'n');    
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
    subplot('Position',[0.03 0.03 0.4 0.9]);
    quiver3(starts(:,1),starts(:,2),starts(:,3),ends(:,1),ends(:,2),ends(:,3));
    hold on
    plot3([1 0],[0 0],[0 0],'k')
    plot3([0 0],[1 0],[0 0],'k')
    plot3([0 0],[0 0],[1 0],'k')
    text(1.1,0,0,'X_R')
    text(0,1.1,0,'Y_R')
    text(0,0,1.1,'Z_R')
    hold off
    set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
    title('Vector Diagram Simulation');
    axis ([-1 1 -1 1 -1 1]);
    subplot('Position',[0.5 0.76 0.47 0.2]);
    plot(time,xsum);
    ylabel('M_X');
    xlabel('time (s)');
    axis([min(time) max(time) -1 1]);
    subplot('Position',[0.5 0.43 0.47 0.2]);
    plot(time,ysum);
    ylabel('M_Y');
    xlabel('time (s)');
    axis([min(time) max(time) -1 1]);
    subplot('Position',[0.5 0.1 0.47 0.2]);
    plot(time,zsum);
    ylabel('M_Z');
    xlabel('time (s)');
    axis([min(time) max(time) -1 1]);
    pause(1/Sim.framerate);
end

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

elseif strcmp(Sim.video,'gif');
%% Performs the incremental rotation and plots each frame with animated GIF output

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
    subplot('Position',[0.03 0.03 0.4 0.9]);
    quiver3(starts(:,1),starts(:,2),starts(:,3),ends(:,1),ends(:,2),ends(:,3));
    hold on
    plot3([1 0],[0 0],[0 0],'k')
    plot3([0 0],[1 0],[0 0],'k')
    plot3([0 0],[0 0],[1 0],'k')
    text(1.1,0,0,'X_R')
    text(0,1.1,0,'Y_R')
    text(0,0,1.1,'Z_R')
    hold off
    set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
    title('Vector Diagram Simulation');
    axis ([-1 1 -1 1 -1 1]);
    subplot('Position',[0.5 0.76 0.47 0.2]);
    plot(time,xsum);
    ylabel('M_X');
    xlabel('time (s)');
    axis([min(time) max(time) -1 1]);
    subplot('Position',[0.5 0.43 0.47 0.2]);
    plot(time,ysum);
    ylabel('M_Y');
    xlabel('time (s)');
    axis([min(time) max(time) -1 1]);
    subplot('Position',[0.5 0.1 0.47 0.2]);
    plot(time,zsum);
    ylabel('M_Z');
    xlabel('time (s)');
    axis([min(time) max(time) -1 1]);
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
        relaxation = zeros(numel(mdegrees(1)),3);
        ends  = zeros(numel(mdegrees(1)),3);
        aviframes = zeros(1,numel(time));

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
            subplot('Position',[0.03 0.03 0.4 0.9]);
            quiver3(starts(:,1),starts(:,2),starts(:,3),ends(:,1),ends(:,2),ends(:,3));
            hold on
            plot3([1 0],[0 0],[0 0],'k')
            plot3([0 0],[1 0],[0 0],'k')
            plot3([0 0],[0 0],[1 0],'k')
            text(1.1,0,0,'X_R')
            text(0,1.1,0,'Y_R')
            text(0,0,1.1,'Z_R')
            hold off
            set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
            title('Vector Diagram Simulation');
            axis ([-1 1 -1 1 -1 1]);
            subplot('Position',[0.5 0.76 0.47 0.2]);
            plot(time,xsum);
            ylabel('M_X');
            xlabel('time (s)');
            axis([min(time) max(time) -1 1]);
            subplot('Position',[0.5 0.43 0.47 0.2]);
            plot(time,ysum);
            ylabel('M_Y');
            xlabel('time (s)');
            axis([min(time) max(time) -1 1]);
            subplot('Position',[0.5 0.1 0.47 0.2]);
            plot(time,zsum);
            ylabel('M_Z');
            xlabel('time (s)');
            axis([min(time) max(time) -1 1]);
            aviframes(i) = getframe(gcf);
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
            subplot('Position',[0.03 0.03 0.4 0.9]);
            quiver3(starts(:,1),starts(:,2),starts(:,3),ends(:,1),ends(:,2),ends(:,3));
            hold on
            plot3([1 0],[0 0],[0 0],'k')
            plot3([0 0],[1 0],[0 0],'k')
            plot3([0 0],[0 0],[1 0],'k')
            text(1.1,0,0,'X_R')
            text(0,1.1,0,'Y_R')
            text(0,0,1.1,'Z_R')
            hold off
            set(gca,'XDir','rev','YDir','rev','GridLineStyle','none','XTick',[],'YTick',[],'ZTick',[]);
            title('Vector Diagram Simulation');
            axis ([-1 1 -1 1 -1 1]);
            subplot('Position',[0.5 0.76 0.47 0.2]);
            plot(time,xsum);
            ylabel('M_X');
            xlabel('time (s)');
            axis([min(time) max(time) -1 1]);
            subplot('Position',[0.5 0.43 0.47 0.2]);
            plot(time,ysum);
            ylabel('M_Y');
            xlabel('time (s)');
            axis([min(time) max(time) -1 1]);
            subplot('Position',[0.5 0.1 0.47 0.2]);
            plot(time,zsum);
            ylabel('M_Z');
            xlabel('time (s)');
            axis([min(time) max(time) -1 1]);
            aviframes(i+framenumber(2)) = getframe(gcf);
        end        

        writeVideo(aviobj,aviframes)
        close(aviobj)
        close(gcf)

        % Converts the rectangular coordinates back to spherical
        % coordinates
        
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

