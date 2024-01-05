%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A highly modified version of Rob Mckarrick's powersat.m script     %
% that allows the air, nickel, and nitrogen conditions to be plotted %
% together. You must also have powersat_multiload.m and              %
% powersat_multifit.m in your matlab folder                          %                                 %
%                                                                    %
% By Andrew Morris, Lorigan Lab                                      %
% akmorris1@gmail.com or morri361@miamioh.edu                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DPPH_constant
DPPH_constant = 14.0437;
% This value is the half saturation value of crystalline DPPH powder in a
% TPX tube equilibrated with nitrogen divided by the peak to peak central
% linewidth of those spectra at pre-saturating power and should be a
% universal constant so long as the same instrument an resonator is used.
% it can be updated if this value needs to be measured again later.

%% powersat multiplot

condition = powersat_multifit;
% takes the average of fitted p1/2 from each rep in each condition 
% second element in the x variable for each rep in each condition is 
% the fitted p1/2 value

condition_label = {'Nitrogen', 'Air', 'Nickel'};

%% Normalize nickel
avg_1st_nickel_pt = 0;
avg_1st_nitro_pt = 0;

if (numel(condition(3)) ~= 0 && numel(condition(1)) ~=0)
    disp (numel(condition(3)));
% normalization will only execute if there is at least 1 nickel replicate
% and at least 1 nitro replicate

% unfortunately this has to be done this way because matlab does not let
% you index through structures like condition(1).rep.height(1)

    for p = 1:numel(condition(1).rep) % average the first height point in the nitro data reps
   
        avg_1st_nitro_pt = avg_1st_nitro_pt + condition(1).rep(p).height(1); 
    
    end
 
    avg_1st_nitro_pt = avg_1st_nitro_pt/numel(condition(1).rep);

    for p = 1:numel(condition(3).rep) % average the first height point in the nickel data reps
   
        avg_1st_nickel_pt = avg_1st_nickel_pt + condition(3).rep(p).height(1); 
    
    end
 
    avg_1st_nickel_pt = avg_1st_nickel_pt/numel(condition(3).rep);


    for p = 1:numel(condition(3).rep) % scale each nickel replicate based on the 
    % ratio of the average of the highest attenuation point for nitro
    % vs nickel
   
        condition(3).rep(p).height = condition(3).rep(p).height*...
        avg_1st_nitro_pt/avg_1st_nickel_pt; 

        condition(3).rep(p).yfitted = condition(3).rep(p).yfitted*...
        avg_1st_nitro_pt/avg_1st_nickel_pt; 

    end
end


for i = 1:numel(condition) % for each condition
    if isempty(condition(i).rep) == 0 
    
        % make an array containing the fitted p12s for each replicate
        for p = 1:numel(condition(i).rep) 
            condition(i).p12s(p) = condition(i).rep(p).x(2);
        end
    
        condition(i).mean_p12 = mean(condition(i).p12s); %find the average p1/2
        condition(i).p12_stdev = std(condition(i).p12s); % and the std dev
    
    
        for p = 1:numel(condition(i).rep) % plots each replicate
        
            plot(condition(i).rep(p).sqrtpower, condition(i).rep(p).height,...
            'o', condition(i).rep(p).xdata, condition(i).rep(p).yfitted,...
            condition(i).color, 'MarkerEdgeColor',condition(i).color',...
            'MarkerFaceColor', 'none');
            % plots the data and the fit
            
            axis tight; % makes the axis fit the data tightly
        
            xlabel('square root power (mw)','FontSize',14),...
            ylabel('peak to trough height','FontSize',14); 
            % Constructs the axes labels
        
            title('Power Saturation Curve','FontSize',16);
            % Constructs the title
        
            hold on;
        
        end
    
        powerlabel = [sprintf('%s ', string(condition_label(i))),...
            sprintf('P_1_/_2 = %.1f', condition(i).mean_p12)...
            ,' \pm ', sprintf('%.1f', condition(i).p12_stdev), ' mw'];
        
        
        
        
        text(0.95,0.2-0.05*i, powerlabel...
            ,'FontSize', 8, 'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'bottom', 'Units', 'normalized',...
            'color', condition(i).color); 
            % Prints the T1 string on the plot

        % NOTE: the code is slightly different from the original in that 
        % it draws the line to the mean p1/2 for that condition rather than
        % to the p1/2 for that replicate    
        
        line([condition(i).mean_p12^0.5;condition(i).mean_p12^0.5],...
            [0;max(condition(i).rep(p).height)],'Color',condition(i).color,...
            'LineWidth',2);
        % Draws a line at the T1 value 
         
        hold on;
        
    end     
end 


for i=1:numel(condition)
    for p=1:numel(condition(i).rep)
        condition(i).subsat_linewidths(p,1:6) = condition(i).rep(p).delH(1:6);
        % semilogx(condition(i).rep(p).power, condition(i).rep(p).delH);
        hold on;
    end
        condition(i).meandelH = mean(mean(condition(i).subsat_linewidths));
    if i > 1
        condition(i).Pi = (condition(i).mean_p12 - condition(1).mean_p12)...
            /condition(i).meandelH/DPPH_constant;
    end 
end
  
%% find the depth parameter
if numel(condition) ~=0
    try     
        phi = reallog(condition(2).Pi/condition(3).Pi);

        phi_label = [sprintf('\\phi = %.1f', phi)];
        
        % gets phi and prepares a label
    catch
        phi_label = [sprintf('\\phi not calculable')];
        % if air or nickel is below nitrogen, it attempts to take log of a
        % negative number and returns an error, this catches the error and
        % prints that phi cannot be calculated

    end



    text(0.95, 0.4, phi_label...
        ,'FontSize', 20, 'color', 'k', 'HorizontalAlignment', 'right',...
        'VerticalAlignment', 'bottom', 'Units', 'normalized');       
    
        %print the calculated depth parameter
end

%% print Pi

if numel(condition) ~=0
    
    if numel(condition) > 1
        
        for i = 2:numel(condition)  
        
               
            pilabel = sprintf('\\Pi_{%s} = %.3f', condition_label{i}, condition(i).Pi);
            text(0.95, 0.4-0.05*i, pilabel ,'FontSize', 10, 'HorizontalAlignment', 'right',...
        'VerticalAlignment', 'bottom', 'Units', 'normalized', 'color', condition(i).color);  
        
              % labels graph with the calculated accessibility parameter
        end    
    end
end

