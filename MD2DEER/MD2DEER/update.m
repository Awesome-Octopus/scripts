function update (fig, coordinate, numframes)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    W A R N I N G
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if more ui elements are added to the program in the future, then the
    % list of GUI elements in fig.Children also changes, and the way you 
    % reference the selected values from that element will have to change
    % too. eg: if a toggle switch is added it will be accessed through
    % fig.Children (end - 6) and so the handle to the plotting axis will
    % become fig.Children(end-7). You will need to manually update the 
    % indeces of all instances wherein you check a gui element 
    % for a value or input
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % each time we added a site, 5 graphics objects were created and added 
    % to fig.children, in addition to the button ui objects, ui labels, and
    % the axes we made. Therefore there will be (n*5 + 7) elements in the
    % cell array fig.Children and the indices are as follows:
    % 
    % --- slider for last site pair added - 1
    % --- chain select dd for site 2 of last site pair added - 2
    % --- chain select dd for site 1 of last site pair added - 3
    % --- seq select dd for site 2 of last site pair added - 4
    % --- seq select dd for site 1 of last site pair added - 5
    % ...
    % --- slider for first site pair added - (end-11)
    % --- chain select dd for site 2 of first site pair added - (end-9)
    % --- chain select dd for site 1 of first site pair added - (end-10)
    % --- seq select dd for site 2 of first site pair added - (end-7)
    % --- seq select dd for site 1 of first site pair added - (end-8)
    % n is the number of site pairs selected as judged by the number of
    % entries in fig. children based on assignment pattern above
    n = (numel(fig.Children) - 7)/5;
    
    % pre-allocate an array for the positions
    % dimensions are (frame, (x y or z), label 1 or 2, label pair number)
    label_positions = zeros (numframes, 3, 2, n);
   
    % for each site pair
    for i = 1:n
        
            % for each member of the pair
            for p = 1:2
                
        
                % have this function choose a method to calculate the 
                % distances, then hand off to another funtion which does the
                % calculating, and finally return the positions.
    
    
                label_choice = fig.Children(end-3).Value;
                seq = fig.Children(end-7-5*(i-1)-(p-1)).Value;
                chain = fig.Children(end-9-5*(i-1)-(p-1)).Value;
                

              
                if strcmpi(label_choice, 'MTSL') == true

                    label_positions(:,:,p,i) = get_MTSL_position(chain, seq,...
                        coordinate, numframes);

                    % if there is no sidechain at this residue the above
                    % function will return NaN, in which case treat it as
                    % TOAC and notify the user.
                    if isnan(label_positions(:,:,p,i)) == true
                        uialert(fig,"One or more selected residues does not have a side chain. Estimate of label at this position will be based on backbone",...
                        'Warning', 'Icon', 'warning');

                        label_positions(:,:,p,i) = get_TOAC_position(chain, seq, ...
                        coordinate, numframes);

                    end    

                elseif strcmpi(label_choice, 'TOAC') == true

                    label_positions(:,:,p,i) = get_TOAC_position(chain, seq, ...
                        coordinate, numframes);

                    % if the selected residue is the first or last one on 
                    % the chain, tell the user you cant TOAC distances.
                    if sum(isnan(label_positions)) ~= 0
                        uialert(fig,"One or more selected residues does not have a side chain. Estimate of label at this position will be based on backbone",...
                        'Warning', 'Icon', 'warning');

                        label_positions(:,:,p,i) = get_TOAC_position(chain, seq, ...
                        coordinate, numframes);

                    end    


                elseif strcmpi(label_choice, 'Backbone') == true

                    label_positions(:,:,p,i) = get_xyz(coordinate, chain, seq,...
                        'BB', numframes);
               
                else
                    uialert(fig,"Oops! This shouldn't ever happen! (function: update)",...
                        'Error', 'Icon', 'error');
                end
            end
    end

    
    distances = zeros (numframes, n); 

    for i = 1:n
        for p = 1:numframes
            dv = label_positions(p, :, 2, i) - label_positions(p, :, 1, i);
            
            distances(p,i) = norm(dv); 
        end    
    end
    
    fig.UserData = distances;


    %% change these to affect the way histograms are displayed
    
    % the longest distance among all distances in all sites, rounded up
    max_dist = max(ceil(max(distances)));
    
    % arbitrary
    nbins = 40;
    binSz = max_dist/nbins;
    edges = 0:binSz:max_dist;
   
    %%
    %take the distances and bin them for a histogram, for each site, counts
    %and bins for each label site are stored in a cell array where they can
    %be weighted seperately from input from the sliders
    counts = cell(n, 1);
    

    display_component_dist = fig.Children(end-5).Value;
    %if the user wants to see the underlying distance components from spin
    %label pairs overlayed without summing
    if display_component_dist == 1
        
        for i=1:n
             counts{i} = histcounts(distances(:,i), edges, ...
                "Normalization",'pdf');
             counts{i} = counts{i}*fig.Children(1+((i-1)*5)).Value;
             bar(fig.Children(end), edges(1:end-1), counts{i}, ...
                 'EdgeAlpha', 0.7, 'FaceAlpha', 0.7, 'BarWidth', 1);
             hold(fig.CurrentAxes, "on");
        end
    
    % otherwise we want to see the combined probability density which is
    % what the resulting distance distribution will actually look like,
    % this should be the default behaviour
    else
        combined_cnts = zeros(1, nbins); % <--there should be a better way 


        for i=1:n
            counts{i} = histcounts(distances(:,i), edges, ...
                "Normalization",'pdf');
            counts{i} = counts{i}*fig.Children(1+((i-1)*5)).Value;
            for s=1:numel(combined_cnts)
                combined_cnts(s) = combined_cnts(s) + counts{i}(s);
            end
        end
            
        bar(fig.Children(end), edges(1:end-1), combined_cnts, 'hist');
            hold(fig.CurrentAxes, "on");
            
    end

    hold(fig.CurrentAxes, "off");
    fig.CurrentAxes.XTickMode = "auto";
end