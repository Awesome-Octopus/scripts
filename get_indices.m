function indices = get_indices (fig, coordinate, n)

% an index vector will be made containing the index of the coordinate element
% that matches the given tags for selecting site 1 and then site 2 

    indices = zeros(n, 4);
    % indices is a matrix with each row being a label pair
    % [site1bb site2bb site1sc (optional) site2sc (optional)]

   

    % for each site pair
    for p=1:n

        c1 = fig.Children(end-6-(p-1)*5).Value;
        c2 = fig.Children(end-7-(p-1)*5).Value;
        s1 = fig.Children(end-4-(p-1)*5).Value;
        s2 = fig.Children(end-5-(p-1)*5).Value;
        % check the tags on each coordinate

        for i=1:numel(coordinate)

            if strcmpi({coordinate(i).chain}, c1) == true;

                if strcmpi({coordinate(i).seq}, s1) == true;

                    if strcmpi({coordinate(i).type}, 'BB') == true;
                        indices(p,1) = i;

                    else strcmpi({coordinate(i).type}, 'SC1') == true;
                        indices(p,3) = i;
                    end
                end
            end
            if strcmpi({coordinate(i).chain}, c2) == true;

                if strcmpi({coordinate(i).seq}, s2) == true;
                    if strcmpi({coordinate(i).type}, 'BB') == true;
                        indices(p,2) = i;
                    else strcmpi({coordinate(i).type}, 'SC1') == true;
                        indices(p,4) = i;
                    end
                end
            end
        end
    end
end
