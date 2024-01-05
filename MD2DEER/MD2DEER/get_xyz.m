function xyz = get_xyz (coordinate, chain, seq, type, numframes)

    % given a string seq, and a string chain, and a string type ('BB' or
    % 'SC1'
    % return the all xyz coodinates that match those fields in coordinate
    % across all frames in a N x 3 matrix
    xyz = zeros(numframes, 3);

    % if an assignment is made this is changed to true
    matched_position = false;

    for i=1:numel(coordinate)
        if strcmpi(coordinate(i).seq, seq) == true
            if strcmpi(coordinate(i).chain, chain) == true
                if strcmpi(coordinate(i).type, type) == true
                    for p=1:numframes
                        
                        xyz(p,:) = coordinate(i).xyz{p};
                    end
                    matched_position = true;
                    return;                  
                end    
            end
        end
    end
    if matched_position == false
        xyz(:,:) = NaN;
    end    
end