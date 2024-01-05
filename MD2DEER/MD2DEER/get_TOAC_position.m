function  toac_pos = get_TOAC_position (chain, seq, coordinate, numframes)
    
    s = str2double(seq);
    
    %% these parameters will need to be adjusted emperically

    theta = pi*3/5;
    scalar = 1.5;


    %% check and ensure that s is not 1 or the C-term residue 
    
     % implement this later
    

    % we first need to get the vectors to the backbone of the bead
    % preceding it and succeeding it. Makes a numframes x 3 x 3 matrix
    prev = s - 1;
    prev = num2str(prev);
    next = s + 1;
    next = num2str(next);
   
    bb_pos = zeros(numframes, 3, 3);
    bb_pos(:,:,1) = get_xyz(coordinate, chain, prev, 'BB', numframes);
    bb_pos(:,:,2) = get_xyz(coordinate, chain, seq, 'BB', numframes);
    bb_pos(:,:,3) = get_xyz(coordinate, chain, next, 'BB', numframes);
    

    toac_pos = zeros(numframes, 3);
    %now for each set of three vectors per frame, apply a rotation matrix
    for i = 1:numframes
        
        toac_pos(i,:) = rotate_in_plane(bb_pos(i,:,1), bb_pos(i,:,2), bb_pos(i,:,3), theta, scalar);
    end
end 