function mtsl_pos = get_MTSL_position (chain, seq, coordinate, numframes)

    % this follows a very simple aproximation of MTSL as a lengthening of
    % the vector from the backbone to the sidechain by a scalar. The best
    % scalar value to use needs to eb determined emperically
    scalar = 2.5;


    sc_pos = get_xyz(coordinate, chain, seq, 'SC1', numframes);
    

    % array for the x y z of the label at each frame 
    mtsl_pos = zeros(numframes, 3);

    % if the user selected a glycine or otherwise the selected residue
    % doesnt have an SC1 bead, send an error then extrapolate from sidechain,
    % treating it as though it were toac
    if sum(isnan(sc_pos)) ~= 0
                
        mtsl_pos = NaN;
    else
        
        bb_pos = get_xyz(coordinate, chain, seq, 'BB', numframes);
        mtsl_pos = (sc_pos - bb_pos)*scalar;
        mtsl_pos = bb_pos + mtsl_pos;
    end
end