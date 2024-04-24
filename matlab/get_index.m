function index = get_index (seq, chain, coordinate, type)
% return the index for the coordinate 
% seq is the selected value from a dropdown component handle
% chain is the selected value from a dropdown component handle
% type is either 'BB' or 'SC1' (could be anything though)
    
    for i = 1:numel(coordinate)
        if strcmpi({coordinate(i).seq}, seq) == true
            if strcmpi({coordinate(i).chain}, chain) == true
                if strcmpi({coordinate(i).type}, type) == true
                        index = i;
                else
                    uialert(uifigure(),"Oops, this isn't supposed to happen coordinate(i).type does not match (function get_index)", 'Error','CloseFcn', 'Icon', 'error'); 
                end
            end
        end
    end    
end