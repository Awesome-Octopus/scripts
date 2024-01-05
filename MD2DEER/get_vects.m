function position_vectors = get_vects (index, coordinate)
% returns an n X 3 matrix where n is the number of frames
% containing xyz for each frame of coordinate(index)
        position_vectors = zeros(numel(coordinate(index).xyz), 3)
        for i=1:numel(position_vectors)
            position_vectors(i,:) = {coordinate(index).xyz(i)};
        end
end