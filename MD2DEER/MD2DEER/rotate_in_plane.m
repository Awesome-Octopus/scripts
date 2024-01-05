function result_vect = rotate_in_plane( v1, v2, v3, theta, scalar)

    % rotates the 3d point v3 arround v2 in the plane formed by v1,v2,v3
    % by the angle theta (in radians) according to the right hand rule
    
    A = v3 - v2;
    B = v1 - v2;
    len_A = sqrt(dot(A, A));
    len_B = sqrt(dot(B, B));
    norm = cross((A/len_A), (B/len_B));
    
    % this will shift our coordinate system so that norm is on the 
    % xz plane by rotating arround the z axis
    
    % the length of the projection of the unit normal vector onto the xy
    % plane
    proj_len = sqrt(norm(1)^2 + norm(2)^2);
    
    xy_rot_m = [ norm(1)/proj_len -1*norm(2)/proj_len 0; ...
        norm(2)/proj_len norm(1)/proj_len 0; 0 0 1];
    
    norm = norm*xy_rot_m;
    
    % apply the same matrix to your vector A
    A = A*xy_rot_m;
    
    % the normal vector is now in the xz plane. now rotate arround the y
    % axis so that norm = [0 0 1]
    xz_rot_m = [ norm(3) 0 norm(1); 0 1 0; -1*norm(1) 0 norm(3) ];
  
    
    A = A*xz_rot_m;

    % A is now in the xy plane. we can now rotate by theta about the normal
    % axis (colinear with z)

    theta_rot_m = [cos(theta) -1*sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    
    % scalar being the factor we want to extend our rotated vector by
    A = A*theta_rot_m*scalar;
    % now we just invert the initial two rotation matrices we used.
    
    A = A*inv(xz_rot_m);
    result_vect = A*inv(xy_rot_m);
    
    %add back the original placement of the backbone position
    result_vect = result_vect + v2;

end