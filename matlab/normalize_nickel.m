%% normalize nickel data to nitrogen


avg_1st_nitro_pt = 0;
avg_1st_nickel_pt = 0;

if (numel(condition(3)) ~=0 & numel(condition(1)) ~=0)
    
%normalization will only execute if there is at least 1 nickel replicate
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
