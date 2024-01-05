function addsite (~, fig, seqlist, chainlist)
    % this will add to the GUI 4 dropdown menus which allow
    % for the selection of chain and residue at 2 sites, along with a
    % slider that can scale down the histogram counts to be used to change
    % weighting of contribution of spin label pairs to the modelled
    % distance distribution.

    %keep track of how many times this function has been called
    persistent Npair;
    if isempty(Npair)
        Npair = 1;
    else    
        Npair = Npair + 1;
    end
    
    uidropdown (fig,"Items",seqlist,"ItemsData",seqlist,...
        'Position', [150 (480 - Npair*80) 60 22]);
    uidropdown (fig,"Items",seqlist,"ItemsData",seqlist,...
        'Position', [375 (480 - Npair*80) 60 22]);

    uidropdown (fig,"Items",chainlist,"ItemsData",chainlist,...
        'Position', [75 (480 - Npair*80) 60 22]);
    uidropdown (fig,"Items",chainlist,"ItemsData",chainlist,...
        'Position', [300 (480 - Npair*80) 60 22]);

    uislider(fig,"Limits",[0 1],"Position",[90, (470-Npair*80), 325 3],...
        'Value', 1, 'Tooltip', 'Adjust the weighting of the distance distribution relative to other label pairs');


end