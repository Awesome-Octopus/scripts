fig = uifigure('Position',[500 500 300 200]);
dd = uidropdown(fig, ...
    'Editable','on', ...
    'Items',{'Red','Green','Blue'}, '', @(src, event)disp("we got there boys"));
