function IND = sub_get_ARCgrid;
% Subset topo/grid for
% specified region
% from ARCc grid
fprintf('Subsetting region from ARCc ...\n');
% Indices for Arctic & N.Atl subregions in ARCc grid:
IND.i1=213;
IND.i2=1580;
IND.j1=150;
IND.j2=1980;

return