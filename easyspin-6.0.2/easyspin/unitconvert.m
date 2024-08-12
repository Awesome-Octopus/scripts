% unitconvert   Unit conversion tool
%
%   output = unitconvert(input, units)
%   output = unitconvert(input, units, g)
%
%   Converts the value in input, through the conversion specified in 
%   units returning the result in output.
%
%   units take the form of: 'unit_a->unit_b' where unit_a and unit_b are:
%   cm^-1, eV, K, mT or MHz
%
%   e.g. 'cm^-1->MHz' for conversion of wavenumbers to megahertz 
%
%                'cm^-1->eV'  'cm^-1->K'  'cm^-1->mT'  'cm^-1->MHz'       
%   'eV->cm^-1'               'eV->K'     'eV->mT'     'eV->MHz'
%   'K->cm^-1'   'K->eV'                  'K->mT'      'K->MHz'
%   'mT->cm^-1'  'mT->eV'     'mT->K'                  'mT->MHz'
%   'MHz->cm^-1' 'MHz->eV'    'MHz->K'    'MHz->mT'
%
%   When converting into or from magnetic field units, a g factor given 
%   as the third parameter is used. If it is not given, the g factor
%   of the free electron (gfree) is used.
%
%   input can be a vector of values. In this case, g
%   can be a scalar or a vector of the same size as input.
%
%   Example:
%     value_MHz = unitconvert(value_wn,'cm^-1->MHz')
%     value_mT = unitconvert(value_wn,'cm^-1->mT',2.005)
