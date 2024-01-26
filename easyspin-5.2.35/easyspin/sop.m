% sop  Spin operator matrices
%
%   SpinOp = sop(SpinSystem, Comps)
%   [SpinOp1,SpinOp2,...] = sop(SpinSystem,Comps1,Comps2,...)
%   ... = sop(...,'sparse')
%
%   Spin operator matrix of the spin system
%   SpinSystem in the standard |mS,mI,...>
%   representation.
%
%   If more than one component is given, a matrix
%   is computed for each component.
%
%   Input:
%   - SpinSystem: vector of spin quantum numbers
%     or a spin system specification structure
%   - Comps: string specifying the operator, with several syntaxes
%      - specify component 'e','x','y','z','+','-' for each spin, indicating
%        E,Sx,Sy,Sz,S+,S-
%      - specify component and spin index, e.g. 'x2,z3'
%      - specify transition after component, e.g. 'x(1|3)' or 'x(1|3)2,z3' or
%        '+(1|2)1,e(3)2'
%
%   Output:
%   - SpinOp: operator matrix as requested
%
%   Examples:
%     SxIy = sop([1/2 1],'xy')    % returns SxIy for a S=1/2, I=1 system.
%
%     SeIp = sop([1/2 1/2],'e+')  % returns SeI+ for a S=I=1/2 system.
%
%     [Sx,Sy,Sz] = sop(1/2,'x','y','z')  % computes three matrices in one call.
%
%     Sxc = sop(5/2,'x(3|4)') % Sx on central transition -1/2<->+1/2
