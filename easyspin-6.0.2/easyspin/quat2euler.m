%  quat2euler  Convert a unit quaternion into a sequence of Euler angles
%              using the z-y-z convention.
%
%   [alpha, beta, gamma] = quat2euler(q);
%
%   Input:
%     q              numeric, size = (4,...)
%                    normalized quaternion
%
%   Output:
%    alpha          double or numeric, size = (1,...)
%                   first Euler angle, in radians
%
%    beta           double or numeric, size = (1,...)
%                   second Euler angle, in radians
%
%    gamma          double or numeric, size = (1,...)
%                   third Euler angle, in radians
