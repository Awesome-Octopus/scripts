%  euler2quat  Convert Euler angles to unit quaternions.
%
%  q = euler2quat(Angles);
%  q = euler2quat(alpha, beta, gamma);
%
%  Input:
%     Angles         numeric array, size = (3,...)
%                    the rows correspond to first, second and third Euler
%                    angle(s) (alpha, beta, and gamma, respectively), in radians
%     alpha          double or numeric array, size = (1,...)
%                    first Euler angle(s), in radians
%     beta           double or numeric array, size = (1,...)
%                    second Euler angle(s), in radians
%     gamma          double or numeric array, size = (1,...)
%                    third Euler angle(s), in radians
%
%  Output:
%     q              numeric, size = (4,...)
%                    normalized quaternion(s)
