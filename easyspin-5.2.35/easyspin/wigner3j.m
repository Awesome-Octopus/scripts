% wigner3j   Wigner 3-j symbol 
%
%   v = wigner3j(j1,j2,j3,m1,m2,m3)
%   v = wigner3j(jm1,jm2,jm3)
%   v = wigner3j(jjj,mmm)
%   v = wigner3j(jjjmmm)
%
%   Computes the value of the Wigner 3-j symbol
%
%      / j1  j2  j3 \
%      |            |
%      \ m1  m2  m3 /
%
%   Definitions for alternative input forms
%   a)  jm1 = [j1 m2], jm2 = [j2 m2], jm3 = [j3 m3]
%   b)  jjj = [j1 j2 j3], mmm = [m1 m2 m3]
%   c)  jjjmmm = [j1 j2 j3 m1 m2 m3]
