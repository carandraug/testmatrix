function A = hanowa(n, d)
%HANOWA  A matrix whose eigenvalues lie on a vertical line in the complex plane.
%        HANOWA(N, d) is the N-by-N block 2x2 matrix (thus N = 2M must be even)
%                      [d*EYE(M)   -DIAG(1:M)
%                       DIAG(1:M)   d*EYE(M)]
%        It has complex eigenvalues lambda(k) = d +/- k*i  (1 <= k <= M).
%        Parameter d defaults to -1.

%        Reference:
%        E. Hairer, S.P. Norsett and G. Wanner, Solving Ordinary
%        Differential Equations I: Nonstiff Problems, Springer-Verlag,
%        Berlin, 1987. (pp. 86-87)

if nargin == 1, d = -1; end

m = n/2;
if round(m) ~= m
   error('N must be even.')
end

A = [ d*eye(m) -diag(1:m)
      diag(1:m)  d*eye(m)];
