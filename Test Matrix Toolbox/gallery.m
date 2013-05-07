## Copyright (C) 1989-1995 Nicholas .J. Higham
## Copyright (C) 2013 Carnë Draug
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {Function File} {} gallery (@var{name})
## @deftypefnx {Function File} {} gallery (@var{name}, @var{args})
## Create interesting matrices for testing.
##
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("cauchy", @var{x})
## @deftypefnx {Function File} {@var{c} =} gallery ("cauchy", @var{x}, @var{y})
## Create a Cauchy matrix.
##
## For vectors @var{x} and @var{y} of length N, creates an N-by-N
## matrix with @code{@var{c}(i,j) = 1/(@var{x}(i)+@var{y}(j))}.
## If unspecified, @var{y} defaults to @var{x}.  As a special case,
## if @var{x} is a scalar, it is interpreted as @code{1:@var{x}}.
##
## Explicit formulas are known for @code{det (@var{c})} (which is nonzero
## if @var{x} and @var{y} both have distinct elements) and the elements
## of @code{inv (@var{c})}.  @var{c} is totally positive if
## @code{0 < @var{x}(1) < @dots{} < @var{x}(N)} and
## @code{0 < @var{y}(1) < @dots{} < @var{x}(N)}.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("chebspec", @var{n})
## @deftypefnx {Function File} {@var{c} =} gallery ("chebspec", @var{n}, @var{k})
## Create a Chebyshev spectral differentiation matrix.
##
## Returns the Chebyshev spectral differentiation matrix of order @var{n}.
## The optional argument @var{k} can be:
##
## @table @asis
## @item 0 (default)
## No boundary conditions. @var{c} is nilpotent, with
## @code{@var{c}^@var{n}N = 0} and it has the null vector
## @code{ones (@var{n},1)}.
## @var{c} is nonsingular and well-conditioned, and its eigenvalues
## have negative real parts.
## @end table
##
## The computed eigenvector matrix @var{x} from @code{eig} is
## ill-conditioned (@code{mesh (real (@var{x}))} is interesting).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("chebvand", @var{p})
## @deftypefnx {Function File} {@var{c} =} gallery ("chebvand", @var{m}, @var{p})
## Create a Vandermonde-like matrix for the Chebyshev polynomials.
##
## If @var{p} is a vector, produces the (primal) Chebyshev Vandermonde
## matrix based on the points @var{p}, i.e.,
## @var{c}(i,j) = @math{T_{i-1}}(@var{p}(j)), where @math{T_{i-1}} is
## the Chebyshev polynomial of degree i-1.  If @var{p} is a scalar, then
## @var{p} equally spaced points on [0,1] are used.
##
## If the optional argument @var{m} is specified, returns a rectangular
## version of the matrix with @var{m} rows.
## @end deftypefn

## Code for the individual matrices by
## Nicholas .J. Higham <Nicholas.J.Higham@manchester.ac.uk>
## Adapted for Octave and into single gallery function by Carnë Draug

function [matrix, varargout] = gallery (name, varargin)

  if (nargin < 1)
    print_usage ();
  elseif (! ischar (name))
    error ("gallery: NAME must be a string.");
  endif

  switch (tolower (name))
    case "binomial"
      error ("gallery: matrix %s not implemented.", name);
    case "cauchy",      matrix = cauchy      (varargin{:});
    case "chebspec",    matrix = chebspec    (varargin{:});
    case "chebvand",    matrix = chebvand    (varargin{:});
    case "chow"
      error ("gallery: matrix %s not implemented.", name);
    case "circul"
      error ("gallery: matrix %s not implemented.", name);
    case "clement"
      error ("gallery: matrix %s not implemented.", name);
    case "compar"
      error ("gallery: matrix %s not implemented.", name);
    case "condex"
      error ("gallery: matrix %s not implemented.", name);
    case "cycol"
      error ("gallery: matrix %s not implemented.", name);
    case "dorr"
      error ("gallery: matrix %s not implemented.", name);
    case "dramadah"
      error ("gallery: matrix %s not implemented.", name);
    case "fiedler"
      error ("gallery: matrix %s not implemented.", name);
    case "forsythe"
      error ("gallery: matrix %s not implemented.", name);
    case "frank"
      error ("gallery: matrix %s not implemented.", name);
    case "gearmat"
      error ("gallery: matrix %s not implemented.", name);
    case "gcdmat"
      error ("gallery: matrix %s not implemented.", name);
    case "grcar"
      error ("gallery: matrix %s not implemented.", name);
    case "hanowa"
      error ("gallery: matrix %s not implemented.", name);
    case "house"
      error ("gallery: matrix %s not implemented.", name);
    case "integerdata"
      error ("gallery: matrix %s not implemented.", name);
    case "invhess"
      error ("gallery: matrix %s not implemented.", name);
    case "invol"
      error ("gallery: matrix %s not implemented.", name);
    case "ipjfact"
      error ("gallery: matrix %s not implemented.", name);
    case "jordbloc"
      error ("gallery: matrix %s not implemented.", name);
    case "kahan"
      error ("gallery: matrix %s not implemented.", name);
    case "kms"
      error ("gallery: matrix %s not implemented.", name);
    case "krylov"
      error ("gallery: matrix %s not implemented.", name);
    case "lauchli"
      error ("gallery: matrix %s not implemented.", name);
    case "lehmer"
      error ("gallery: matrix %s not implemented.", name);
    case "leslie"
      error ("gallery: matrix %s not implemented.", name);
    case "lesp"
      error ("gallery: matrix %s not implemented.", name);
    case "lotkin"
      error ("gallery: matrix %s not implemented.", name);
    case "minij"
      error ("gallery: matrix %s not implemented.", name);
    case "moler"
      error ("gallery: matrix %s not implemented.", name);
    case "neumann"
      error ("gallery: matrix %s not implemented.", name);
    case "normaldata"
      error ("gallery: matrix %s not implemented.", name);
    case "orthog"
      error ("gallery: matrix %s not implemented.", name);
    case "parter"
      error ("gallery: matrix %s not implemented.", name);
    case "pei"
      error ("gallery: matrix %s not implemented.", name);
    case "poisson"
      error ("gallery: matrix %s not implemented.", name);
    case "prolate"
      error ("gallery: matrix %s not implemented.", name);
    case "randcolu"
      error ("gallery: matrix %s not implemented.", name);
    case "randcorr"
      error ("gallery: matrix %s not implemented.", name);
    case "randhess"
      error ("gallery: matrix %s not implemented.", name);
    case "randjorth"
      error ("gallery: matrix %s not implemented.", name);
    case "rando"
      error ("gallery: matrix %s not implemented.", name);
    case "randsvd"
      error ("gallery: matrix %s not implemented.", name);
    case "redheff"
      error ("gallery: matrix %s not implemented.", name);
    case "riemann"
      error ("gallery: matrix %s not implemented.", name);
    case "ris"
      error ("gallery: matrix %s not implemented.", name);
    case "sampling"
      error ("gallery: matrix %s not implemented.", name);
    case "smoke"
      error ("gallery: matrix %s not implemented.", name);
    case "toeppd"
      error ("gallery: matrix %s not implemented.", name);
    case "tridiag"
      error ("gallery: matrix %s not implemented.", name);
    case "triw"
      error ("gallery: matrix %s not implemented.", name);
    case "uniformdata"
      error ("gallery: matrix %s not implemented.", name);
    case "wathen"
      error ("gallery: matrix %s not implemented.", name);
    case "wilk"
      error ("gallery: matrix %s not implemented.", name);
    otherwise
      error ("gallery: unknown matrix with NAME %s", name);
  endswitch

endfunction

function C = cauchy (x, y)
  ##CAUCHY  Cauchy matrix.
  ##        C = CAUCHY(X, Y), where X, Y are N-vectors, is the N-by-N matrix
  ##        with C(i,j) = 1/(X(i)+Y(j)).   By default, Y = X.
  ##        Special case: if X is a scalar CAUCHY(X) is the same as CAUCHY(1:X).
  ##        Explicit formulas are known for DET(C) (which is nonzero if X and Y
  ##        both have distinct elements) and the elements of INV(C).
  ##        C is totally positive if 0 < X(1) < ... < X(N) and
  ##        0 < Y(1) < ... < Y(N).
  ##
  ##        References:
  ##        N.J. Higham, Accuracy and Stability of Numerical Algorithms,
  ##          Society for Industrial and Applied Mathematics, Philadelphia, PA,
  ##          USA, 1996; sec. 26.1.
  ##        D.E. Knuth, The Art of Computer Programming, Volume 1,
  ##          Fundamental Algorithms, second edition, Addison-Wesley, Reading,
  ##          Massachusetts, 1973, p. 36.
  ##        E.E. Tyrtyshnikov, Cauchy-Toeplitz matrices and some applications,
  ##          Linear Algebra and Appl., 149 (1991), pp. 1-18.
  ##          O. Taussky and M. Marcus, Eigenvalues of finite matrices, in
  ##          Survey of Numerical Analysis, J. Todd, ed., McGraw-Hill, New York,
  ##          pp. 279-313, 1962. (States the totally positive property on p. 295.)

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for Cauchy matrix.");
  elseif (! isnumeric (x))
    error ("gallery: X must be numeric for Cauchy matrix.");
  elseif (nargin == 2 && ! isnumeric (y))
    error ("gallery: Y must be numeric for Cauchy matrix.");
  endif

  n = length (x);
  if (n == 1)
     n = x;
     x = 1:n;
  endif

  if (nargin == 1)
    y = x;
  endif

  ## Ensure x and y are column vectors
  x = x(:);
  y = y(:);
  if (numel (x) != numel (y))
    error ("gallery: X and Y must be vectors of same length for cauchy matrix.");
  endif

  C = x * ones (1, n) + ones (n, 1) * y.';
  C = ones (n) ./ C;
endfunction

function C = chebspec (n, k = 0)
  ## CHEBSPEC  Chebyshev spectral differentiation matrix.
  ##           C = CHEBSPEC(N, K) is a Chebyshev spectral differentiation
  ##           matrix of order N.  K = 0 (the default) or 1.
  ##           For K = 0 (`no boundary conditions'), C is nilpotent, with
  ##               C^N = 0 and it has the null vector ONES(N,1).
  ##               C is similar to a Jordan block of size N with eigenvalue zero.
  ##           For K = 1, C is nonsingular and well-conditioned, and its eigenvalues
  ##               have negative real parts.
  ##           For both K, the computed eigenvector matrix X from EIG is
  ##               ill-conditioned (MESH(REAL(X)) is interesting).
  ##
  ##           References:
  ##           C. Canuto, M.Y. Hussaini, A. Quarteroni and T.A. Zang, Spectral
  ##              Methods in Fluid Dynamics, Springer-Verlag, Berlin, 1988; p. 69.
  ##           L.N. Trefethen and M.R. Trummer, An instability phenomenon in
  ##              spectral methods, SIAM J. Numer. Anal., 24 (1987), pp. 1008-1023.
  ##           D. Funaro, Computing the inverse of the Chebyshev collocation
  ##              derivative, SIAM J. Sci. Stat. Comput., 9 (1988), pp. 1050-1057.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for Chebyshev matrix.");
  elseif (! isnumeric (n) || ! isscalar (n))
    error ("gallery: N must be a numeric scalar.");
  endif

  ## No check for k. The original code did not make this check and
  ## apparently matlab did not add it either. So while undocumented,
  ## for matlab compatibility K is 0 for ANYTHING other than a value of 1

  ## k = 1 case obtained from k = 0 case with one bigger n.
  if (k == 1)
    n = n + 1;
  endif

  n = n-1;
  C = zeros (n+1);

  one    = ones (n+1, 1);
  x      = cos ((0:n)' * (pi/n));
  d      = ones (n+1, 1);
  d(1)   = 2;
  d(n+1) = 2;

  ## eye(size(C)) on next line avoids div by zero.
  C = (d * (one./d)') ./ (x*one'-one*x' + eye (size (C)));

  ##  Now fix diagonal and signs.
  C(1,1) = (2*n^2+1)/6;
  for i = 2:n+1
      if (rem (i, 2) == 0)
         C(:,i) = -C(:,i);
         C(i,:) = -C(i,:);
      endif
      if (i < n+1)
         C(i,i) = -x(i)/(2*(1-x(i)^2));
      else
         C(n+1,n+1) = -C(1,1);
      endif
  endfor

  if (k == 1)
     C = C(2:n+1,2:n+1);
  endif
endfunction

function C = chebvand (m, p)
  ## CHEBVAND Vandermonde-like matrix for the Chebyshev polynomials.
  ##          C = CHEBVAND(P), where P is a vector, produces the (primal)
  ##          Chebyshev Vandermonde matrix based on the points P,
  ##          i.e., C(i,j) = T_{i-1}(P(j)), where T_{i-1} is the Chebyshev
  ##          polynomial of degree i-1.
  ##          CHEBVAND(M,P) is a rectangular version of CHEBVAND(P) with M rows.
  ##          Special case: If P is a scalar then P equally spaced points on
  ##                        [0,1] are used.
  ##
  ##          Reference:
  ##          N.J. Higham, Stability analysis of algorithms for solving confluent
  ##            Vandermonde-like systems, SIAM J. Matrix Anal. Appl., 11 (1990),
  ##            pp. 23-41.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for Chebvand matrix.");
  endif

  if (nargin == 1)
    p = m;
  endif
  n = length (p);

  if (n == 1)
     n = p;
     p = seqa (0, 1, n);
  endif

  if (nargin == 1)
    m = n;
  endif

  p = p(:).'; # Ensure p is a row vector.
  C = ones (m, n);
  if (m != 1)
    C(2,:) = p;
    ##      Use Chebyshev polynomial recurrence.
    for i = 3:m
        C(i,:) = 2.*p.*C(i-1,:) - C(i-2,:);
    endfor
  endif
endfunction


## subfunction used by the some of the matrix generating functions
function y = seqa (a, b, n = 10)
  ## SEQA   Additive sequence.
  ##        Y = SEQA(A, B, N) produces a row vector comprising N equally
  ##        spaced numbers starting at A and finishing at B.
  ##        If N is omitted then 10 points are generated.
  if (n <= 1)
     y = a;
  else
    y = [a+(0:n-2)*(b-a)/(n-1), b];
  endif
endfunction
