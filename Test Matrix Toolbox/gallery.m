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
##
## @deftypefn  {Function File} {@var{c} =} gallery ("chow", @var{n})
## @deftypefnx {Function File} {@var{c} =} gallery ("chow", @var{x}, @var{alpha})
## @deftypefnx {Function File} {@var{c} =} gallery ("chow", @var{x}, @var{alpha}, @var{delta})
## Create a Chow matrix -- a singular Toeplitz lower Hessenberg matrix.
##
## The matrix @var{c} is a Toeplitz lower Hessenberg matrix equal to
## @code{H(@var{alpha}) + @var{delta} * eye}, where
## @math{H(i,j) = @var{alpha}^(i-j+1)}.  H(@var{alpha}) has
## @code{p = floor (N/2)} zero eigenvalues, the rest being 
## @code{4*@var{alpha}*cos( k*pi/(N+2) )^2}, for @code{k = 1:N-p}.
##
## The optional arguments @var{alpha} and @var{delta} default to 1 and 0,
## respectively.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("circul", @var{v})
## Create a circulant matrix.
##
## @var{c} is circulant matrix whose first row is the vector @var{v}.
## If @var{v} is a scalar, it is interpreted as @code{1:@var{v}}.
##
## A circulant matrix has the property that each row is obtained
## from the previous one by cyclically permuting the entries one step
## forward; it is a special Toeplitz matrix in which the diagonals
## `wrap round'.
##
## The eigensystem of @var{c} (N-by-N) is known explicitly.  If t is an Nth
## root of unity, then the inner product of @var{v} with
## @code{W = [1 t t^2 @dots{} t^N]} is an eigenvalue of @var{c}, and
## @code{W(N:-1:1)} is an eigenvector of @var{c}.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("clement", @var{n})
## @deftypefnx {Function File} {@var{c} =} gallery ("clement", @var{n}, @var{k})
## Create a tridiagonal matrix with zero diagonal entries.
##
## @var{c} is a tridiagonal matrix with zero diagonal entries
## and known eigenvalues.  It is singular if @var{n} is odd.  About 64
## percent of the entries of the inverse are zero.  The eigenvalues
## are plus and minus the numbers N-1, N-3, N-5, @dots{}, (1 or 0).
##
## For @var{k} = 0 (default) the matrix is unsymmetric, while for
## @var{k} = 1 it is symmetric. @code{gallery ("clement", @var{n}, 1)} is
## diagonally similar to @code{gallery ("clement", @var{n})}.
##
## Similar properties hold for @code{gallery ("tridiag", X, Y, Z)] where
## @code{y = zeros (@var{n}, 1)}.  The eigenvalues still come in plus/minus
## pairs but they are not known explicitly
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("compar", @var{a})
## @deftypefnx {Function File} {@var{c} =} gallery ("compar", @var{a}, 1)
## Create a comparison matrix.
##
## @var{c} is @code{diag(@var{b}) - tril (@var{b}, -1) - triu (@var{b} ,1),
## where @code{@var{b} = abs (@var{a})}.  This is often denoted by M(@var{a})
## in the literature.
##
## If the second argument is 1, @var{c} is @var{a} with each diagonal
## element replaced by its absolute value, and each off-diagonal element
## replaced by minus the absolute value of the largest element in absolute
## value in its row.  However, if @var{a} is triangular @var{c} will be too.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("condex", @var{n})
## @deftypefnx {Function File} {@var{c} =} gallery ("condex", @var{n}, @var{k})
## @deftypefnx {Function File} {@var{c} =} gallery ("condex", @var{n}, @var{k}, @var{theta})
## Create a `counterexample' matrix to a condition estimator.
##
## @var{c} is a `conterexample' matrix of order @var{n}, and scalar parameter
## @var{theta} (default to 100).  If @var{n} is not equal to the `natural'
## size of the matrix then the matrix is padded out with an identity matrix
## to order @var{n}.  The matrix, its natural size, and the estimator to which
## it applies are specified by @var{K} as follows:
##
## @table @asis
## @item 1
## 4-by-4,     LINPACK (RCOND)
## @item 2
## 3-by-3,     LINPACK (RCOND)
## @item 3
## arbitrary,  LINPACK (RCOND) (independent of @var{theta})
## @item 4 (default)
## N >= 4,     SONEST (Higham 1988)
## @end table
##
## Note that in practice the K = 4 matrix is not usually a
## counterexample because of the rounding errors in forming it.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("cycol", [@var{m} @var{n}])
## @deftypefnx {Function File} {@var{c} =} gallery ("cycol", @var{n})
## @deftypefnx {Function File} {@var{c} =} gallery (@dots{}, @var{k})
## Create a matrix whose columns repeat cyclically.
##
## @var{c} is an @var{m}-by-@var{n} matrix of the form
## @var{c} = B(1:@var{m},1:@var{n}), where B = [A A A @dots{}] and
## @code{A = randn (@var{m}, @var{k})}.  Thus @var{c}'s columns repeat
## cyclically, and @var{c} has rank at most @var{K}.   @var{k} need not
## divide @var{n}.  @var{k} defaults to @code{round (@var{n}/4)}.  If
## @var{n} is a scalar, @var{c} is a @var{n}-by-@var{n} matrix.
##
## This type of matrix can lead to underflow problems for Gaussian
## elimination: see NA Digest Volume 89, Issue 3 (January 22, 1989).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("triw", @var{n})
## @deftypefnx {Function File} {@var{c} =} gallery ("triw", @var{n}, @var{alpha})
## @deftypefnx {Function File} {@var{c} =} gallery ("triw", @var{n}, @var{alpha}, @var{k})
## Create an upper triangular matrix discussed by Kahan, Golub and Wilkinson.
##
## @var{c} is the upper triangular matrix with ones on the diagonal
## and @var{alpha}s on the first @var{k} >= 0 superdiagonals.  Both @var{alpha}
## and @var{k} default to -1 (full upper triangle).
##
## @var{n} may be a 2 element vector in which case the matrix is
## @var{n}(1)-by-@var{n}(2) and upper trapezoidal.
##
## Ostrowski (1954) shows that
## @code{cond (gallery ("triw", @var{n}, 2)) = cot (pi/(4*@var{n}))^2}
## and for large @code{abs (@var{alpha})},
## @code{cond (gallery ("triw", @var{n}, @var{alpha}))} is approximately
## @code{abs (@var{alpha})^@var{n}*sin (pi/(4 * @var{n}-2)).
##
## Adding @code{-2^(2-@var{n})} to the (@var{n},1) element makes
## @var{c} singular, as does adding @code{-2^(1-@var{n})} to all elements
## in the first column.
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

  ## NOTE: there isn't a lot of input check in the individual functions
  ## that actually build the functions.  This is by design. The original
  ## code by Higham did not perform it and was propagated to Matlab, so
  ## for compatibility, we also don't make it. For example, arguments
  ## that behave as switches, and in theory accepting a value of 0 or 1,
  ## will use a value of 0, for any value other than 1 (only check made
  ## is if the value is equal to 1). It will often also accept string
  ## values instead of numeric. Only input check added was where it
  ## would be causing an error anyway.

  switch (tolower (name))
    case "binomial"
      error ("gallery: matrix %s not implemented.", name);
    case "cauchy",      matrix = cauchy      (varargin{:});
    case "chebspec",    matrix = chebspec    (varargin{:});
    case "chebvand",    matrix = chebvand    (varargin{:});
    case "chow",        matrix = chow        (varargin{:});
    case "circul",      matrix = circul      (varargin{:});
    case "clement",     matrix = clement     (varargin{:});
    case "compar",      matrix = compar      (varargin{:});
    case "condex",      matrix = condex      (varargin{:});
    case "cycol",       matrix = cycol       (varargin{:});
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
    case "triw",        matrix = triw        (varargin{:});
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

function A = chow (n, alpha = 1, delta = 0)
  ## CHOW    Chow matrix - a singular Toeplitz lower Hessenberg matrix.
  ##         A = CHOW(N, ALPHA, DELTA) is a Toeplitz lower Hessenberg matrix
  ##         A = H(ALPHA) + DELTA*EYE, where H(i,j) = ALPHA^(i-j+1).
  ##         H(ALPHA) has p = FLOOR(N/2) zero eigenvalues, the rest being
  ##         4*ALPHA*COS( k*PI/(N+2) )^2, k=1:N-p.
  ##         Defaults: ALPHA = 1, DELTA = 0.
  ##
  ##         References:
  ##         T.S. Chow, A class of Hessenberg matrices with known
  ##            eigenvalues and inverses, SIAM Review, 11 (1969), pp. 391-395.
  ##         G. Fairweather, On the eigenvalues and eigenvectors of a class of
  ##            Hessenberg matrices, SIAM Review, 13 (1971), pp. 220-221.

  if (nargin < 1 || nargin > 3)
    error ("gallery: 1 to 3 arguments are required for Chow matrix.");
  elseif (! isnumeric (n) || ! isscalar (n))
    error ("gallery: N must be a numeric scalar");
  elseif (! isnumeric (alpha) || ! isscalar (alpha))
    error ("gallery: ALPHA must be a numeric scalar");
  elseif (! isnumeric (delta) || ! isscalar (delta))
    error ("gallery: DELTA must be a numeric scalar");
  endif

  A = toeplitz (alpha.^(1:n), [alpha 1 zeros(1, n-2)]) + delta * eye (n);
endfunction

function C = circul (v)
  ## CIRCUL  Circulant matrix.
  ##         C = CIRCUL(V) is the circulant matrix whose first row is V.
  ##         (A circulant matrix has the property that each row is obtained
  ##         from the previous one by cyclically permuting the entries one step
  ##         forward; it is a special Toeplitz matrix in which the diagonals
  ##         `wrap round'.)
  ##         Special case: if V is a scalar then C = CIRCUL(1:V).
  ##         The eigensystem of C (N-by-N) is known explicitly.   If t is an Nth
  ##         root of unity, then the inner product of V with W = [1 t t^2 ... t^N]
  ##         is an eigenvalue of C, and W(N:-1:1) is an eigenvector of C.
  ##
  ##         Reference:
  ##         P.J. Davis, Circulant Matrices, John Wiley, 1977.

  if (nargin != 1)
    error ("gallery: only 1 argument required for circul matrix.");
  endif

  n = length (v);
  if n == 1
     n = v;
     v = 1:n;
  endif

  v = v(:).';   # Make sure v is a row vector
  C = toeplitz ([v(1) v(n:-1:2)], v);
endfunction

function A = clement (n, k = 0)
  ## CLEMENT   Clement matrix - tridiagonal with zero diagonal entries.
  ##           CLEMENT(N, K) is a tridiagonal matrix with zero diagonal entries
  ##           and known eigenvalues.  It is singular if N is odd.  About 64
  ##           percent of the entries of the inverse are zero.  The eigenvalues
  ##           are plus and minus the numbers N-1, N-3, N-5, ..., (1 or 0).
  ##           For K = 0 (the default) the matrix is unsymmetric, while for
  ##           K = 1 it is symmetric.
  ##           CLEMENT(N, 1) is diagonally similar to CLEMENT(N).
  ##
  ##           Similar properties hold for TRIDIAG(X,Y,Z) where Y = ZEROS(N,1).
  ##           The eigenvalues still come in plus/minus pairs but they are not
  ##           known explicitly.
  ##
  ##           References:
  ##           P.A. Clement, A class of triple-diagonal matrices for test
  ##              purposes, SIAM Review, 1 (1959), pp. 50-52.
  ##           A. Edelman and E. Kostlan, The road from Kac's matrix to Kac's
  ##              random polynomials. In John~G. Lewis, editor, Proceedings of
  ##              the Fifth SIAM Conference on Applied Linear Algebra Society
  ##              for Industrial and Applied Mathematics, Philadelphia, 1994,
  ##              pp. 503-507.
  ##           O. Taussky and J. Todd, Another look at a matrix of Mark Kac,
  ##              Linear Algebra and Appl., 150 (1991), pp. 341-360.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for clement matrix.");
  endif

  n = n-1;
  x = n:-1:1;
  z = 1:n;

  if (k == 0)
     A = diag (x, -1) + diag (z, 1);
  else
     y = sqrt (x.*z);
     A = diag (y, -1) + diag (y, 1);
  endif
endfunction

function C = compar (A, k = 0)
  ## COMP    Comparison matrices.
  ##         COMP(A) is DIAG(B) - TRIL(B,-1) - TRIU(B,1), where B = ABS(A).
  ##         COMP(A, 1) is A with each diagonal element replaced by its
  ##         absolute value, and each off-diagonal element replaced by minus
  ##         the absolute value of the largest element in absolute value in
  ##         its row.  However, if A is triangular COMP(A, 1) is too.
  ##         COMP(A, 0) is the same as COMP(A).
  ##         COMP(A) is often denoted by M(A) in the literature.
  ##
  ##         Reference (e.g.):
  ##         N.J. Higham, A survey of condition number estimation for
  ##         triangular matrices, SIAM Review, 29 (1987), pp. 575-596.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for compar matrix.");
  endif

  [m, n] = size (A);
  p = min (m, n);

  if (k == 0)
    ## This code uses less temporary storage than
    ## the `high level' definition above.
    C = -abs (A);
    for j = 1:p
      C(j,j) = abs (A(j,j));
    end

  elseif (k == 1)
    C = A';
    for j = 1:p
      C(k,k) = 0;
    endfor
    mx = max (abs (C));
    C  = -mx'*ones (1, n);
    for j = 1:p
      C(j,j) = abs (A(j,j));
    endfor
    if (all (A == tril (A))), C = tril(C); endif
    if (all (A == triu (A))), C = triu(C); endif

  else
    error ("gallery: K must be 0 or 1 for compar matrix.");
  endif

endfunction

function A = condex (n, k = 4, theta = 100)
  ## CONDEX   `Counterexamples' to matrix condition number estimators.
  ##          CONDEX(N, K, THETA) is a `counterexample' matrix to a condition
  ##          estimator.  It has order N and scalar parameter THETA (default 100).
  ##          If N is not equal to the `natural' size of the matrix then
  ##          the matrix is padded out with an identity matrix to order N.
  ##          The matrix, its natural size, and the estimator to which it applies
  ##          are specified by K (default K = 4) as follows:
  ##              K = 1:   4-by-4,     LINPACK (RCOND)
  ##              K = 2:   3-by-3,     LINPACK (RCOND)
  ##              K = 3:   arbitrary,  LINPACK (RCOND) (independent of THETA)
  ##              K = 4:   N >= 4,     SONEST (Higham 1988)
  ##          (Note that in practice the K = 4 matrix is not usually a
  ##           counterexample because of the rounding errors in forming it.)
  ##
  ##          References:
  ##          A.K. Cline and R.K. Rew, A set of counter-examples to three
  ##             condition number estimators, SIAM J. Sci. Stat. Comput.,
  ##             4 (1983), pp. 602-611.
  ##          N.J. Higham, FORTRAN codes for estimating the one-norm of a real or
  ##             complex matrix, with applications to condition estimation
  ##             (Algorithm 674), ACM Trans. Math. Soft., 14 (1988), pp. 381-396.

  if (nargin < 1 || nargin > 3)
    error ("gallery: 1 to 3 arguments are required for condex matrix.");
  elseif (! isnumeric (n))
    error ("gallery: N must be numeric for condex matrix.");
  elseif (! isnumeric (k) || ! isscalar (k))
    error ("gallery: K must be a numeric scalar for condex matrix.);
  endif

  if (k == 1)       # Cline and Rew (1983), Example B.
     A = [1  -1  -2*theta     0
          0   1     theta  -theta
          0   1   1+theta  -(theta+1)
          0   0   0         theta];

  elseif (k == 2)   # Cline and Rew (1983), Example C.
    A = [1   1-2/theta^2  -2
         0   1/theta      -1/theta
         0   0             1];

  elseif (k == 3)   # Cline and Rew (1983), Example D.
    A = gallery ("triw", n, -1)';
    A(n,n) = -1;

  elseif (k == 4)   # Higham (1988), p. 390.
    x = ones (n, 3);            #  First col is e
    x(2:n,2) = zeros (n-1, 1);  #  Second col is e(1)

    ## Third col is special vector b in SONEST
    x(:, 3) = (-1).^[0:n-1]' .* ( 1 + [0:n-1]'/(n-1) );

    Q = orth (x);  #  Q*Q' is now the orthogonal projector onto span(e(1),e,b)).
    P = eye (n) - Q*Q';
    A = eye (n) + theta*P;

  else
    error ("gallery: unknown estimator K for condex matrix.");
  endif

  ## Pad out with identity as necessary.
  [m, m] = size (A);
  if m < n
    for i = n:-1:m+1
      A(i,i) = 1;
    endfor
  endif
endfunction

function A = cycol (n, k)
  ## CYCOL   Matrix whose columns repeat cyclically.
  ##         A = CYCOL([M N], K) is an M-by-N matrix of the form A = B(1:M,1:N)
  ##         where B = [C C C...] and C = RANDN(M, K).  Thus A's columns repeat
  ##         cyclically, and A has rank at most K.   K need not divide N.
  ##         K defaults to ROUND(N/4).
  ##         CYCOL(N, K), where N is a scalar, is the same as CYCOL([N N], K).
  ##
  ##         This type of matrix can lead to underflow problems for Gaussian
  ##         elimination: see NA Digest Volume 89, Issue 3 (January 22, 1989).

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for cycol matrix.");
  elseif (! isscalar (k))
    error ("gallery: K must be a scalar for cycol matrix.);
  endif

  m = n(1);              # Parameter n specifies dimension: m-by-n.
  n = n(length (n));

  if (nargin < 2)
    k = max (round (n/4), 1);
  endif

  A = randn (m, k);
  for i = 2:ceil (n/k)
    A = [A A(:, 1:k)];
  endfor
  A = A(:, 1:n);
endfunction

function t = triw (n, alpha = -1, k = -1)
  ## TRIW   Upper triangular matrix discussed by Wilkinson and others.
  ##        TRIW(N, ALPHA, K) is the upper triangular matrix with ones on
  ##        the diagonal and ALPHAs on the first K >= 0 superdiagonals.
  ##        N may be a 2-vector, in which case the matrix is N(1)-by-N(2) and
  ##        upper trapezoidal.
  ##        Defaults: ALPHA = -1,
  ##                  K = N - 1     (full upper triangle).
  ##        TRIW(N) is a matrix discussed by Kahan, Golub and Wilkinson.
  ##
  ##        Ostrowski (1954) shows that
  ##          COND(TRIW(N,2)) = COT(PI/(4*N))^2,
  ##        and for large ABS(ALPHA),
  ##          COND(TRIW(N,ALPHA)) is approximately ABS(ALPHA)^N*SIN(PI/(4*N-2)).
  ##
  ##        Adding -2^(2-N) to the (N,1) element makes TRIW(N) singular,
  ##        as does adding -2^(1-N) to all elements in the first column.
  ##
  ##        References:
  ##        G.H. Golub and J.H. Wilkinson, Ill-conditioned eigensystems and the
  ##           computation of the Jordan canonical form, SIAM Review,
  ##           18(4), 1976, pp. 578-619.
  ##        W. Kahan, Numerical linear algebra, Canadian Math. Bulletin,
  ##           9 (1966), pp. 757-801.
  ##        A.M. Ostrowski, On the spectrum of a one-parametric family of
  ##           matrices, J. Reine Angew. Math., 193 (3/4), 1954, pp. 143-160.
  ##        J.H. Wilkinson, Singular-value decomposition---basic aspects,
  ##           in D.A.H. Jacobs, ed., Numerical Software---Needs and Availability,
  ##           Academic Press, London, 1978, pp. 109-135.

  if (nargin < 1 || nargin > 3)
    error ("gallery: 1 to 3 arguments are required for triw matrix.");
  elseif (! isscalar (alpha))
     error("gallery: ALPHA must be a scalar for triw matrix.")
  endif

  m = n(1);              # Parameter n specifies dimension: m-by-n.
  n = n(length (n));

  t = tril (eye (m, n) + alpha * triu (ones (m, n), 1), k);
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
