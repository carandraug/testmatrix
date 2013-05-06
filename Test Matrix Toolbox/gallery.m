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
    case "chebspec"
      error ("gallery: matrix %s not implemented.", name);
    case "chebvand"
      error ("gallery: matrix %s not implemented.", name);
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
