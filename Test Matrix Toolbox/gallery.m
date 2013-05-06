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
    case "cauchy"
      error ("gallery: matrix %s not implemented.", name);
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
