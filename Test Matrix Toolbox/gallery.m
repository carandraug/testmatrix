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
## Similar properties hold for @code{gallery ("tridiag", X, Y, Z)]} where
## @code{y = zeros (@var{n}, 1)}.  The eigenvalues still come in plus/minus
## pairs but they are not known explicitly
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("compar", @var{a})
## @deftypefnx {Function File} {@var{c} =} gallery ("compar", @var{a}, 1)
## Create a comparison matrix.
##
## @var{c} is @code{diag(@var{b}) - tril (@var{b}, -1) - triu (@var{b} ,1)},
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
## @deftypefn  {Function File} {[@var{c},@var{d}, @var{e}] =} gallery ("dorr", @var{n})
## @deftypefnx {Function File} {[@var{c},@var{d}, @var{e}] =} gallery ("dorr", @var{n}, @var{theta})
## @deftypefnx {Function File} {@var{a} =} gallery ("dorr", @dots{})
## Create a diagonally dominant, ill conditioned, tridiagonal matrix.
##
## @var{c}, @var{d}, and @var{e} are the vectors defining a row diagonally
## dominant, tridiagonal M-matrix that is ill conditioned for small
## values of the parameter @var{theta} >= 0 (defaults is 0.01).
##
## Alternatively, @var{a} is
## @code{full (gallery ("tridiag", @var{c}, @var{d}, @var{e}))}, i.e., the
## matrix iself is returned.
##
## The columns of @code{inv (@var{c})} vary greatly in norm.
## The amount of diagonal dominance is given by (ignoring rounding errors):
## @code{gallery ("compar", @var{c}) * ones (@var{n},1) = @var{theta}*(@var{n}+1)^2 * [1 0 0 @dots{} 0 1]'}.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("dramadah", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("dramadah", @var{n}, @var{k})
## Create a (0, 1) matrix whose inverse has large integer entries.
##
## An anti-Hadamard matrix @var{a} is a matrix with elements 0 or 1 for
## which MU(A) := norm (inv (A),"fro") is maximal.
## The matrix @var{a} is an @var{n}-by-@var{n} (0,1) matrix for which
## MU(@var{a}) is relatively large, although not necessarily maximal.
##
## The optional value @var{k}, defines the type:
## @table @asis
## @item 1 (default)
## @var{a} is Toeplitz, with @code{abs (det (@var{a}))} equal to 1,
## and MU(@var{a}) > c(1.75)^@var{n}, where c is a constant.
## @item 2
## @var{a} is upper triangular and Toeplitz.
## The inverses of both types have integer entries.
## @item 3
## @var{a} has maximal determinant among (0,1) lower Hessenberg
## matrices: @code{det (@var{a})} equal to the n'th Fibonacci number.
## @var{a} is Toeplitz.  The eigenvalues have an interesting distribution
## in the complex plane.
## @end table
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("fiedler", @var{c})
## Create a symmetric Fiedler matrix.
##
## @var{a} is a N-by-N symmetric matrix where N is the length of
## the vector @var{c}, and with elements @code{abs (@var{c}(i) - @var{c}(j))}.
## If @var{c} is a scalar, the vector 1:@var{c} is used insted, (i.e.
## @code{@var{a}(i,j) = abs (i - j)}.
##
## Properties:
## @itemize @bullet
## @item
## @var{a} has a dominant positive eigenvalue and all the other
## eigenvalues are negative (Szego, 1936).
## @item
## Explicit formulas for @code{inv (@var{a})} and @var{det (@var{a})} are
## given by Todd (1977) and attributed to Fiedler.  These indicate that
## @code{inv (@var{a})} is tridiagonal except for nonzero (1,n) and
## (n,1) elements.
## @end itemize
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("forsythe", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("forsythe", @var{n}, @var{alpha})
## @deftypefnx {Function File} {@var{a} =} gallery ("forsythe", @var{n}, @var{alpha}, @var{lambda})
## Create a Forsythe matrix (a perturbed Jordan block).
##
## The matrix @var{a} is the @var{n}-by-@var{n} matrix equal to
## @code{gallery ("jordbloc", @var{n}, @var{lambda})} except it has
## an @var{alpha} in the (@var{n},1) position.
##
## It has the characteristic polynomial
## @code{det (@var{a}-t*eye) = (@var{lambda}-t)^@var{n} - (-1)^@var{n} @var{alpha}}.
##
## @var{alpha} defaults to @code{sqrt (eps)} and @var{lambda} to 0.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{f} =} gallery ("frank", @var{n})
## @deftypefnx {Function File} {@var{f} =} gallery ("frank", @var{n}, @var{k})
## Create a Frank matrix (ill conditioned eigenvalues).
##
## @var{f} is the Frank matrix of order @var{n}.  It is upper
## Hessenberg with determinant 1.  @var{k} = 0 is the default; if
## @var{k} = 1 the elements are reflected about the anti-diagonal
## (1,@var{n})--(@var{n},1).
##
## @var{f} has all positive eigenvalues and they occur in reciprocal pairs
## (so that 1 is an eigenvalue if @var{n} is odd).
## The eigenvalues of @var{f} may be obtained in terms of the zeros of the
## Hermite polynomials.
##
## The @code{floor (@var{n}/2)} smallest eigenvalues of @var{f} are ill
## conditioned, the more so for bigger @var{n}.
##
## @code{det (gallery ("frank", @var{n})')} comes out far from 1 for
## large @var{n} --- see Frank (1958) and Wilkinson (1960) for discussions.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{c} =} gallery ("gcdmat", @var{n})
## Create a greatest common divisor matrix.
##
## @var{c} is an @var{n}-by-@var{n} matrix whose values correspond to the
## greatest common divisor of its coordinate values, i.e., @var{c}(i,j)
## correspond @code{gcd (i, j)}.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("gearmat", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("gearmat", @var{n}, @var{i})
## @deftypefnx {Function File} {@var{a} =} gallery ("gearmat", @var{n}, @var{i}, @var{j})
## Create a Gear matrix.
##
## @var{a} is the @var{n}-by-@var{n} matrix with ones on the sub- and
## super-diagonals, @code{sign (@var{i})} in the @code{(1, abs (@var{i}))}
## position, @code{sign (@var{j})} in the @code{(@var{n}, @var{n}+1-abs (@var{j}))}
## position, and zeros everywhere else.
##
## @var{i} and @var{j} default to @var{n} and -@var{n}, respectively, in
## which case, @var{a} is singular.
##
## All eigenvalues are of the form @code{2*cos(a)} and the eigenvectors
## are of the form @code{[sin(w+a), sin(w+2a), @dots{}, sin(w+Na)]}.
## The values of a and w are given in C.W. Gear, A simple set of test
## matrices for eigenvalue programs, Math. Comp., 23 (1969), pp. 119-125.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{g} =} gallery ("grcar", @var{n})
## @deftypefnx {Function File} {@var{g} =} gallery ("grcar", @var{n}, @var{k})
## Create a Toeplitz matrix with sensitive eigenvalues.
##
## GRCAR(N, K) is an N-by-N matrix with -1s on the
## subdiagonal, 1s on the diagonal, and K superdiagonals of 1s.
## The default is K = 3.  The eigenvalues of this matrix form an
## interesting pattern in the complex plane (try PS(GRCAR(32))).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("hanowa", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("hanowa", @var{n}, @var{d})
## Create a matrix whose eigenvalues lie on a vertical line in the complex plane.
##
## HANOWA(N, d) is the N-by-N block 2x2 matrix (thus N = 2M must be even)
##               [d*EYE(M)   -DIAG(1:M)
##                DIAG(1:M)   d*EYE(M)]
## It has complex eigenvalues lambda(k) = d +/- k*i  (1 <= k <= M).
## Parameter d defaults to -1.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{v} =} gallery ("house", @var{x})
## @deftypefnx {Function File} {[@var{v}, @var{beta}] =} gallery ("house", @var{x})
## Create a householder matrix.
##
## If [v, beta] = HOUSE(x) then H = EYE - beta*v*v' is a Householder
## matrix such that Hx = -sign(x(1))*norm(x)*e_1.
## NB: If x = 0 then v = 0, beta = 1 is returned.
##     x can be real or complex.
##     sign(x) := exp(i*arg(x)) ( = x./abs(x) when x ~= 0).
##
## Theory: (textbook references Golub & Van Loan 1989, 38-43;
##          Stewart 1973, 231-234, 262; Wilkinson 1965, 48-50).
## Hx = y: (I - beta*v*v')x = -s*e_1.
## Must have |s| = norm(x), v = x+s*e_1, and
## x'y = x'Hx =(x'Hx)' real => arg(s) = arg(x(1)).
## So take s = sign(x(1))*norm(x) (which avoids cancellation).
## v'v = (x(1)+s)^2 + x(2)^2 + ... + x(n)^2
##     = 2*norm(x)*(norm(x) + |x(1)|).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("invhess", @var{x})
## @deftypefnx {Function File} {@var{a} =} gallery ("invhess", @var{x}, @var{y})
## Create the inverse of an upper Hessenberg matrix.
##
## INVHESS(X, Y), where X is an N-vector and Y an N-1 vector,
## is the matrix whose lower triangle agrees with that of
## ONES(N,1)*X' and whose strict upper triangle agrees with
## that of [1 Y]*ONES(1,N).
## The matrix is nonsingular if X(1) ~= 0 and X(i+1) ~= Y(i)
## for all i, and its inverse is an upper Hessenberg matrix.
## If Y is omitted it defaults to -X(1:N-1).
## Special case: if X is a scalar INVHESS(X) is the same as
## INVHESS(1:X).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("invol", @var{n})
## Create an involutory matrix.
##
## A = INVOL(N) is an N-by-N involutory (A*A = EYE(N)) and
## ill-conditioned matrix.
## It is a diagonally scaled version of HILB(N).
## NB: B = (EYE(N)-A)/2 and B = (EYE(N)+A)/2 are idempotent (B*B = B).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("ipjfact", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("ipjfact", @var{n}, @var{k})
## Create an Hankel matrix with factorial elements.
##
## A = IPJFACT(N, K) is the matrix with
##           A(i,j) = (i+j)!    (K = 0, default)
##           A(i,j) = 1/(i+j)!  (K = 1)
## Both are Hankel matrices.
## The determinant and inverse are known explicitly.
## If a second output argument is present, d = DET(A) is returned:
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("jordbloc", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("jordbloc", @var{n}, @var{lambda})
## Create a Jordan block.
##
## JORDBLOC(N, LAMBDA) is the N-by-N Jordan block with eigenvalue
## LAMBDA.  LAMBDA = 1 is the default.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{u} =} gallery ("kahan", @var{n})
## @deftypefnx {Function File} {@var{u} =} gallery ("kahan", @var{n}, @var{theta})
## @deftypefnx {Function File} {@var{u} =} gallery ("kahan", @var{n}, @var{theta}, @var{pert})
## Create a Kahan matrix (upper trapezoidal).
##
##        KAHAN(N, THETA) is an upper trapezoidal matrix
##        that has some interesting properties regarding estimation of
##        condition and rank.
##        The matrix is N-by-N unless N is a 2-vector, in which case it
##        is N(1)-by-N(2).
##        The parameter THETA defaults to 1.2.
##        The useful range of THETA is 0 < THETA < PI.
##
##        To ensure that the QR factorization with column pivoting does not
##        interchange columns in the presence of rounding errors, the diagonal
##        is perturbed by PERT*EPS*diag( [N:-1:1] ).
##        The default is PERT = 25, which ensures no interchanges for KAHAN(N)
##        up to at least N = 90 in IEEE arithmetic.
##        KAHAN(N, THETA, PERT) uses the given value of PERT.
##
##        The inverse of KAHAN(N, THETA) is known explicitly: see
##        Higham (1987, p. 588), for example.
##        The diagonal perturbation was suggested by Christian Bischof.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("kms", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("kms", @var{n}, @var{rho})
## Create a Kac-Murdock-Szego Toeplitz matrix.
##
##       A = KMS(N, RHO) is the N-by-N Kac-Murdock-Szego Toeplitz matrix with
##       A(i,j) = RHO^(ABS((i-j))) (for real RHO).
##       If RHO is complex, then the same formula holds except that elements
##       below the diagonal are conjugated.
##       RHO defaults to 0.5.
##       Properties:
##          A has an LDL' factorization with
##                   L = INV(TRIW(N,-RHO,1)'),
##                   D(i,i) = (1-ABS(RHO)^2)*EYE(N) except D(1,1) = 1.
##          A is positive definite if and only if 0 < ABS(RHO) < 1.
##          INV(A) is tridiagonal.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{b} =} gallery ("krylov", @var{a})
## @deftypefnx {Function File} {@var{b} =} gallery ("krylov", @var{a}, @var{x})
## @deftypefnx {Function File} {@var{b} =} gallery ("krylov", @var{a}, @var{x}, @var{j})
## Create a Krylov matrix.
##
##           KRYLOV(A, x, j) is the Krylov matrix
##                [x, Ax, A^2x, ..., A^(j-1)x],
##           where A is an n-by-n matrix and x is an n-vector.
##           Defaults: x = ONES(n,1), j = n.
##           KRYLOV(n) is the same as KRYLOV(RANDN(n)).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("lauchli", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("lauchli", @var{n}, @var{mu})
## Create a Lauchli matrix (rectangular).
##
##           LAUCHLI(N, MU) is the (N+1)-by-N matrix [ONES(1,N); MU*EYE(N))].
##           It is a well-known example in least squares and other problems
##           that indicates the dangers of forming A'*A.
##           MU defaults to SQRT(EPS).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("lehmer", @var{n})
## Create a Lehmer matrix (symmetric positive definite).
##
##         A = LEHMER(N) is the symmetric positive definite N-by-N matrix with
##                          A(i,j) = i/j for j >= i.
##         A is totally nonnegative.  INV(A) is tridiagonal, and explicit
##         formulas are known for its entries.
##         N <= COND(A) <= 4*N*N.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{t} =} gallery ("lesp", @var{nx})
## Create a tridiagonal matrix with real, sensitive eigenvalues.
##
##        LESP(N) is an N-by-N matrix whose eigenvalues are real and smoothly
##        distributed in the interval approximately [-2*N-3.5, -4.5].
##        The sensitivities of the eigenvalues increase exponentially as
##        the eigenvalues grow more negative.
##        The matrix is similar to the symmetric tridiagonal matrix with
##        the same diagonal entries and with off-diagonal entries 1,
##        via a similarity transformation with D = diag(1!,2!,...,N!).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("lotkin", @var{n})
## Create a Lotkin matrix.
##
##         A = LOTKIN(N) is the Hilbert matrix with its first row altered to
##         all ones.  A is unsymmetric, ill-conditioned, and has many negative
##         eigenvalues of small magnitude.
##         The inverse has integer entries and is known explicitly.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("minij", @var{n})
## Create a symmetric positive definite matrix MIN(i,j).
##
##         A = MINIJ(N) is the N-by-N symmetric positive definite matrix with
##         A(i,j) = MIN(i,j).
##         Properties, variations:
##         INV(A) is tridiagonal: it is minus the second difference matrix
##                     except its (N,N) element is 1.
##         2*A-ONES(N) (Givens' matrix) has tridiagonal inverse and
##                     eigenvalues .5*sec^2([2r-1)PI/4N], r=1:N.
##         (N+1)*ONES(N)-A also has a tridiagonal inverse.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("moler", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("moler", @var{n}, @var{alpha})
## Create a Moler matrix (symmetric positive definite).
##
##         A = MOLER(N, ALPHA) is the symmetric positive definite N-by-N matrix
##         U'*U where U = TRIW(N, ALPHA).
##         For ALPHA = -1 (the default) A(i,j) = MIN(i,j)-2, A(i,i) = i.
##         A has one small eigenvalue.
##
##         Nash (1990) attributes the ALPHA = -1 matrix to Moler.
## @end deftypefn
##
## @deftypefn  {Function File} {[@var{a}, @var{t}] =} gallery ("neumann", @var{n})
## Create a singular matrix from the discrete Neumann problem (sparse).
##
##          NEUMANN(N) is the singular, row diagonally dominant matrix resulting
##          from discretizing the Neumann problem with the usual five point
##          operator on a regular mesh.
##          It has a one-dimensional null space with null vector ONES(N,1).
##          The dimension N should be a perfect square, or else a 2-vector,
##          in which case the dimension of the matrix is N(1)*N(2).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{q} =} gallery ("orthog", @var{n})
## @deftypefnx {Function File} {@var{q} =} gallery ("orthog", @var{n}, @var{k})
## Create orthogonal and nearly orthogonal matrices.
##
##        Q = ORTHOG(N, K) selects the K'th type of matrix of order N.
##        K > 0 for exactly orthogonal matrices, K < 0 for diagonal scalings of
##        orthogonal matrices.
##        Available types: (K = 1 is the default)
##        K = 1:  Q(i,j) = SQRT(2/(n+1)) * SIN( i*j*PI/(n+1) )
##                Symmetric eigenvector matrix for second difference matrix.
##        K = 2:  Q(i,j) = 2/SQRT(2*n+1)) * SIN( 2*i*j*PI/(2*n+1) )
##                Symmetric.
##        K = 3:  Q(r,s) = EXP(2*PI*i*(r-1)*(s-1)/n) / SQRT(n)  (i=SQRT(-1))
##                Unitary, the Fourier matrix.  Q^4 is the identity.
##                This is essentially the same matrix as FFT(EYE(N))/SQRT(N)!
##        K = 4:  Helmert matrix: a permutation of a lower Hessenberg matrix,
##                whose first row is ONES(1:N)/SQRT(N).
##        K = 5:  Q(i,j) = SIN( 2*PI*(i-1)*(j-1)/n ) + COS( 2*PI*(i-1)*(j-1)/n ).
##                Symmetric matrix arising in the Hartley transform.
##        K = -1: Q(i,j) = COS( (i-1)*(j-1)*PI/(n-1) )
##                Chebyshev Vandermonde-like matrix, based on extrema of T(n-1).
##        K = -2: Q(i,j) = COS( (i-1)*(j-1/2)*PI/n) )
##                Chebyshev Vandermonde-like matrix, based on zeros of T(n).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("parter", @var{n})
## Create a Parter matrix (a Toeplitz matrix with singular values near pi).
##
##           PARTER(N) is the matrix with (i,j) element 1/(i-j+0.5).
##           It is a Cauchy matrix and a Toeplitz matrix.
##
##           At the Second SIAM Conference on Linear Algebra, Raleigh, N.C.,
##           1985, Cleve Moler noted that most of the singular values of
##           PARTER(N) are very close to PI.  An explanation of the phenomenon
##           was given by Parter; see also the paper by Tyrtyshnikov.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{p} =} gallery ("pei", @var{n})
## @deftypefnx {Function File} {@var{p} =} gallery ("pei", @var{n}, @var{alpha})
## Create a Pei matrix.
##
##        PEI(N, ALPHA), where ALPHA is a scalar, is the symmetric matrix
##        ALPHA*EYE(N) + ONES(N).
##        If ALPHA is omitted then ALPHA = 1 is used.
##        The matrix is singular for ALPHA = 0, -N.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("poisson", @var{n})
## Create a block tridiagonal matrix from Poisson's equation (sparse).
##
##           POISSON(N) is the block tridiagonal matrix of order N^2
##           resulting from discretizing Poisson's equation with the
##           5-point operator on an N-by-N mesh.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("prolate", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("prolate", @var{n}, @var{w})
## Create a prolate matrix (symmetric, ill-conditioned Toeplitz matrix).
##
##           A = PROLATE(N, W) is the N-by-N prolate matrix with parameter W.
##           It is a symmetric Toeplitz matrix.
##           If 0 < W < 0.5 then
##              - A is positive definite
##              - the eigenvalues of A are distinct, lie in (0, 1), and
##                tend to cluster around 0 and 1.
##           W defaults to 0.25.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{h} =} gallery ("randhess", @var{x})
## Create a random, orthogonal upper Hessenberg matrix.
##
##        H = RANDHESS(N) is an N-by-N real, random, orthogonal
##        upper Hessenberg matrix.
##        Alternatively, H = RANDHESS(X), where X is an arbitrary real
##        N-vector (N > 1) constructs H non-randomly using the elements
##        of X as parameters.
##        In both cases H is constructed via a product of N-1 Givens rotations.
##
##        Note: See Gragg (1986) for how to represent an N-by-N (complex)
##        unitary Hessenberg matrix with positive subdiagonal elements in terms
##        of 2N-1 real parameters (the Schur parametrization).
##        This M-file handles the real case only and is intended simply as a
##        convenient way to generate random or non-random orthogonal Hessenberg
##        matrices.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("rando", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("rando", @var{n}, @var{k})
## Create a random matrix with elements -1, 0 or 1.
##
##         A = RANDO(N, K) is a random N-by-N matrix with elements from
##         one of the following discrete distributions (default K = 1):
##           K = 1:  A(i,j) =  0 or 1    with equal probability,
##           K = 2:  A(i,j) = -1 or 1    with equal probability,
##           K = 3:  A(i,j) = -1, 0 or 1 with equal probability.
##         N may be a 2-vector, in which case the matrix is N(1)-by-N(2).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("randsvd", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("randsvd", @var{n}, @var{kappa})
## @deftypefnx {Function File} {@var{a} =} gallery ("randsvd", @var{n}, @var{kappa}, @var{mode})
## @deftypefnx {Function File} {@var{a} =} gallery ("randsvd", @var{n}, @var{kappa}, @var{mode}, @var{kl})
## @deftypefnx {Function File} {@var{a} =} gallery ("randsvd", @var{n}, @var{kappa}, @var{mode}, @var{kl}, @var{ku})
## Create a random matrix with pre-assigned singular values.
##
##       RANDSVD(N, KAPPA, MODE, KL, KU) is a (banded) random matrix of order N
##       with COND(A) = KAPPA and singular values from the distribution MODE.
##       N may be a 2-vector, in which case the matrix is N(1)-by-N(2).
##       Available types:
##              MODE = 1:   one large singular value,
##              MODE = 2:   one small singular value,
##              MODE = 3:   geometrically distributed singular values,
##              MODE = 4:   arithmetically distributed singular values,
##              MODE = 5:   random singular values with unif. dist. logarithm.
##       If omitted, MODE defaults to 3, and KAPPA defaults to SQRT(1/EPS).
##       If MODE < 0 then the effect is as for ABS(MODE) except that in the
##       original matrix of singular values the order of the diagonal entries
##       is reversed: small to large instead of large to small.
##       KL and KU are the lower and upper bandwidths respectively; if they
##       are omitted a full matrix is produced.
##       If only KL is present, KU defaults to KL.
##       Special case: if KAPPA < 0 then a random full symmetric positive
##                     definite matrix is produced with COND(A) = -KAPPA and
##                     eigenvalues distributed according to MODE.
##                     KL and KU, if present, are ignored.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("redheff", @var{n})
## Create a zero and ones matrix of Redheffer associated with the Riemann
## hypothesis.
##
##            A = REDHEFF(N) is an N-by-N matrix of 0s and 1s defined by
##                A(i,j) = 1 if j = 1 or if i divides j,
##                A(i,j) = 0 otherwise.
##            It has N - FLOOR(LOG2(N)) - 1 eigenvalues equal to 1,
##            a real eigenvalue (the spectral radius) approximately SQRT(N),
##            a negative eigenvalue approximately -SQRT(N),
##            and the remaining eigenvalues are provably ``small''.
##            Barrett and Jarvis (1992) conjecture that
##              ``the small eigenvalues all lie inside the unit circle
##                ABS(Z) = 1'',
##            and a proof of this conjecture, together with a proof that some
##            eigenvalue tends to zero as N tends to infinity, would yield
##            a new proof of the prime number theorem.
##            The Riemann hypothesis is true if and only if
##            DET(A) = O( N^(1/2+epsilon) ) for every epsilon > 0
##                                              (`!' denotes factorial).
##            See also RIEMANN.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("riemann", @var{n})
## Create a matrix associated with the Riemann hypothesis.
##
##            A = RIEMANN(N) is an N-by-N matrix for which the
##            Riemann hypothesis is true if and only if
##            DET(A) = O( N! N^(-1/2+epsilon) ) for every epsilon > 0
##                                              (`!' denotes factorial).
##            A = B(2:N+1, 2:N+1), where
##            B(i,j) = i-1 if i divides j and -1 otherwise.
##            Properties include, with M = N+1:
##               Each eigenvalue E(i) satisfies ABS(E(i)) <= M - 1/M.
##               i <= E(i) <= i+1 with at most M-SQRT(M) exceptions.
##               All integers in the interval (M/3, M/2] are eigenvalues.
##
##            See also REDHEFF.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("ris", @var{n})
## Create a symmetric Hankel matrix.
##
##           A = RIS(N) is the symmetric N-by-N Hankel matrix with
##                          A(i,j) = 0.5/(N-i-j+1.5).
##           The eigenvalues of A cluster around PI/2 and -PI/2.
##
##           Invented by F.N. Ris.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("smoke", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("smoke", @var{n}, @var{k})
## Create a complex matrix, with a `smoke ring' pseudospectrum.
##
##           SMOKE(N) is an N-by-N matrix with 1s on the
##           superdiagonal, 1 in the (N,1) position, and powers of
##           roots of unity along the diagonal.
##           SMOKE(N, 1) is the same except for a zero (N,1) element.
##           The eigenvalues of SMOKE(N, 1) are the N'th roots of unity;
##           those of SMOKE(N) are the N'th roots of unity times 2^(1/N).
##
##           Try PS(SMOKE(32)).  For SMOKE(N, 1) the pseudospectrum looks
##           like a sausage folded back on itself.
##           GERSH(SMOKE(N, 1)) is interesting.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{t} =} gallery ("toeppd", @var{n})
## @deftypefnx {Function File} {@var{t} =} gallery ("toeppd", @var{n}, @var{m})
## @deftypefnx {Function File} {@var{t} =} gallery ("toeppd", @var{n}, @var{m}, @var{w})
## @deftypefnx {Function File} {@var{t} =} gallery ("toeppd", @var{n}, @var{m}, @var{w}, @var{theta})
## Create a symmetric positive definite Toeplitz matrix.
##
##          TOEPPD(N, M, W, THETA) is an N-by-N symmetric positive (semi-)
##          definite (SPD) Toeplitz matrix, comprised of the sum of M rank 2
##          (or, for certain THETA, rank 1) SPD Toeplitz matrices.
##          Specifically,
##                  T = W(1)*T(THETA(1)) + ... + W(M)*T(THETA(M)),
##          where T(THETA(k)) has (i,j) element COS(2*PI*THETA(k)*(i-j)).
##          Defaults: M = N, W = RAND(M,1), THETA = RAND(M,1).
## @end deftypefn
##
## @deftypefn  {Function File} {@var{p} =} gallery ("toeppen", @var{n})
## @deftypefnx {Function File} {@var{p} =} gallery ("toeppen", @var{n}, @var{a})
## @deftypefnx {Function File} {@var{p} =} gallery ("toeppen", @var{n}, @var{a}, @var{b})
## @deftypefnx {Function File} {@var{p} =} gallery ("toeppen", @var{n}, @var{a}, @var{b}, @var{c})
## @deftypefnx {Function File} {@var{p} =} gallery ("toeppen", @var{n}, @var{a}, @var{b}, @var{c}, @var{d})
## @deftypefnx {Function File} {@var{p} =} gallery ("toeppen", @var{n}, @var{a}, @var{b}, @var{c}, @var{d}, @var{e})
## Create a pentadiagonal Toeplitz matrix (sparse).
##
##           P = TOEPPEN(N, A, B, C, D, E) is the N-by-N pentadiagonal
##           Toeplitz matrix with diagonals composed of the numbers
##           A =: P(3,1), B =: P(2,1), C =: P(1,1), D =: P(1,2), E =: P(1,3).
##           Default: (A,B,C,D,E) = (1,-10,0,10,1) (a matrix of Rutishauser).
##                     This matrix has eigenvalues lying approximately on
##                     the line segment 2*cos(2*t) + 20*i*sin(t).
##
##           Interesting plots are
##           PS(FULL(TOEPPEN(32,0,1,0,0,1/4)))  - `triangle'
##           PS(FULL(TOEPPEN(32,0,1/2,0,0,1)))  - `propeller'
##           PS(FULL(TOEPPEN(32,0,1/2,1,1,1)))  - `fish'
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("tridiag", @var{x}, @var{y}, @var{z})
## @deftypefnx {Function File} {@var{a} =} gallery ("tridiag", @var{n})
## @deftypefnx {Function File} {@var{a} =} gallery ("tridiag", @var{n}, @var{c}, @var{d}, @var{e})
## Create a tridiagonal matrix (sparse).
##
## @var{a} is the tridiagonal matrix with subdiagonal @var{x},
## diagonal @var{y}, and superdiagonal @var{z}.   @var{x} and @var{z}
## must be vectors of length one less than @var{y}.
##
## Alternatively if @var{c}, @var{d}, and @var{e} are all scalars, yields
## the Toeplitz tridiagonal matrix of order @var{n} with subdiagonal
## elements @var{c}, diagonal elements @var{d}, and superdiagonal
## elements @var{e}.  This matrix has eigenvalues (Todd 1977)
## @code{@var{d} + 2*sqrt(@var{c}*@var{e})*cos(k*pi/(@var{n}+1))}, for k=1:N.
## If unspecified, @var{c}, @var{d}, @var{e}, default to -1, 2, and -1
## respectively, which generates a symmetric positive definite
## M-matrix (the negative of the second difference matrix).
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
## @code{abs (@var{alpha})^@var{n}*sin (pi/(4 * @var{n}-2))}.
##
## Adding @code{-2^(2-@var{n})} to the (@var{n},1) element makes
## @var{c} singular, as does adding @code{-2^(1-@var{n})} to all elements
## in the first column.
## @end deftypefn
##
## @deftypefn  {Function File} {@var{a} =} gallery ("wathen", @var{nx}, @var{ny})
## @deftypefnx {Function File} {@var{a} =} gallery ("wathen", @var{nx}, @var{ny}, @var{k})
## Create the Wathen matrix.
##
##     The Wathen matrix is a finite element matrix which is sparse.
##
##     The entries of the matrix depend in part on a physical quantity
##     related to density.  That density is here assigned random values between
##     0 and 100.
##
##     A = WATHEN ( NX, NY ) is a sparse random N-by-N finite element matrix
##     where N = 3*NX*NY + 2*NX + 2*NY + 1.
##
##     A is the consistent mass matrix for a regular NX-by-NY
##     grid of 8-node (serendipity) elements in 2 space dimensions.
##
##     Here is an illustration for NX = 3, NX = 2:
##
##      23-24-25-26-27-28-29
##       |     |     |     |
##      19    20    21    22
##       |     |     |     |
##      12-13-14-15-16-17-18
##       |     |     |     |
##       8     9    10    11
##       |     |     |     |
##       1--2--3--4--5--6--7
##
##     For this example, the total number of nodes is, as expected,
##
##       N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29.
##
##     A is symmetric positive definite for any (positive) values of
##     the density, RHO(NX,NY), which is chosen randomly in this routine.
##
##     In particular, if D = DIAG(DIAG(A)), then
##       0.25 <= EIG(INV(D)*A) <= 4.5
##     for any positive integers NX and NY and any densities RHO(NX,NY).
##
##     A = WATHEN ( NX, NY, 1 ) returns the diagonally scaled matrix.
##
##   Parameters:
##
##     Input, integer NX, NY, the number of elements in the X and Y directions
##     of the finite element grid.  NX and NY must each be at least 1.
##
##     Optional input, integer K, is used to request that the diagonally scaled
##     version of the matrix be returned.  This happens if K is specified with
##     the value 1.
##
##     Output, sparse real A(N,N), the matrix.  The dimension N is determined by
##     NX and NY, as described above.  A is stored in the MATLAB sparse matrix
##     format.  
## @end deftypefn
##
## @deftypefn  {Function File} {[@var{a}, @var{b}] =} gallery ("wilk", @var{n})
## Create various specific matrices devised/discussed by Wilkinson.
##
##        [A, b] = WILK(N) is the matrix or system of order N.
##        N = 3: upper triangular system Ux=b illustrating inaccurate solution.
##        N = 4: lower triangular system Lx=b, ill-conditioned.
##        N = 5: HILB(6)(1:5,2:6)*1.8144.  Symmetric positive definite.
##        N = 21: W21+, tridiagonal.   Eigenvalue problem.
## @end deftypefn

## Code for most of the individual matrices (except binomial, gcdmat,
## integerdata, leslie, normaldata, randcolu, randcorr, randjorth, sampling,
## uniformdata) by Nicholas .J. Higham <Nicholas.J.Higham@manchester.ac.uk>
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
    case "cauchy"     , matrix = cauchy      (varargin{:});
    case "chebspec"   , matrix = chebspec    (varargin{:});
    case "chebvand"   , matrix = chebvand    (varargin{:});
    case "chow"       , matrix = chow        (varargin{:});
    case "circul"     , matrix = circul      (varargin{:});
    case "clement"    , matrix = clement     (varargin{:});
    case "compar"     , matrix = compar      (varargin{:});
    case "condex"     , matrix = condex      (varargin{:});
    case "cycol"      , matrix = cycol       (varargin{:});
    case "dorr"       , [matrix, varargout{1:nargout-1}] = dorr (varargin{:});
    case "dramadah"   , matrix = dramadah    (varargin{:});
    case "fiedler"    , matrix = fiedler     (varargin{:});
    case "forsythe"   , matrix = forsythe    (varargin{:});
    case "frank"      , matrix = frank       (varargin{:});
    case "gearmat"    , matrix = gearmat     (varargin{:});
    case "gcdmat"     , matrix = gcdmat      (varargin{:});
    case "grcar"      , matrix = grcar       (varargin{:});
    case "hanowa"     , matrix = hanowa      (varargin{:});
    case "house"      , matrix = house       (varargin{:});
    case "integerdata"
      error ("gallery: matrix %s not implemented.", name);
    case "invhess"    , matrix = invhess     (varargin{:});
    case "invol"      , matrix = invol       (varargin{:});
    case "ipjfact"    , [matrix, varargout{1}] = ipjfact (varargin{:});
    case "jordbloc"   , matrix = jordbloc    (varargin{:});
    case "kahan"      , matrix = kahan       (varargin{:});
    case "kms"        , matrix = kms         (varargin{:});
    case "krylov"     , matrix = krylov      (varargin{:});
    case "lauchli"    , matrix = lauchli     (varargin{:});
    case "lehmer"     , matrix = lehmer      (varargin{:});
    case "leslie"
      error ("gallery: matrix %s not implemented.", name);
    case "lesp"       , matrix = lesp        (varargin{:});
    case "lotkin"     , matrix = lotkin      (varargin{:});
    case "minij"      , matrix = minij       (varargin{:});
    case "moler"      , matrix = moler       (varargin{:});
    case "neumann"    , [matrix, varargout{1}] = neumann (varargin{:});
    case "normaldata"
      error ("gallery: matrix %s not implemented.", name);
    case "orthog"     , matrix = orthog      (varargin{:});
    case "parter"     , matrix = parter      (varargin{:});
    case "pei"        , matrix = pei         (varargin{:});
    case "poisson"    , matrix = poisson     (varargin{:});
    case "prolate"    , matrix = prolate     (varargin{:});
    case "randcolu"
      error ("gallery: matrix %s not implemented.", name);
    case "randcorr"
      error ("gallery: matrix %s not implemented.", name);
    case "randhess"    , matrix = randhess    (varargin{:});
    case "randjorth"
      error ("gallery: matrix %s not implemented.", name);
    case "rando"       , matrix = rando       (varargin{:});
    case "randsvd"     , matrix = randsvd     (varargin{:});
    case "redheff"     , matrix = redheff     (varargin{:});
    case "riemann"     , matrix = riemann     (varargin{:});
    case "ris"         , matrix = ris         (varargin{:});
    case "sampling"
      error ("gallery: matrix %s not implemented.", name);
    case "smoke"       , matrix = smoke       (varargin{:});
    case "toeppd"      , matrix = toeppd      (varargin{:});
    case "toeppen"     , matrix = toeppen     (varargin{:});
    case "tridiag"     , matrix = tridiag     (varargin{:});
    case "triw"        , matrix = triw        (varargin{:});
    case "uniformdata"
      error ("gallery: matrix %s not implemented.", name);
    case "wathen"      , matrix = wathen      (varargin{:});
    case "wilk"        , [matrix, varargout{1}] = wilk (varargin{:});
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
    p = linspace (0, 1, n);
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
    endfor

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
    error ("gallery: K must be a numeric scalar for condex matrix.");
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
    error ("gallery: K must be a scalar for cycol matrix.");
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

function [c, d, e] = dorr (n, theta = 0.01)
  ## DORR  Dorr matrix - diagonally dominant, ill conditioned, tridiagonal.
  ##       [C, D, E] = DORR(N, THETA) returns the vectors defining a row diagonally
  ##       dominant, tridiagonal M-matrix that is ill conditioned for small
  ##       values of the parameter THETA >= 0.
  ##       If only one output parameter is supplied then
  ##       C = FULL(TRIDIAG(C,D,E)), i.e., the matrix iself is returned.
  ##       The columns of INV(C) vary greatly in norm.  THETA defaults to 0.01.
  ##       The amount of diagonal dominance is given by (ignoring rounding errors):
  ##             COMP(C)*ONES(N,1) = THETA*(N+1)^2 * [1 0 0 ... 0 1]'.
  ##
  ##       Reference:
  ##       F.W. Dorr, An example of ill-conditioning in the numerical
  ##       solution of singular perturbation problems, Math. Comp., 25 (1971),
  ##       pp. 271-283.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 or 2 arguments are required for dorr matrix.");
  elseif (! isscalar (n) || ! isnumeric (n))
    error ("gallery: N must be a numeric scalar for dorr matrix.");
  elseif (! isscalar (theta) || ! isnumeric (theta))
    error ("gallery: THETA must be a numeric scalar for dorr matrix.");
  endif

  c = zeros (n,1);
  e = c;
  d = c;
  ##  All length n for convenience.  Make c, e of length n-1 later.

  h = 1/(n+1);
  m = floor ((n+1)/2);
  term = theta/h^2;

  i = (1:m)';
  c(i) = -term * ones (m,1);
  e(i) = c(i) - (0.5-i*h)/h;
  d(i) = -(c(i) + e(i));

  i = (m+1:n)';
  e(i) = -term * ones (n-m,1);
  c(i) = e(i) + (0.5-i*h)/h;
  d(i) = -(c(i) + e(i));

  c = c(2:n);
  e = e(1:n-1);

  if (nargout <= 1)
    c = tridiag (c, d, e);
  endif
endfunction

function A = dramadah (n, k = 1)
  ## DRAMADAH  A (0,1) matrix whose inverse has large integer entries.
  ##           An anti-Hadamard matrix A is a matrix with elements 0 or 1 for
  ##           which MU(A) := NORM(INV(A),'FRO') is maximal.
  ##           A = DRAMADAH(N, K) is an N-by-N (0,1) matrix for which MU(A) is
  ##           relatively large, although not necessarily maximal.
  ##           Available types (the default is K = 1):
  ##           K = 1: A is Toeplitz, with ABS(DET(A)) = 1, and MU(A) > c(1.75)^N,
  ##                  where c is a constant.
  ##           K = 2: A is upper triangular and Toeplitz.
  ##           The inverses of both types have integer entries.
  ##
  ##           Another interesting (0,1) matrix:
  ##           K = 3: A has maximal determinant among (0,1) lower Hessenberg
  ##           matrices: det(A) = the n'th Fibonacci number.  A is Toeplitz.
  ##           The eigenvalues have an interesting distribution in the complex
  ##           plane.
  ##
  ##           References:
  ##           R.L. Graham and N.J.A. Sloane, Anti-Hadamard matrices,
  ##              Linear Algebra and Appl., 62 (1984), pp. 113-137.
  ##           L. Ching, The maximum determinant of an nxn lower Hessenberg
  ##              (0,1) matrix, Linear Algebra and Appl., 183 (1993), pp. 147-153.

  if (k == 1)  # Toeplitz
    c = ones(n,1);
    for i=2:4:n
      m = min(1,n-i);
      c(i:i+m) = zeros(m+1,1);
    endfor
    r = zeros(n,1);
    r(1:4) = [1 1 0 1];
    if (n < 4)
      r = r(1:n);
    endif
    A = toeplitz(c,r);

  elseif (k == 2)  # Upper triangular and Toeplitz
    c = zeros(n,1);
    c(1) = 1;
    r = ones(n,1);
    for i=3:2:n
      r(i) = 0;
    endfor
    A = toeplitz(c,r);

  elseif k == 3  # Lower Hessenberg.
    c = ones(n,1);
    for i=2:2:n
      c(i)=0;
    endfor
    A = toeplitz(c, [1 1 zeros(1,n-2)]);
  else
    error ("gallery: invalid K for dramadah matrix.")
  endif
endfunction

function A = fiedler (c)
  ## FIEDLER  Fiedler matrix - symmetric.
  ##          A = FIEDLER(C), where C is an n-vector, is the n-by-n symmetric
  ##          matrix with elements ABS(C(i)-C(j)).
  ##          Special case: if C is a scalar, then A = FIEDLER(1:C)
  ##                        (i.e. A(i,j) = ABS(i-j)).
  ##          Properties:
  ##            FIEDLER(N) has a dominant positive eigenvalue and all the other
  ##                       eigenvalues are negative (Szego, 1936).
  ##            Explicit formulas for INV(A) and DET(A) are given by Todd (1977)
  ##            and attributed to Fiedler.  These indicate that INV(A) is
  ##            tridiagonal except for nonzero (1,n) and (n,1) elements.
  ##            [I think these formulas are valid only if the elements of
  ##            C are in increasing or decreasing order---NJH.]
  ##
  ##            References:
  ##            G. Szego, Solution to problem 3705, Amer. Math. Monthly,
  ##               43 (1936), pp. 246-259.
  ##            J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
  ##               Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.

  n = max(size(c));

  ##   Handle scalar c.
  if n == 1
    n = c;
    c = 1:n;
  endif

  c = c(:).';                    % Ensure c is a row vector.

  A = ones(n,1)*c;
  A = abs(A - A.');              % NB. array transpose.
endfunction

function A = forsythe (n, alpha = sqrt (eps), lambda = 0)
  ## FORSYTHE  Forsythe matrix - a perturbed Jordan block.
  ##           FORSYTHE(N, ALPHA, LAMBDA) is the N-by-N matrix equal to
  ##           JORDBLOC(N, LAMBDA) except it has an ALPHA in the (N,1) position.
  ##           It has the characteristic polynomial
  ##                   DET(A-t*EYE) = (LAMBDA-t)^N - (-1)^N ALPHA.
  ##           ALPHA defaults to SQRT(EPS) and LAMBDA to 0.

  A = jordbloc (n, lambda);
  A(n,1) = alpha;
endfunction

function F = frank (n, k = 0)
  ## FRANK   Frank matrix---ill conditioned eigenvalues.
  ##         F = FRANK(N, K) is the Frank matrix of order N.  It is upper
  ##         Hessenberg with determinant 1.  K = 0 is the default; if K = 1 the
  ##         elements are reflected about the anti-diagonal (1,N)--(N,1).
  ##         F has all positive eigenvalues and they occur in reciprocal pairs
  ##         (so that 1 is an eigenvalue if N is odd).
  ##         The eigenvalues of F may be obtained in terms of the zeros of the
  ##         Hermite polynomials.
  ##         The FLOOR(N/2) smallest eigenvalues of F are ill conditioned,
  ##         the more so for bigger N.
  ##
  ##         DET(FRANK(N)') comes out far from 1 for large N---see Frank (1958)
  ##         and Wilkinson (1960) for discussions.
  ##
  ##         This version incorporates improvements suggested by W. Kahan.
  ##
  ##         References:
  ##         W.L. Frank, Computing eigenvalues of complex matrices by determinant
  ##            evaluation and by methods of Danilewski and Wielandt, J. Soc.
  ##            Indust. Appl. Math., 6 (1958), pp. 378-392 (see pp. 385, 388).
  ##         G.H. Golub and J.H. Wilkinson, Ill-conditioned eigensystems and the
  ##            computation of the Jordan canonical form, SIAM Review, 18 (1976),
  ##              pp. 578-619 (Section 13).
  ##         H. Rutishauser, On test matrices, Programmation en Mathematiques
  ##            Numeriques, Editions Centre Nat. Recherche Sci., Paris, 165,
  ##            1966, pp. 349-365.  Section 9.
  ##         J.H. Wilkinson, Error analysis of floating-point computation,
  ##            Numer. Math., 2 (1960), pp. 319-340 (Section 8).
  ##         J.H. Wilkinson, The Algebraic Eigenvalue Problem, Oxford University
  ##            Press, 1965 (pp. 92-93).
  ##         The next two references give details of the eigensystem, as does
  ##         Rutishauser (see above).
  ##         P.J. Eberlein, A note on the matrices denoted by B_n, SIAM J. Appl.
  ##            Math., 20 (1971), pp. 87-92.
  ##         J.M. Varah, A generalization of the Frank matrix, SIAM J. Sci. Stat.
  ##            Comput., 7 (1986), pp. 835-839.

  p = n:-1:1;
  F = triu( p( ones(n,1), :) - diag( ones(n-1,1), -1), -1 );
  if (k != 0)
    F = F(p,p)';
  endif
endfunction

function c = gcdmat (n)
  if (nargin != 1)
    error ("gallery: unspecified N for gcdmat matrix.");
  elseif (! isscalar (n) || ! isnumeric (n))
    error ("gallery: N must be a numeric scalar for gcdmat matrix.");
  endif
  c = gcd (repmat ((1:n)', [1 n]), repmat (1:n, [n 1]));
endfunction

function A = gearmat (n, i = n, j = -n)
  ## NOTE: this function was named gearm in the original Test Matrix Toolbox
  ## GEARMAT   Gear matrix.
  ##         A = GEARMAT(N,I,J) is the N-by-N matrix with ones on the sub- and
  ##         super-diagonals, SIGN(I) in the (1,ABS(I)) position, SIGN(J)
  ##         in the (N,N+1-ABS(J)) position, and zeros everywhere else.
  ##         Defaults: I = N, j = -N.
  ##         All eigenvalues are of the form 2*COS(a) and the eigenvectors
  ##         are of the form [SIN(w+a), SIN(w+2a), ..., SIN(w+Na)].
  ##         The values of a and w are given in the reference below.
  ##         A can have double and triple eigenvalues and can be defective.
  ##         GEARMAT(N) is singular.
  ##
  ##         (GEAR is a Simulink function, hence GEARMAT for Gear matrix.)
  ##         Reference:
  ##         C.W. Gear, A simple set of test matrices for eigenvalue programs,
  ##         Math. Comp., 23 (1969), pp. 119-125.

  if ~(i~=0 & abs(i)<=n & j~=0 & abs(j)<=n)
    error('Invalid I and J parameters')
  endif

  A = diag(ones(n-1,1),-1) + diag(ones(n-1,1),1);
  A(1, abs(i)) = sign(i);
  A(n, n+1-abs(j)) = sign(j);
endfunction

function G = grcar (n, k = 3)
  ## GRCAR     Grcar matrix - a Toeplitz matrix with sensitive eigenvalues.
  ##           GRCAR(N, K) is an N-by-N matrix with -1s on the
  ##           subdiagonal, 1s on the diagonal, and K superdiagonals of 1s.
  ##           The default is K = 3.  The eigenvalues of this matrix form an
  ##           interesting pattern in the complex plane (try PS(GRCAR(32))).
  ##
  ##           References:
  ##           J.F. Grcar, Operator coefficient methods for linear equations,
  ##                Report SAND89-8691, Sandia National Laboratories, Albuquerque,
  ##                New Mexico, 1989 (Appendix 2).
  ##           N.M. Nachtigal, L. Reichel and L.N. Trefethen, A hybrid GMRES
  ##                algorithm for nonsymmetric linear systems, SIAM J. Matrix Anal.
  ##                Appl., 13 (1992), pp. 796-825.

  G = tril(triu(ones(n)), k) - diag(ones(n-1,1), -1);
endfunction

function A = hanowa (n, d = -1)
  ## HANOWA  A matrix whose eigenvalues lie on a vertical line in the complex plane.
  ##         HANOWA(N, d) is the N-by-N block 2x2 matrix (thus N = 2M must be even)
  ##                       [d*EYE(M)   -DIAG(1:M)
  ##                        DIAG(1:M)   d*EYE(M)]
  ##         It has complex eigenvalues lambda(k) = d +/- k*i  (1 <= k <= M).
  ##         Parameter d defaults to -1.
  ##
  ##         Reference:
  ##         E. Hairer, S.P. Norsett and G. Wanner, Solving Ordinary
  ##         Differential Equations I: Nonstiff Problems, Springer-Verlag,
  ##         Berlin, 1987. (pp. 86-87)

  m = n/2;
  if round(m) ~= m
    error('N must be even.')
  end

  A = [ d*eye(m) -diag(1:m)
        diag(1:m)  d*eye(m)];
endfunction

function [v, beta] = house (x)
  ## HOUSE   Householder matrix.
  ##         If [v, beta] = HOUSE(x) then H = EYE - beta*v*v' is a Householder
  ##         matrix such that Hx = -sign(x(1))*norm(x)*e_1.
  ##         NB: If x = 0 then v = 0, beta = 1 is returned.
  ##             x can be real or complex.
  ##             sign(x) := exp(i*arg(x)) ( = x./abs(x) when x ~= 0).
  ##
  ##         Theory: (textbook references Golub & Van Loan 1989, 38-43;
  ##                  Stewart 1973, 231-234, 262; Wilkinson 1965, 48-50).
  ##         Hx = y: (I - beta*v*v')x = -s*e_1.
  ##         Must have |s| = norm(x), v = x+s*e_1, and
  ##         x'y = x'Hx =(x'Hx)' real => arg(s) = arg(x(1)).
  ##         So take s = sign(x(1))*norm(x) (which avoids cancellation).
  ##         v'v = (x(1)+s)^2 + x(2)^2 + ... + x(n)^2
  ##             = 2*norm(x)*(norm(x) + |x(1)|).
  ## 
  ##         References:
  ##         G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
  ##            Johns Hopkins University Press, Baltimore, Maryland, 1989.
  ##         G.W. Stewart, Introduction to Matrix Computations, Academic Press,
  ##            New York, 1973,
  ##         J.H. Wilkinson, The Algebraic Eigenvalue Problem, Oxford University
  ##            Press, 1965.

  [m, n] = size(x);
  if (n > 1)
    error('Argument must be a column vector.')
  endif

  s = norm(x) * (sign(x(1)) + (x(1)==0));    % Modification for sign(0)=1.
  v = x;
  if (s == 0)
    ## Quit if x is the zero vector.
    beta = 1;
  else
    v(1) = v(1) + s;
    beta = 1/(s'*v(1));                        % NB the conjugated s.
    ##  beta = 1/(abs(s)*(abs(s)+abs(x(1)) would guarantee beta real.
    ##  But beta as above can be non-real (due to rounding) only when x is complex.
  endif
endfunction

function A = invhess (x, y)
  ## INVHESS  Inverse of an upper Hessenberg matrix.
  ##          INVHESS(X, Y), where X is an N-vector and Y an N-1 vector,
  ##          is the matrix whose lower triangle agrees with that of
  ##          ONES(N,1)*X' and whose strict upper triangle agrees with
  ##          that of [1 Y]*ONES(1,N).
  ##          The matrix is nonsingular if X(1) ~= 0 and X(i+1) ~= Y(i)
  ##          for all i, and its inverse is an upper Hessenberg matrix.
  ##          If Y is omitted it defaults to -X(1:N-1).
  ##          Special case: if X is a scalar INVHESS(X) is the same as
  ##          INVHESS(1:X).
  ##
  ##          References:
  ##          F.N. Valvi and V.S. Geroyannis, Analytic inverses and
  ##              determinants for a class of matrices, IMA Journal of Numerical
  ##              Analysis, 7 (1987), pp. 123-128.
  ##          W.-L. Cao and W.J. Stewart, A note on inverses of Hessenberg-like
  ##              matrices, Linear Algebra and Appl., 76 (1986), pp. 233-240.
  ##          Y. Ikebe, On inverses of Hessenberg matrices, Linear Algebra and
  ##              Appl., 24 (1979), pp. 93-97.
  ##          P. Rozsa, On the inverse of band matrices, Integral Equations and
  ##              Operator Theory, 10 (1987), pp. 82-95.

  n = max(size(x));
  ##   Handle scalar x.
  if n == 1
    n = x;
    x = 1:n;
  endif
  x = x(:);

  if nargin < 2, y = -x; end
  y = y(:);

  ##  On next line, z = x'; A = z(ones(n,1),:) would be more efficient.
  A = ones(n,1)*x';  
  for j=2:n
    A(1:j-1,j) = y(1:j-1);
  endfor
endfunction

function A = invol (n)
  ## INVOL   An involutory matrix.
  ##         A = INVOL(N) is an N-by-N involutory (A*A = EYE(N)) and
  ##         ill-conditioned matrix.
  ##         It is a diagonally scaled version of HILB(N).
  ##         NB: B = (EYE(N)-A)/2 and B = (EYE(N)+A)/2 are idempotent (B*B = B).
  ##
  ##         Reference:
  ##         A.S. Householder and J.A. Carpenter, The singular values
  ##         of involutory and of idempotent matrices, Numer. Math. 5 (1963),
  ##         pp. 234-237.

  A = hilb(n);

  d = -n;
  A(:, 1) = d*A(:, 1);

  for i = 1:n-1
    d = -(n+i)*(n-i)*d/(i*i);
    A(i+1, :) = d*A(i+1, :);
  endfor
endfunction

function [A, detA] = ipjfact (n, k = 0)
  ## IPJFACT   A Hankel matrix with factorial elements.
  ##           A = IPJFACT(N, K) is the matrix with
  ##                     A(i,j) = (i+j)!    (K = 0, default)
  ##                     A(i,j) = 1/(i+j)!  (K = 1)
  ##           Both are Hankel matrices.
  ##           The determinant and inverse are known explicitly.
  ##           If a second output argument is present, d = DET(A) is returned:
  ##           [A, d] = IPJFACT(N, K);
  ##
  ##           Suggested by P. R. Graves-Morris.
  ## 
  ##           Reference:
  ##           M.J.C. Gover, The explicit inverse of factorial Hankel matrices,
  ##           Dept. of Mathematics, University of Bradford, 1993.

  c = cumprod(2:n+1);
  d = cumprod(n+1:2*n) * c(n-1);

  A = hankel(c, d);

  if k == 1
    A = ones(n)./A;
  endif

  if nargout == 2
    d = 1;
    if k == 0
      for i=1:n-1
        d = d*prod(1:i+1)*prod(1:n-i);
      endfor
      d = d*prod(1:n+1);
    else
      for i=0:n-1
        d = d*prod(1:i)/prod(1:n+1+i);
      endfor
      if rem(n*(n-1)/2,2), d = -d; endif
    endif
    detA = d;
  endif
endfunction

function J = jordbloc (n, lambda = 1)
  ## JORDBLOC  Jordan block.
  ##           JORDBLOC(N, LAMBDA) is the N-by-N Jordan block with eigenvalue
  ##           LAMBDA.  LAMBDA = 1 is the default.

  J = lambda*eye(n) + diag(ones(n-1,1),1);
endfunction

function U = kahan (n, theta = 1.2, pert = 25)
  ## KAHAN  Kahan matrix - upper trapezoidal.
  ##        KAHAN(N, THETA) is an upper trapezoidal matrix
  ##        that has some interesting properties regarding estimation of
  ##        condition and rank.
  ##        The matrix is N-by-N unless N is a 2-vector, in which case it
  ##        is N(1)-by-N(2).
  ##        The parameter THETA defaults to 1.2.
  ##        The useful range of THETA is 0 < THETA < PI.
  ##
  ##        To ensure that the QR factorization with column pivoting does not
  ##        interchange columns in the presence of rounding errors, the diagonal
  ##        is perturbed by PERT*EPS*diag( [N:-1:1] ).
  ##        The default is PERT = 25, which ensures no interchanges for KAHAN(N)
  ##        up to at least N = 90 in IEEE arithmetic.
  ##        KAHAN(N, THETA, PERT) uses the given value of PERT.
  ##
  ##        The inverse of KAHAN(N, THETA) is known explicitly: see
  ##        Higham (1987, p. 588), for example.
  ##        The diagonal perturbation was suggested by Christian Bischof.
  ##
  ##        References:
  ##        W. Kahan, Numerical linear algebra, Canadian Math. Bulletin,
  ##           9 (1966), pp. 757-801.
  ##        N.J. Higham, A survey of condition number estimation for
  ##           triangular matrices, SIAM Review, 29 (1987), pp. 575-596.

  r = n(1);              % Parameter n specifies dimension: r-by-n.
  n = n(max(size(n)));

  s = sin(theta);
  c = cos(theta);

  U = eye(n) - c*triu(ones(n), 1);
  U = diag(s.^[0:n-1])*U + pert*eps*diag( [n:-1:1] );
  if r > n
    U(r,n) = 0;         % Extend to an r-by-n matrix.
  else
    U = U(1:r,:);       % Reduce to an r-by-n matrix.
  endif
endfunction

function A = kms (n, rho = 0.5)
  ## KMS   Kac-Murdock-Szego Toeplitz matrix.
  ##       A = KMS(N, RHO) is the N-by-N Kac-Murdock-Szego Toeplitz matrix with
  ##       A(i,j) = RHO^(ABS((i-j))) (for real RHO).
  ##       If RHO is complex, then the same formula holds except that elements
  ##       below the diagonal are conjugated.
  ##       RHO defaults to 0.5.
  ##       Properties:
  ##          A has an LDL' factorization with
  ##                   L = INV(TRIW(N,-RHO,1)'),
  ##                   D(i,i) = (1-ABS(RHO)^2)*EYE(N) except D(1,1) = 1.
  ##          A is positive definite if and only if 0 < ABS(RHO) < 1.
  ##          INV(A) is tridiagonal.
  ##
  ##        Reference:
  ##        W.F. Trench, Numerical solution of the eigenvalue problem
  ##        for Hermitian Toeplitz matrices, SIAM J. Matrix Analysis and Appl.,
  ##        10 (1989), pp. 135-146 (and see the references therein).

  A = (1:n)'*ones(1,n);
  A = abs(A - A');
  A = rho .^ A;
  if imag(rho)
    A = conj(tril(A,-1)) + triu(A);
  endif
endfunction

function B = krylov (A, x, j)
  ## KRYLOV    Krylov matrix.
  ##           KRYLOV(A, x, j) is the Krylov matrix
  ##                [x, Ax, A^2x, ..., A^(j-1)x],
  ##           where A is an n-by-n matrix and x is an n-vector.
  ##           Defaults: x = ONES(n,1), j = n.
  ##           KRYLOV(n) is the same as KRYLOV(RANDN(n)).
  ##
  ##           Reference:
  ##           G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
  ##           Johns Hopkins University Press, Baltimore, Maryland, 1989, p. 369.

  [n, n] = size(A);

  if n == 1   % Handle special case A = scalar.
    n = A;
    A = randn(n);
  endif

  if nargin < 3, j = n; end
  if nargin < 2, x = ones(n,1); end


  B = ones(n,j);
  B(:,1) = x(:);
  for i=2:j
    B(:,i) = A*B(:,i-1);
  endfor
endfunction

function A = lauchli (n, mu = sqrt (eps))
  ## LAUCHLI   Lauchli matrix - rectangular.
  ##           LAUCHLI(N, MU) is the (N+1)-by-N matrix [ONES(1,N); MU*EYE(N))].
  ##           It is a well-known example in least squares and other problems
  ##           that indicates the dangers of forming A'*A.
  ##           MU defaults to SQRT(EPS).
  ##
  ##           Reference:
  ##           P. Lauchli, Jordan-Elimination und Ausgleichung nach
  ##           kleinsten Quadraten, Numer. Math, 3 (1961), pp. 226-240.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for lauchli matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error("gallery: N must be an integer for lauchli matrix.")
  elseif (! isscalar (mu))
    error("gallery: MU must be a scalar for lauchli matrix.")
  endif

  A = [ones(1, n)
       mu*eye(n) ];
endfunction

function A = lehmer (n)
  ## LEHMER  Lehmer matrix - symmetric positive definite.
  ##         A = LEHMER(N) is the symmetric positive definite N-by-N matrix with
  ##                          A(i,j) = i/j for j >= i.
  ##         A is totally nonnegative.  INV(A) is tridiagonal, and explicit
  ##         formulas are known for its entries.
  ##         N <= COND(A) <= 4*N*N.
  ##
  ##         References:
  ##         M. Newman and J. Todd, The evaluation of matrix inversion
  ##            programs, J. Soc. Indust. Appl. Math., 6 (1958), pp. 466-476.
  ##         Solutions to problem E710 (proposed by D.H. Lehmer): The inverse
  ##            of a matrix, Amer. Math. Monthly, 53 (1946), pp. 534-535.
  ##         J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
  ##            Birkhauser, Basel, and Academic Press, New York, 1977, p. 154.

  if (nargin != 1)
    error ("gallery: 1 argument is required for lehmer matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for lehmer matrix.");
  endif

  A = ones (n, 1) * (1:n);
  A = A./A';
  A = tril (A) + tril (A, -1)';
endfunction

function T = lesp (n)
  ## LESP   A tridiagonal matrix with real, sensitive eigenvalues.
  ##        LESP(N) is an N-by-N matrix whose eigenvalues are real and smoothly
  ##        distributed in the interval approximately [-2*N-3.5, -4.5].
  ##        The sensitivities of the eigenvalues increase exponentially as
  ##        the eigenvalues grow more negative.
  ##        The matrix is similar to the symmetric tridiagonal matrix with
  ##        the same diagonal entries and with off-diagonal entries 1,
  ##        via a similarity transformation with D = diag(1!,2!,...,N!).
  ##
  ##        References:
  ##        H.W.J. Lenferink and M.N. Spijker, On the use of stability regions in
  ##             the numerical analysis of initial value problems,
  ##             Math. Comp., 57 (1991), pp. 221-237.
  ##        L.N. Trefethen, Pseudospectra of matrices, in Numerical Analysis 1991,
  ##             Proceedings of the 14th Dundee Conference,
  ##             D.F. Griffiths and G.A. Watson, eds, Pitman Research Notes in
  ##             Mathematics, volume 260, Longman Scientific and Technical, Essex,
  ##             UK, 1992, pp. 234-266.

  if (nargin != 1)
    error ("gallery: 1 argument is required for lesp matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for lesp matrix.");
  endif

  x = 2:n;
  T = full (tridiag (ones (size (x)) ./x, -(2*[x n+1]+1), x));
endfunction

function A = lotkin (n)
  ## LOTKIN  Lotkin matrix.
  ##         A = LOTKIN(N) is the Hilbert matrix with its first row altered to
  ##         all ones.  A is unsymmetric, ill-conditioned, and has many negative
  ##         eigenvalues of small magnitude.
  ##         The inverse has integer entries and is known explicitly.
  ##
  ##         Reference:
  ##         M. Lotkin, A set of test matrices, MTAC, 9 (1955), pp. 153-161.

  if (nargin != 1)
    error ("gallery: 1 argument is required for lotkin matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for lotkin matrix.");
  endif

  A = hilb (n);
  A(1,:) = ones (1, n);
endfunction

function A = minij (n)
  ## MINIJ   Symmetric positive definite matrix MIN(i,j).
  ##         A = MINIJ(N) is the N-by-N symmetric positive definite matrix with
  ##         A(i,j) = MIN(i,j).
  ##         Properties, variations:
  ##         INV(A) is tridiagonal: it is minus the second difference matrix
  ##                     except its (N,N) element is 1.
  ##         2*A-ONES(N) (Givens' matrix) has tridiagonal inverse and
  ##                     eigenvalues .5*sec^2([2r-1)PI/4N], r=1:N.
  ##         (N+1)*ONES(N)-A also has a tridiagonal inverse.
  ##
  ##         References:
  ##         J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
  ##            Birkhauser, Basel, and Academic Press, New York, 1977, p. 158.
  ##         D.E. Rutherford, Some continuant determinants arising in physics and
  ##            chemistry---II, Proc. Royal Soc. Edin., 63, A (1952), pp. 232-241.
  ##            (For the eigenvalues of Givens' matrix.)

  if (nargin != 1)
    error ("gallery: 1 argument is required for minij matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for minij matrix.");
  endif

  A = min (ones (n, 1) * (1:n), (1:n)' * ones (1, n));
endfunction

function A = moler (n, alpha = -1)
  ## MOLER   Moler matrix - symmetric positive definite.
  ##         A = MOLER(N, ALPHA) is the symmetric positive definite N-by-N matrix
  ##         U'*U where U = TRIW(N, ALPHA).
  ##         For ALPHA = -1 (the default) A(i,j) = MIN(i,j)-2, A(i,i) = i.
  ##         A has one small eigenvalue.
  ##
  ##         Nash (1990) attributes the ALPHA = -1 matrix to Moler.
  ##
  ##         Reference:
  ##         J.C. Nash, Compact Numerical Methods for Computers: Linear
  ##         Algebra and Function Minimisation, second edition, Adam Hilger,
  ##         Bristol, 1990 (Appendix 1).

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for moler matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error("gallery: N must be an integer for moler matrix.")
  elseif (! isscalar (alpha))
    error("gallery: ALPHA must be a scalar for moler matrix.")
  endif

  A = triw (n, alpha)' * triw (n, alpha);
endfunction

function [A, T] = neumann (n)
  ## NEUMANN  Singular matrix from the discrete Neumann problem (sparse).
  ##          NEUMANN(N) is the singular, row diagonally dominant matrix resulting
  ##          from discretizing the Neumann problem with the usual five point
  ##          operator on a regular mesh.
  ##          It has a one-dimensional null space with null vector ONES(N,1).
  ##          The dimension N should be a perfect square, or else a 2-vector,
  ##          in which case the dimension of the matrix is N(1)*N(2).
  ##
  ##          Reference:
  ##          R.J. Plemmons, Regular splittings and the discrete Neumann
  ##          problem, Numer. Math., 25 (1976), pp. 153-161.


  if (nargin != 1)
    error ("gallery: 1 argument is required for neumann matrix.");
  elseif (! isnumeric (n) || all (numel (n) != [1 2]) || fix (n) != n)
    error ("gallery: N must be a 1 or 2 element integer for neumann matrix.");
  endif

  if (isscalar (n))
    m = sqrt (n);
    if (m^2 != n)
      error ("gallery: N must be a perfect square for neumann matrix.");
    endif
    n(1) = m;
    n(2) = m;
  endif

  T = tridiag (n(1), -1, 2, -1);
  T(1,2) = -2;
  T(n(1),n(1)-1) = -2;

  A = kron (T, eye (n(2))) + kron (eye (n(2)), T);
endfunction

function Q = orthog (n, k = 1)
  ## ORTHOG Orthogonal and nearly orthogonal matrices.
  ##        Q = ORTHOG(N, K) selects the K'th type of matrix of order N.
  ##        K > 0 for exactly orthogonal matrices, K < 0 for diagonal scalings of
  ##        orthogonal matrices.
  ##        Available types: (K = 1 is the default)
  ##        K = 1:  Q(i,j) = SQRT(2/(n+1)) * SIN( i*j*PI/(n+1) )
  ##                Symmetric eigenvector matrix for second difference matrix.
  ##        K = 2:  Q(i,j) = 2/SQRT(2*n+1)) * SIN( 2*i*j*PI/(2*n+1) )
  ##                Symmetric.
  ##        K = 3:  Q(r,s) = EXP(2*PI*i*(r-1)*(s-1)/n) / SQRT(n)  (i=SQRT(-1))
  ##                Unitary, the Fourier matrix.  Q^4 is the identity.
  ##                This is essentially the same matrix as FFT(EYE(N))/SQRT(N)!
  ##        K = 4:  Helmert matrix: a permutation of a lower Hessenberg matrix,
  ##                whose first row is ONES(1:N)/SQRT(N).
  ##        K = 5:  Q(i,j) = SIN( 2*PI*(i-1)*(j-1)/n ) + COS( 2*PI*(i-1)*(j-1)/n ).
  ##                Symmetric matrix arising in the Hartley transform.
  ##        K = -1: Q(i,j) = COS( (i-1)*(j-1)*PI/(n-1) )
  ##                Chebyshev Vandermonde-like matrix, based on extrema of T(n-1).
  ##        K = -2: Q(i,j) = COS( (i-1)*(j-1/2)*PI/n) )
  ##                Chebyshev Vandermonde-like matrix, based on zeros of T(n).
  ##
  ##        References:
  ##        N.J. Higham and D.J. Higham, Large growth factors in Gaussian
  ##             elimination with pivoting, SIAM J. Matrix Analysis and  Appl.,
  ##             10 (1989), pp. 155-164.
  ##        P. Morton, On the eigenvectors of Schur's matrix, J. Number Theory,
  ##             12 (1980), pp. 122-127. (Re. ORTHOG(N, 3))
  ##        H.O. Lancaster, The Helmert Matrices, Amer. Math. Monthly, 72 (1965),
  ##             pp. 4-12.
  ##        D. Bini and P. Favati, On a matrix algebra related to the discrete
  ##             Hartley transform, SIAM J. Matrix Anal. Appl., 14 (1993),
  ##             pp. 500-507.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for orthog matrix.");
  elseif (! isnumeric (n) || all (numel (n) != [1 2]) || fix (n) != n)
    error ("gallery: N must be a 1 or 2 element integer for orthog matrix.");
  elseif (! isnumeric (k) || ! isscalar (k))
    error ("gallery: K must be a numeric scalar for orthog matrix.");
  endif

  switch (k)
    case (1)
      ## E'vectors second difference matrix
      m = (1:n)'*(1:n) * (pi/(n+1));
      Q = sin (m) * sqrt (2/(n+1));

    case (2)
      m = (1:n)'*(1:n) * (2*pi/(2*n+1));
      Q = sin (m) * (2/ sqrt (2*n+1));

    case (3)
      ## Vandermonde based on roots of unity
      m = 0:n-1;
      Q = exp (m'*m*2*pi* sqrt (-1) / n) / sqrt (n);

    case (4)
      ## Helmert matrix
      Q = tril (ones (n));
      Q(1,2:n) = ones (1, n-1);
      for i = 2:n
        Q(i,i) = -(i-1);
      end
      Q = diag (sqrt ([n 1:n-1] .* [1:n])) \ Q;

    case (5)
      ## Hartley matrix
      m = (0:n-1)'*(0:n-1) * (2*pi/n);
      Q = (cos (m) + sin (m)) / sqrt (n);

    case (-1)
      ##  extrema of T(n-1)
      m = (0:n-1)'*(0:n-1) * (pi/(n-1));
      Q = cos (m);

    case (-2)
      ## zeros of T(n)
      m = (0:n-1)'*(.5:n-.5) * (pi/n);
      Q = cos (m);

    otherwise
      error("gallery: unknown K `%d' for orthog matrix.", k);
  endswitch
endfunction

function A = parter (n)
  ## PARTER    Parter matrix - a Toeplitz matrix with singular values near PI.
  ##           PARTER(N) is the matrix with (i,j) element 1/(i-j+0.5).
  ##           It is a Cauchy matrix and a Toeplitz matrix.
  ##
  ##           At the Second SIAM Conference on Linear Algebra, Raleigh, N.C.,
  ##           1985, Cleve Moler noted that most of the singular values of
  ##           PARTER(N) are very close to PI.  An explanation of the phenomenon
  ##           was given by Parter; see also the paper by Tyrtyshnikov.
  ## 
  ##           References:
  ##           The MathWorks Newsletter, Volume 1, Issue 1, March 1986, page 2.
  ##           S.V. Parter, On the distribution of the singular values of Toeplitz
  ##                matrices, Linear Algebra and Appl., 80 (1986), pp. 115-130.
  ##           E.E. Tyrtyshnikov, Cauchy-Toeplitz matrices and some applications,
  ##                Linear Algebra and Appl., 149 (1991), pp. 1-18.

  if (nargin != 1)
    error ("gallery: 1 argument is required for parter matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for parter matrix.");
  endif

  A = cauchy ((1:n) + 0.5, -(1:n));
endfunction

function P = pei (n, alpha = 1)
  ## PEI    Pei matrix.
  ##        PEI(N, ALPHA), where ALPHA is a scalar, is the symmetric matrix
  ##        ALPHA*EYE(N) + ONES(N).
  ##        If ALPHA is omitted then ALPHA = 1 is used.
  ##        The matrix is singular for ALPHA = 0, -N.
  ##
  ##        Reference:
  ##        M.L. Pei, A test matrix for inversion procedures,
  ##        Comm. ACM, 5 (1962), p. 508.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for pei matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for pei matrix.");
  elseif (! isnumeric (w) || ! isscalar (w))
    error ("gallery: ALPHA must be a scalar for pei matrix.");
  endif

  P = alpha * eye (n) + ones (n);
endfunction

function A = poisson (n)
  ## POISSON   Block tridiagonal matrix from Poisson's equation (sparse).
  ##           POISSON(N) is the block tridiagonal matrix of order N^2
  ##           resulting from discretizing Poisson's equation with the
  ##           5-point operator on an N-by-N mesh.
  ##
  ##           Reference:
  ##           G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
  ##           Johns Hopkins University Press, Baltimore, Maryland, 1989
  ##           (Section 4.5.4).

  if (nargin != 1)
    error ("gallery: 1 argument is required for poisson matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for poisson matrix.");
  endif

  S = tridiag (n, -1, 2, -1);
  I = speye (n);
  A = kron (I, S) + kron (S, I);
endfunction

function A = prolate (n, w = 0.25)
  ## PROLATE   Prolate matrix - symmetric, ill-conditioned Toeplitz matrix.
  ##           A = PROLATE(N, W) is the N-by-N prolate matrix with parameter W.
  ##           It is a symmetric Toeplitz matrix.
  ##           If 0 < W < 0.5 then
  ##              - A is positive definite
  ##              - the eigenvalues of A are distinct, lie in (0, 1), and
  ##                tend to cluster around 0 and 1.
  ##           W defaults to 0.25.
  ##
  ##           Reference:
  ##           J.M. Varah. The Prolate matrix. Linear Algebra and Appl.,
  ##           187:269--278, 1993.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for prolate matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for prolate matrix.");
  elseif (! isnumeric (w) || ! isscalar (w))
    error ("gallery: W must be a scalar for prolate matrix.");
  endif

  a      = zeros (n, 1);
  a(1)   = 2*w;
  a(2:n) = sin (2*pi*w*(1:n-1)) ./ (pi*(1:n-1));

  A = toeplitz(a);
endfunction

function H = randhess (x)
  ## NOTE: this function was named ohess in the original Test Matrix Toolbox
  ## RANDHESS  Random, orthogonal upper Hessenberg matrix.
  ##        H = RANDHESS(N) is an N-by-N real, random, orthogonal
  ##        upper Hessenberg matrix.
  ##        Alternatively, H = RANDHESS(X), where X is an arbitrary real
  ##        N-vector (N > 1) constructs H non-randomly using the elements
  ##        of X as parameters.
  ##        In both cases H is constructed via a product of N-1 Givens rotations.
  ##
  ##        Note: See Gragg (1986) for how to represent an N-by-N (complex)
  ##        unitary Hessenberg matrix with positive subdiagonal elements in terms
  ##        of 2N-1 real parameters (the Schur parametrization).
  ##        This M-file handles the real case only and is intended simply as a
  ##        convenient way to generate random or non-random orthogonal Hessenberg
  ##        matrices.
  ## 
  ##        Reference:
  ##        W.B. Gragg, The QR algorithm for unitary Hessenberg matrices,
  ##        J. Comp. Appl. Math., 16 (1986), pp. 1-8.

  if (nargin != 1)
    error ("gallery: 1 argument is required for randhess matrix.");
  elseif (! isnumeric (x) || ! isreal (x))
    error ("gallery: N or X must be numeric real values for randhess matrix.");
  endif

  if (isscalar (x))
    n = x;
    x = rand (n-1, 1) * 2*pi;
    H = eye (n);
    H(n,n) = sign (randn);
  elseif (isvector (x))
    n = numel (x);
    H = eye (n);
    H(n,n) = sign (x(n)) + (x(n) == 0); # Second term ensures H(n,n) nonzero.
  else
    error ("gallery: N or X must be a scalar or a vector for randhess matrix.");
  endif

  for i = n:-1:2
    ## Apply Givens rotation through angle x(i-1).
    theta = x(i-1);
    c = cos (theta);
    s = sin (theta);
    H([i-1 i], :) = [ c*H(i-1,:)+s*H(i,:)
                     -s*H(i-1,:)+c*H(i,:) ];
  endfor
endfunction

function A = rando (n, k = 1)
  ## RANDO   Random matrix with elements -1, 0 or 1.
  ##         A = RANDO(N, K) is a random N-by-N matrix with elements from
  ##         one of the following discrete distributions (default K = 1):
  ##           K = 1:  A(i,j) =  0 or 1    with equal probability,
  ##           K = 2:  A(i,j) = -1 or 1    with equal probability,
  ##           K = 3:  A(i,j) = -1, 0 or 1 with equal probability.
  ##         N may be a 2-vector, in which case the matrix is N(1)-by-N(2).

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for rando matrix.");
  elseif (! isnumeric (n) || all (numel (n) != [1 2]) || fix (n) != n)
    error ("gallery: N must be an integer for rando matrix.");
  elseif (! isnumeric (k) || ! isscalar (k))
    error ("gallery: K must be a numeric scalar for smoke matrix.");
  endif

  ## Parameter n specifies dimension: m-by-n.
  m = n(1);
  n = n(end);

  switch (k)
    case (1), A =   floor (  rand(m, n) + 0.5);     # {0, 1}
    case (2), A = 2*floor (  rand(m, n) + 0.5) -1;  # {-1, 1}
    case (3), A =   round (3*rand(m, n) - 1.5);     # {-1, 0, 1}
    otherwise
      error ("gallery: unknown K `%d' for smoke matrix.", k);
  endswitch

endfunction

function A = randsvd (n, kappa = sqrt (1/eps), mode = 3, kl = n-1, ku = kl)
  ## RANDSVD  Random matrix with pre-assigned singular values.
  ##       RANDSVD(N, KAPPA, MODE, KL, KU) is a (banded) random matrix of order N
  ##       with COND(A) = KAPPA and singular values from the distribution MODE.
  ##       N may be a 2-vector, in which case the matrix is N(1)-by-N(2).
  ##       Available types:
  ##              MODE = 1:   one large singular value,
  ##              MODE = 2:   one small singular value,
  ##              MODE = 3:   geometrically distributed singular values,
  ##              MODE = 4:   arithmetically distributed singular values,
  ##              MODE = 5:   random singular values with unif. dist. logarithm.
  ##       If omitted, MODE defaults to 3, and KAPPA defaults to SQRT(1/EPS).
  ##       If MODE < 0 then the effect is as for ABS(MODE) except that in the
  ##       original matrix of singular values the order of the diagonal entries
  ##       is reversed: small to large instead of large to small.
  ##       KL and KU are the lower and upper bandwidths respectively; if they
  ##       are omitted a full matrix is produced.
  ##       If only KL is present, KU defaults to KL.
  ##       Special case: if KAPPA < 0 then a random full symmetric positive
  ##                     definite matrix is produced with COND(A) = -KAPPA and
  ##                     eigenvalues distributed according to MODE.
  ##                     KL and KU, if present, are ignored.
  ##
  ##       Reference:
  ##       N.J. Higham, Accuracy and Stability of Numerical Algorithms,
  ##          Society for Industrial and Applied Mathematics, Philadelphia, PA,
  ##          USA, 1996; sec. 26.3.
  ##
  ##       This routine is similar to the more comprehensive Fortran routine xLATMS
  ##       in the following reference:
  ##       J.W. Demmel and A. McKenney, A test matrix generation suite,
  ##       LAPACK Working Note #9, Courant Institute of Mathematical Sciences,
  ##       New York, 1989.

  if (nargin < 1 || nargin > 5)
    error ("gallery: 1 to 5 arguments are required for randsvd matrix.");
  elseif (! isnumeric (n) || all (numel (n) != [1 2]) || fix (n) != n)
    error ("gallery: N must be a 1 or 2 element integer vector for randsvd matrix.");
  elseif (! isnumeric (kappa) || ! isscalar (kappa))
    error ("gallery: KAPPA must be a numeric scalar for randsvd matrix.");
  elseif (abs (kappa) < 1)
    error ("gallery: KAPPA must larger than or equal to 1 for randsvd matrix.");
  elseif (! isnumeric (mode) || ! isscalar (mode))
    error ("gallery: MODE must be a numeric scalar for randsvd matrix.");
  elseif (! isnumeric (kl) || ! isscalar (kl))
    error ("gallery: KL must be a numeric scalar for randsvd matrix.");
  elseif (! isnumeric (ku) || ! isscalar (ku))
    error ("gallery: KU must be a numeric scalar for randsvd matrix.");
  endif

  posdef = 0;
  if (kappa < 0)
    posdef = 1;
    kappa  = -kappa;
  endif

  ## Parameter n specifies dimension: m-by-n.
  m = n(1);
  n = n(end);
  p = min ([m n]);

  ## If A will be a vector
  if (p == 1)
    A = randn (m, n);
    A = A / norm (A);
    return
  end

  ##  Set up vector sigma of singular values.
  switch (abs (mode))
    case (1)
      sigma = ones (p, 1) ./ kappa;
      sigma(1) = 1;
    case (2)
      sigma = ones (p, 1);
      sigma(p) = 1 / kappa;
    case (3)
      factor = kappa^(-1/(p-1));
      sigma  = factor.^[0:p-1];
    case (4)
      sigma = ones (p, 1) - (0:p-1)'/(p-1)*(1-1/kappa);
    case (5)
      ## In this case cond (A) <= kappa.
      rand ("uniform");
      sigma = exp (-rand (p, 1) * log (kappa));
    otherwise
      error ("gallery: unknown MODE `%d' for randsvd matrix.", mode);
  endswitch

  ##  Convert to diagonal matrix of singular values.
  if (mode < 0)
    sigma = sigma (p:-1:1);
  end
  sigma = diag (sigma);

  if (posdef)
    ## handle case where KAPPA was negative
    Q = qmult (p);
    A = Q' * sigma * Q;
    A = (A + A') / 2;  # Ensure matrix is symmetric.
    return
  endif

  if (m != n)
    ## Expand to m-by-n diagonal matrix
    sigma(m, n) = 0;
  end

  if (kl == 0 && ku == 0)
    ## Diagonal matrix requested - nothing more to do.
    A = sigma;
  else
    ##  A = U*sigma*V, where U, V are random orthogonal matrices from the
    ##  Haar distribution.
    A = qmult (sigma');
    A = qmult (A');

    if (kl < n-1 || ku < n-1)
      ## Bandwidth reduction
      A = bandred (A, kl, ku);
    endif
  endif
endfunction

function A = redheff (n)
  ## REDHEFF    A (0,1) matrix of Redheffer associated with the Riemann hypothesis.
  ##            A = REDHEFF(N) is an N-by-N matrix of 0s and 1s defined by
  ##                A(i,j) = 1 if j = 1 or if i divides j,
  ##                A(i,j) = 0 otherwise.
  ##            It has N - FLOOR(LOG2(N)) - 1 eigenvalues equal to 1,
  ##            a real eigenvalue (the spectral radius) approximately SQRT(N),
  ##            a negative eigenvalue approximately -SQRT(N),
  ##            and the remaining eigenvalues are provably ``small''.
  ##            Barrett and Jarvis (1992) conjecture that
  ##              ``the small eigenvalues all lie inside the unit circle
  ##                ABS(Z) = 1'',
  ##            and a proof of this conjecture, together with a proof that some
  ##            eigenvalue tends to zero as N tends to infinity, would yield
  ##            a new proof of the prime number theorem.
  ##            The Riemann hypothesis is true if and only if
  ##            DET(A) = O( N^(1/2+epsilon) ) for every epsilon > 0
  ##                                              (`!' denotes factorial).
  ##            See also RIEMANN.
  ##
  ##            Reference:
  ##            W.W. Barrett and T.J. Jarvis,
  ##            Spectral Properties of a Matrix of Redheffer,
  ##            Linear Algebra and Appl., 162 (1992), pp. 673-683.

  if (nargin != 1)
    error ("gallery: 1 argument is required for redheff matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for redheff matrix.");
  endif

  i = (1:n)' * ones (1, n);
  A = ! rem (i', i);
  A(:,1) = ones (n, 1);
endfunction

function A = riemann (n)
  ## RIEMANN    A matrix associated with the Riemann hypothesis.
  ##            A = RIEMANN(N) is an N-by-N matrix for which the
  ##            Riemann hypothesis is true if and only if
  ##            DET(A) = O( N! N^(-1/2+epsilon) ) for every epsilon > 0
  ##                                              (`!' denotes factorial).
  ##            A = B(2:N+1, 2:N+1), where
  ##            B(i,j) = i-1 if i divides j and -1 otherwise.
  ##            Properties include, with M = N+1:
  ##               Each eigenvalue E(i) satisfies ABS(E(i)) <= M - 1/M.
  ##               i <= E(i) <= i+1 with at most M-SQRT(M) exceptions.
  ##               All integers in the interval (M/3, M/2] are eigenvalues.
  ##
  ##            See also REDHEFF.
  ##
  ##            Reference:
  ##            F. Roesler, Riemann's hypothesis as an eigenvalue problem,
  ##            Linear Algebra and Appl., 81 (1986), pp. 153-198.

  if (nargin != 1)
    error ("gallery: 1 argument is required for riemann matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for riemann matrix.");
  endif

  n = n+1;
  i = (2:n)' * ones (1, n-1);
  j = i';
  A = i .* (! rem (j, i)) - ones (n-1);
endfunction

function A = ris (n)
  ## NOTE: this function was named dingdong in the original Test Matrix Toolbox
  ## RIS       Dingdong matrix - a symmetric Hankel matrix.
  ##           A = RIS(N) is the symmetric N-by-N Hankel matrix with
  ##                          A(i,j) = 0.5/(N-i-j+1.5).
  ##           The eigenvalues of A cluster around PI/2 and -PI/2.
  ##
  ##           Invented by F.N. Ris.
  ##
  ##           Reference:
  ##           J.C. Nash, Compact Numerical Methods for Computers: Linear
  ##           Algebra and Function Minimisation, second edition, Adam Hilger,
  ##           Bristol, 1990 (Appendix 1).

  if (nargin != 1)
    error ("gallery: 1 argument is required for ris matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for ris matrix.");
  endif

  p = -2*(1:n) + (n+1.5);
  A = cauchy (p);
endfunction

function A = smoke (n, k = 0)
  ## SMOKE     Smoke matrix - complex, with a `smoke ring' pseudospectrum.
  ##           SMOKE(N) is an N-by-N matrix with 1s on the
  ##           superdiagonal, 1 in the (N,1) position, and powers of
  ##           roots of unity along the diagonal.
  ##           SMOKE(N, 1) is the same except for a zero (N,1) element.
  ##           The eigenvalues of SMOKE(N, 1) are the N'th roots of unity;
  ##           those of SMOKE(N) are the N'th roots of unity times 2^(1/N).
  ##
  ##           Try PS(SMOKE(32)).  For SMOKE(N, 1) the pseudospectrum looks
  ##           like a sausage folded back on itself.
  ##           GERSH(SMOKE(N, 1)) is interesting.
  ##
  ##           Reference:
  ##           L. Reichel and L.N. Trefethen, Eigenvalues and pseudo-eigenvalues of
  ##           Toeplitz matrices, Linear Algebra and Appl., 162-164:153-185, 1992.

  if (nargin < 1 || nargin > 2)
    error ("gallery: 1 to 2 arguments are required for smoke matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be an integer for smoke matrix.");
  elseif (! isnumeric (n) || ! isscalar (n))
    error ("gallery: K must be a numeric scalar for smoke matrix.");
  endif

  w = exp(2*pi*i/n);
  A = diag( [w.^(1:n-1) 1] ) + diag(ones(n-1,1),1);

  switch (k)
    case (0), A(n,1) = 1;
    case (1), # do nothing
    otherwise,
      error ("gallery: K must have a value of 0 or 1 for smoke matrix.");
  endswitch
endfunction

function T = toeppd (n, m = n, w = rand (m,1), theta = rand (m,1))
  ## NOTE: this function was named pdtoep in the original Test Matrix Toolbox
  ## TOEPPD   Symmetric positive definite Toeplitz matrix.
  ##          TOEPPD(N, M, W, THETA) is an N-by-N symmetric positive (semi-)
  ##          definite (SPD) Toeplitz matrix, comprised of the sum of M rank 2
  ##          (or, for certain THETA, rank 1) SPD Toeplitz matrices.
  ##          Specifically,
  ##                  T = W(1)*T(THETA(1)) + ... + W(M)*T(THETA(M)),
  ##          where T(THETA(k)) has (i,j) element COS(2*PI*THETA(k)*(i-j)).
  ##          Defaults: M = N, W = RAND(M,1), THETA = RAND(M,1).
  ##
  ##          Reference:
  ##          G. Cybenko and C.F. Van Loan, Computing the minimum eigenvalue of
  ##          a symmetric positive definite Toeplitz matrix, SIAM J. Sci. Stat.
  ##          Comput., 7 (1986), pp. 123-131.

  if (nargin < 1 || nargin > 4)
    error ("gallery: 1 to 4 arguments are required for toeppd matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be a numeric integer for toeppd matrix.");
  elseif (! isnumeric (m) || ! isscalar (m) || fix (m) != m)
    error ("gallery: M must be a numeric integer for toeppd matrix.");
  elseif (numel (w) != m || numel (theta) != m)
    error ("gallery: W and THETA must be vectors of length M for toeppd matrix.");
  endif

  T = zeros (n);
  E = 2*pi * ((1:n)' * ones (1, n) - ones (n, 1) * (1:n));

  for i = 1:m
    T = T + w(i) * cos (theta(i)*E);
  endfor
endfunction

function P = toeppen (n, a = 1, b = -10, c = 0, d = 10, e = 1)
  ## NOTE: this function was named pentoep in the original Test Matrix Toolbox
  ## TOEPPEN   Pentadiagonal Toeplitz matrix (sparse).
  ##           P = TOEPPEN(N, A, B, C, D, E) is the N-by-N pentadiagonal
  ##           Toeplitz matrix with diagonals composed of the numbers
  ##           A =: P(3,1), B =: P(2,1), C =: P(1,1), D =: P(1,2), E =: P(1,3).
  ##           Default: (A,B,C,D,E) = (1,-10,0,10,1) (a matrix of Rutishauser).
  ##                     This matrix has eigenvalues lying approximately on
  ##                     the line segment 2*cos(2*t) + 20*i*sin(t).
  ##
  ##           Interesting plots are
  ##           PS(FULL(TOEPPEN(32,0,1,0,0,1/4)))  - `triangle'
  ##           PS(FULL(TOEPPEN(32,0,1/2,0,0,1)))  - `propeller'
  ##           PS(FULL(TOEPPEN(32,0,1/2,1,1,1)))  - `fish'
  ##
  ##           References:
  ##           R.M. Beam and R.F. Warming, The asymptotic spectra of
  ##              banded Toeplitz and quasi-Toeplitz matrices, SIAM J. Sci.
  ##              Comput. 14 (4), 1993, pp. 971-1006.
  ##           H. Rutishauser, On test matrices, Programmation en Mathematiques
  ##              Numeriques, Editions Centre Nat. Recherche Sci., Paris, 165,
  ##              1966, pp. 349-365.

  if (nargin < 1 || nargin > 6)
    error ("gallery: 1 to 6 arguments are required for toeppen matrix.");
  elseif (! isnumeric (n) || ! isscalar (n) || fix (n) != n)
    error ("gallery: N must be a numeric integer for toeppen matrix.");
  elseif (any (cellfun (@(x) ! isnumeric (x) || ! isscalar (x), {a b c d e})))
    error ("gallery: A, B, C, D and E must be numeric scalars for toeppen matrix.");
  endif

  P = spdiags ([a*ones(n,1) b*ones(n,1) c*ones(n,1) d*ones(n,1) e*ones(n,1)],
                -2:2, n, n);
endfunction

function T = tridiag (n, x = -1, y = 2, z = -1)
  ## TRIDIAG  Tridiagonal matrix (sparse).
  ##          TRIDIAG(X, Y, Z) is the tridiagonal matrix with subdiagonal X,
  ##          diagonal Y, and superdiagonal Z.
  ##          X and Z must be vectors of dimension one less than Y.
  ##          Alternatively TRIDIAG(N, C, D, E), where C, D, and E are all
  ##          scalars, yields the Toeplitz tridiagonal matrix of order N
  ##          with subdiagonal elements C, diagonal elements D, and superdiagonal
  ##          elements E.   This matrix has eigenvalues (Todd 1977)
  ##                   D + 2*SQRT(C*E)*COS(k*PI/(N+1)), k=1:N.
  ##          TRIDIAG(N) is the same as TRIDIAG(N,-1,2,-1), which is
  ##          a symmetric positive definite M-matrix (the negative of the
  ##          second difference matrix).
  ##
  ##          References:
  ##          J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
  ##            Birkhauser, Basel, and Academic Press, New York, 1977, p. 155.
  ##          D.E. Rutherford, Some continuant determinants arising in physics and
  ##            chemistry---II, Proc. Royal Soc. Edin., 63, A (1952), pp. 232-241.

  if (nargin != 1 && nargin != 3 && nargin != 4)
    error ("gallery: 1, 3, or 4 arguments are required for tridiag matrix.");
  elseif (nargin == 3)
    z = y;
    y = x;
    x = n;
  endif

  ## Force column vectors
  x = x(:);
  y = y(:);
  z = z(:);

  if (isscalar (x) && isscalar (y) && isscalar (z))
    x *= ones (n-1, 1);
    z *= ones (n-1, 1);
    y *= ones (n,   1);
  elseif (numel (y) != numel (x) + 1)
    error ("gallery: X must have one element less than Y for tridiag matrix.");
  elseif (numel (y) != numel (z) + 1)
    error ("gallery: Z must have one element less than Y for tridiag matrix.");
  endif

  ##  T = diag (x, -1) + diag (y) + diag (z, 1);  # For non-sparse matrix.
  n = numel (y);
  T = spdiags ([[x;0] y [0;z]], -1:1, n, n);
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
  elseif (! isnumeric (n) || all (numel (n) != [1 2]))
    error("gallery: N must be a 1 or 2 elements vector for triw matrix.")
  elseif (! isscalar (alpha))
    error("gallery: ALPHA must be a scalar for triw matrix.")
  elseif (! isscalar (k) || ! isnumeric (k) || fix (k) != k)
    error("gallery: K must be a numeric integer for triw matrix.")
  endif

  m = n(1);              # Parameter n specifies dimension: m-by-n.
  n = n(end);

  t = tril (eye (m, n) + alpha * triu (ones (m, n), 1), k);
endfunction

function A = wathen (nx, ny, k = 0)
  ## # WATHEN returns the Wathen matrix.
  ##
  ##   Discussion:
  ##
  ##     The Wathen matrix is a finite element matrix which is sparse.
  ##
  ##     The entries of the matrix depend in part on a physical quantity
  ##     related to density.  That density is here assigned random values between
  ##     0 and 100.
  ##
  ##     A = WATHEN ( NX, NY ) is a sparse random N-by-N finite element matrix
  ##     where N = 3*NX*NY + 2*NX + 2*NY + 1.
  ##
  ##     A is the consistent mass matrix for a regular NX-by-NY
  ##     grid of 8-node (serendipity) elements in 2 space dimensions.
  ##
  ##     Here is an illustration for NX = 3, NX = 2:
  ##
  ##      23-24-25-26-27-28-29
  ##       |     |     |     |
  ##      19    20    21    22
  ##       |     |     |     |
  ##      12-13-14-15-16-17-18
  ##       |     |     |     |
  ##       8     9    10    11
  ##       |     |     |     |
  ##       1--2--3--4--5--6--7
  ##
  ##     For this example, the total number of nodes is, as expected,
  ##
  ##       N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29.
  ##
  ##     A is symmetric positive definite for any (positive) values of
  ##     the density, RHO(NX,NY), which is chosen randomly in this routine.
  ##
  ##     In particular, if D = DIAG(DIAG(A)), then
  ##       0.25 <= EIG(INV(D)*A) <= 4.5
  ##     for any positive integers NX and NY and any densities RHO(NX,NY).
  ##
  ##     A = WATHEN ( NX, NY, 1 ) returns the diagonally scaled matrix.
  ##
  ##   Modified:
  ##
  ##     17 September 2007
  ##
  ##   Author:
  ##
  ##     Nicholas Higham
  ##
  ##   Reference:
  ##
  ##     Nicholas Higham,
  ##     Algorithm 694: A Collection of Test Matrices in MATLAB,
  ##     ACM Transactions on Mathematical Software,
  ##     Volume 17, Number 3, September 1991, pages 289-305.
  ##
  ##     Andrew Wathen,
  ##     Realistic eigenvalue bounds for the Galerkin mass matrix,
  ##     IMA Journal of Numerical Analysis,
  ##     Volume 7, 1987, pages 449-457.
  ##
  ##   Parameters:
  ##
  ##     Input, integer NX, NY, the number of elements in the X and Y directions
  ##     of the finite element grid.  NX and NY must each be at least 1.
  ##
  ##     Optional input, integer K, is used to request that the diagonally scaled
  ##     version of the matrix be returned.  This happens if K is specified with
  ##     the value 1.
  ##
  ##     Output, sparse real A(N,N), the matrix.  The dimension N is determined by
  ##     NX and NY, as described above.  A is stored in the MATLAB sparse matrix
  ##     format.

  if (nargin < 2 || nargin > 3)
    error ("gallery: 2 or 3 arguments are required for wathen matrix.");
  elseif (! isnumeric (nx) || ! isscalar (nx) || nx < 1)
    error ("gallery: NX must be a positive scalar for wathen matrix.");
  elseif (! isnumeric (ny) || ! isscalar (ny) || ny < 1)
    error ("gallery: NY must be a positive scalar for wathen matrix.");
  elseif (! isscalar (k))
    error ("gallery: K must be a scalar for wathen matrix.");
  endif

  e1 = [ 6  -6   2  -8
        -6  32  -6  20
         2  -6   6  -6
        -8  20  -6  32 ];

  e2 = [ 3  -8   2  -6
        -8  16  -8  20
         2  -8   3  -8
        -6  20  -8  16 ];

  e = [ e1  e2
        e2' e1] / 45;

  n = 3*nx*ny + 2*nx + 2*ny + 1;

  A = sparse (n, n);

  rho = 100 * rand (nx, ny);

  for j = 1:ny
    for i = 1:nx
      ## 
      ##   For the element (I,J), determine the indices of the 8 nodes.
      ## 
      nn(1) = 3*j*nx + 2*i + 2*j + 1;
      nn(2) = nn(1) - 1;
      nn(3) = nn(2) - 1;
      nn(4) = (3*j - 1) * nx + 2*j + i - 1;
      nn(5) = 3 * (j-1) * nx + 2*i + 2*j - 3;
      nn(6) = nn(5) + 1;
      nn(7) = nn(6) + 1;
      nn(8) = nn(4) + 1;

      em = e * rho(i,j);

      for krow = 1:8
        for kcol = 1:8
          A(nn(krow),nn(kcol)) = A(nn(krow),nn(kcol)) + em(krow,kcol);
        endfor
      endfor

    endfor
  endfor

  ## If requested, return A with diagonal scaling.
  if (k)
    A = diag (diag (A)) \ A;
  endif
endfunction

function [A, b] = wilk (n)
  ## WILK   Various specific matrices devised/discussed by Wilkinson.
  ##        [A, b] = WILK(N) is the matrix or system of order N.
  ##        N = 3: upper triangular system Ux=b illustrating inaccurate solution.
  ##        N = 4: lower triangular system Lx=b, ill-conditioned.
  ##        N = 5: HILB(6)(1:5,2:6)*1.8144.  Symmetric positive definite.
  ##        N = 21: W21+, tridiagonal.   Eigenvalue problem.
  ##
  ##        References:
  ##        J.H. Wilkinson, Error analysis of direct methods of matrix inversion,
  ##           J. Assoc. Comput. Mach., 8 (1961),  pp. 281-330.
  ##        J.H. Wilkinson, Rounding Errors in Algebraic Processes, Notes on Applied
  ##           Science No. 32, Her Majesty's Stationery Office, London, 1963.
  ##        J.H. Wilkinson, The Algebraic Eigenvalue Problem, Oxford University
  ##           Press, 1965.

  if (nargin != 1)
    error ("gallery: 1 argument is required for wilk matrix.");
  elseif (! isnumeric (n) || ! isscalar (n))
    error ("gallery: N must be a numeric scalar for wilk matrix.");
  endif

  if (n == 3)
    ## Wilkinson (1961) p.323.
    A = [ 1e-10   0.9  -0.4
          0       0.9  -0.4
          0       0     1e-10 ];

    b = [ 0
          0
          1];

  elseif (n == 4)
    ## Wilkinson (1963) p.105.
    A = [0.9143e-4  0          0          0
         0.8762     0.7156e-4  0          0
         0.7943     0.8143     0.9504e-4  0
         0.8017     0.6123     0.7165     0.7123e-4];

    b = [0.6524
         0.3127
         0.4186
         0.7853];

  elseif (n == 5)
    ## Wilkinson (1965), p.234.
    A = hilb (6);
    A = A(1:5, 2:6) * 1.8144;

  elseif (n == 21)
    ## Wilkinson (1965), p.308.
    E = diag (ones (n-1, 1), 1);
    m = (n-1)/2;
    A = diag (abs (-m:m)) + E + E';

  else
    error("gallery: unknown N %d for wilk matrix.", n);
  endif
endfunction

## NOTE: bandred is part of the Test Matrix Toolbox and is used by randsvd()
function A = bandred (A, kl, ku)
  ## BANDRED  Band reduction by two-sided unitary transformations.
  ##          B = BANDRED(A, KL, KU) is a matrix unitarily equivalent to A
  ##          with lower bandwidth KL and upper bandwidth KU
  ##          (i.e. B(i,j) = 0 if i > j+KL or j > i+KU).
  ##          The reduction is performed using Householder transformations.
  ##          If KU is omitted it defaults to KL.
  ##
  ##          Called by RANDSVD.
  ##          This is a `standard' reduction.  Cf. reduction to bidiagonal form
  ##          prior to computing the SVD.  This code is a little wasteful in that
  ##          it computes certain elements which are immediately set to zero!
  ## 
  ##          Reference:
  ##          G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
  ##          Johns Hopkins University Press, Baltimore, Maryland, 1989.
  ##          Section 5.4.3.

  ##  Check for special case where order of left/right transformations matters.
  ##  Easiest approach is to work on the transpose, flipping back at the end.
  flip = false;
  if (ku == 0)
    flip = true;
    A = A';
    [ku, kl] = deal (kl, ku);
  endif

  [m, n] = size (A);

  for j = 1:min (min (m, n), max (m-kl-1, n-ku-1))
    if (j+kl+1 <= m)
      [v, beta] = house (A(j+kl:m,j));
      temp = A(j+kl:m,j:n);
      A(j+kl:m,j:n) = temp - beta*v*(v'*temp);
      A(j+kl+1:m,j) = zeros (m-j-kl, 1);
    endif

    if (j+ku+1 <= n)
      [v, beta] = house (A(j,j+ku:n)');
      temp = A(j:m,j+ku:n);
      A(j:m,j+ku:n) = temp - beta*(temp*v)*v';
      A(j,j+ku+1:n) = zeros (1, n-j-ku);
    endif
  endfor

  if (flip)
    A = A';
  endif
endfunction
