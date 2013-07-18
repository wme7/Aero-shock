 subroutine nsoli(x,fval,tol, ipar, rpar, sol, it_hist, ierr, n, iupar, rupar) !x_hist,n)
!function [sol, it_hist, ierr, x_hist] = nsoli(x,f,tol, parms)
! NSOLI  Newton-Krylov solver, globally convergent 
!        solver for f(x) = 0
!
! Inexact-Newton-Armijo iteration
!
! Eisenstat-Walker forcing term
!
! Parabolic line search via three point interpolation.
!
! C. T. Kelley, April 27, 2001
!
! This code comes with no guarantee or warranty of any kind.
!
! function [sol, it_hist, ierr, x_hist] = nsoli(x,f,tol,parms)
!
! inputs:
!        initial iterate = x
!        function = f
!        tol = [atol, rtol] relative/absolute
!            error tolerances for the nonlinear iteration
!        parms = [maxit, maxitl, etamax, lmeth, restart_limit]
! ipar(1) 
!            maxit = maxmium number of nonlinear iterations
!                default = 40
! ipar(2)
!            maxitl = maximum number of inner iterations before restart
!                in GMRES(m), m = maxitl 
!                default = 40
!                
!                For iterative methods other than GMRES(m) maxitl
!                is the upper bound on linear iterations.
! !rpar
!            |etamax| = Maximum error tolerance for residual in inner
!                iteration. The inner iteration terminates
!                when the relative linear residual is
!                smaller than eta*| F(x_c) |. eta is determined
!                by the modified Eisenstat-Walker formula if etamax > 0.
!                If etamax < 0, then eta = |etamax| for the entire
!                iteration.
!                default: etamax = .9
! ipar(3)
!            lmeth = choice of linear iterative method
!                    1 (GMRES), 2 GMRES(m), 
!                    3 (BICGSTAB), 4 (TFQMR)
!                 default = 1 (GMRES, no restarts)
! ipar(4) 
!            restart_limit = max number of restarts for GMRES if
!                    lmeth = 2
!                  default = 20
!
! output:
!        sol = solution
!        it_hist(maxit,3) = l2 norms of nonlinear residuals
!            for the iteration, number of function evaluations,
!            and number of steplength reductions
!        ierr = 0 upon successful termination
!        ierr = 1 if after maxit iterations
!             the termination criterion is not satsified
!        ierr = 2 failure in the line search. The iteration
!             is terminated if too many steplength reductions
!             are taken.
!
!    x_hist = matrix of the entire interation history.
!             The columns are the nonlinear iterates. This
!             is useful for making movies, for example, but
!             can consume way too much storage. This is an
!             OPTIONAL argument. Storage is only allocated
!             if x_hist is in the output argument list.
!
!
!
! internal parameters:
!       debug = turns on/off iteration statistics display as
!               the iteration progresses
!
!       alpha = 1.d-4, parameter to measure sufficient decrease
!
!       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
!
!       maxarm = 20, maximum number of steplength reductions before
!                    failure is reported
!
implicit none
integer:: n
real(kind=8):: x(n)
external :: fval
real(kind=8), intent(in)::tol(2) 

integer:: ipar(4) 
real(kind=8):: rpar
integer:: iupar(*)
real(kind=8):: rupar(*)
 

real(kind=8)::sol(n)

real(kind=8)::it_hist(ipar(1),3)
real(kind=8):: it_histx(ipar(1),3)

integer:: ierr
!, x_hist

! Functions
real(kind=8):: norm

!Local Variables
integer:: debug, maxarm
integer:: maxit, lmaxit, lmeth, itc
integer:: iarm, restart_limit
integer:: igmparms(3)
integer:: inner_it_count, inner_f_evals
real(kind=8):: alpha, sigma0, sigma1, gamma
real(kind=8):: rgmparms, etamax, etaold, etanew, errstep
real(kind=8):: rtol, atol, fnrm
real(kind=8):: stop_tol, fnrmo, rat
real(kind=8):: lambda, lamm, lamc
real(kind=8):: nft, nf0, ff0, ffc, ffm
real(kind=8):: xold(n), xt(n), step(n), ft(n), f0(n)
real(kind=8):: outstat(ipar(1),4)

real(kind=8):: parab3p
! Set the debug parameter 1 turns display on, otherwise off.

debug = 0

! Set internal parameters.

alpha = 1.d-4 
sigma0 = .1d0
sigma1 = .5d0 
maxarm = 20 
gamma = .9d0

! Initialize it_hist, ierr, x_hist, and set the default values of
! those iteration parameters which are optional inputs.

ierr = 0 
!it_histx = zeros(maxit,3)

!if nargout == 4, x_hist = x end

! Initialize parameters for the iterative methods.
! Check for optional inputs.

!gmparms = [abs(etamax), lmaxit]
if (ipar(1) .gt. 0) then
    maxit = ipar(1)
else
    maxit = 40
endif
if (ipar(2) .gt. 0) then
    lmaxit = ipar(2) 
else
    lmaxit = 40
endif
if (ipar(3) .gt. 0) then
    lmeth = ipar(3)
else
    lmeth = 1 !default GMRES
endif
if (ipar(4) .gt. 0) then
    restart_limit = ipar(4)
else
    restart_limit = 20
endif

!if (rpar .gt. 0d0) then
    etamax = rpar
!else
!    etamax = .9d0
!endif

    rgmparms = abs(etamax)
!    if length(parms) == 5
       igmparms = (/lmaxit, ipar(4), 1/)
!    end

!
rtol = tol(2) 
atol = tol(1) 
!n = length(x) 
fnrm = 1d0 
itc = 0


! Evaluate f at the initial iterate,and
! compute the stop tolerance.
!
! Call subroutine
call fval(x,f0, iupar, rupar)
fnrm = norm(f0,n)

it_histx(itc+1,1) = fnrm
it_histx(itc+1,2) = 0d0 
it_histx(itc+1,3) = 0d0

fnrmo = 1
stop_tol = atol + rtol*fnrm
outstat(itc+1, :) = (/fnrm, 0, 0, 0 /)

! main iteration loop

do while ((fnrm > stop_tol) .and. (itc < maxit))

! Keep track of the ratio (rat = fnrm/frnmo)
! of successive residual norms and 
! the iteration counter (itc).

    rat = fnrm/fnrmo
    fnrmo = fnrm 
    itc = itc+1

    call dkrylov(f0, fval, x, igmparms, rgmparms, lmeth, step, errstep, &
      inner_it_count,inner_f_evals,n, iupar, rupar)
!    [step, errstep, inner_it_count,inner_f_evals] = ...
!         dkrylov(f0, f, x, gmparms, lmeth)

!   The line search starts here.

    xold = x
    lambda = 1 
    lamm = 1 
    lamc = lambda 
    iarm = 0

    xt = x + lambda*step
    call fval(xt,ft, iupar, rupar)

    nft = norm(ft,n) 
    nf0 = norm(f0,n) 
    ff0 = nf0*nf0 
    ffc = nft*nft 
    ffm = nft*nft
   do  while (nft >= (1d0 - alpha*lambda) * nf0)

!   Apply the three point parabolic model.

        if (iarm == 0) then
            lambda = sigma1*lambda 
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm) 
        endif

! Update x keep the books on lambda.

        xt = x+lambda*step
        lamm = lamc
        lamc = lambda

! Keep the books on the function norms.

        call fval(xt,ft, iupar, rupar)
        nft = norm(ft,n)
        ffm = ffc
        ffc = nft*nft
        iarm = iarm+1
        if (iarm > maxarm) then
            write(6,*) "Armijo failure, too many reductions"
            ierr = 2
            write(6,*) "outstat",outstat
            it_hist = it_histx(1:itc+1,:)
      !!!  if nargout == 4, x_hist = [x_hist,x] end
            sol = xold
            return
        endif
    end do
    x = xt
    f0 = ft

!   End of line search.

!!!    if nargout == 4, x_hist = [x_hist,x] end

    fnrm = norm(f0,n)
    it_histx(itc+1,1) = fnrm 

!   How many function evaluations did this iteration require?

    it_histx(itc+1,2) = it_histx(itc,2)+inner_f_evals+iarm+1

    if (itc == 1) then
       it_histx(itc+1,2) = it_histx(itc+1,2)+1 
    endif
    it_histx(itc+1,3) = iarm
!
    rat = fnrm/fnrmo

!   Adjust eta as per Eisenstat-Walker.

    if (etamax > 0) then
        etaold = rgmparms
        etanew = gamma*rat*rat
        if (gamma*etaold*etaold > .1d0) then
            etanew = max(etanew,gamma*etaold*etaold)
        endif

        rgmparms = min(etanew,etamax)
        rgmparms = max(rgmparms, .5d0*stop_tol/fnrm)
    endif
!
    outstat(itc+1, :) = (/fnrm,dble(inner_it_count),rat,dble(iarm)/)
!
enddo !main

sol = x
it_hist = it_histx(1:itc+1,:)

if (debug == 1) then
    write(6,*) "outstat", outstat
    it_hist = it_histx(1:itc+1,:)
endif

! on failure, set the error flag

        call fval(x,ft, iupar, rupar)
        fnrm = norm(ft,n)

if (fnrm > stop_tol) then
   ierr = 1 
endif
        call fval(x,ft, iupar, rupar)
        nft = norm(ft,n)

ierr=0
return
!
end subroutine nsoli

real(kind=8) function parab3p(lambdac, lambdam, ff0, ffc, ffm)
!function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
! Apply three-point safeguarded parabolic model for a line search.
!
! C. T. Kelley, April 1, 2003
!
! This code comes with no guarantee or warranty of any kind.
!
! function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
!
! input:
!       lambdac = current steplength
!       lambdam = previous steplength
!       ff0 = value of \| F(x_c) \|^2
!       ffc = value of \| F(x_c + \lambdac d) \|^2
!       ffm = value of \| F(x_c + \lambdam d) \|^2
!
! output:
!       lambdap = new value of lambda given parabolic model
!
! internal parameters:
!       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
!
implicit none

real(kind=8):: lambdac, lambdam, ff0, ffc, ffm

! internal parameters:

real(kind=8):: sigma0, sigma1
real(kind=8):: c1, c2, lambdap
! set internal parameters

 sigma0 = .1d0
 sigma1 = .5d0
!
! compute coefficients of interpolation polynomial
!
! p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
!
! d1 = (lambdac - lambdam)*lambdac*lambdam < 0
!      so if c2 > 0 we have negative curvature and default to
!      lambdap = sigam1 * lambda
!
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0)
if (c2 >= 0) then
    lambdap = sigma1*lambdac 
    parab3p = lambdap
    return
endif
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0)
lambdap = -c1*.5d0/c2
if (lambdap < sigma0*lambdac) then
 lambdap = sigma0*lambdac 
endif
if (lambdap > sigma1*lambdac) then
 lambdap = sigma1*lambdac 
endif

parab3p = lambdap
return
end
!

subroutine dkrylov(f0, fval, x, ipar,rpar, lmeth, step, errstep, &
total_iters, f_evals,n, iupar, rupar)
!function [step, errstep, total_iters, f_evals] = ...
!    dkrylov(f0, f, x, params, lmeth)

! Krylov linear equation solver for use in nsoli
!
! C. T. Kelley, April 1, 2003
!
!
! This code comes with no guarantee or warranty of any kind.
!
! function [step, errstep, total_iters, f_evals] 
!                              = dkrylov(f0, f, x, params, lmeth)
!
!
! Input:  f0 = function at current point
!         f = nonlinear function
!              the format for f is  function fx = f(x)
!              Note that for Newton-GMRES we incorporate any
!              preconditioning into the function routine.
!         x = current point
!         params = vector to control iteration
! rpar
!              params(1) = relative residual reduction factor
! ipar
!              params(2) = max number of iterations
!              params(3) = max number of restarts for GMRES(m)
!              params(4) (Optional) = reorthogonalization method in GMRES
!                   1 -- Brown/Hindmarsh condition (default)
!                   2 -- Never reorthogonalize (not recommended)
!                   3 -- Always reorthogonalize (not cheap!)
!
!         lmeth = method choice
!              1 GMRES without restarts (default)
!              2 GMRES(m), m = params(2) and the maximum number
!                   of restarts is params(3) 
!              3 Bi-CGSTAB
!              4 TFQMR
!
! Output: x = solution
!         errstep = vector of residual norms for the history of
!                 the iteration
!         total_iters = number of iterations
!
!
implicit none

integer:: n, ipar(3), lmeth
real(kind=8):: f0(n), x(n), step(n), xinit(n)
real(kind=8):: rpar
real(kind=8):: errstep(n)
integer, intent(out):: total_iters, f_evals
external:: fval
integer:: iupar(*)
real(kind=8):: rupar(*)

integer:: restart_limit, lmaxit, kinn
integer::  igmparms(3)
real(kind=8):: rgmparms
real(kind=8):: norm
! initialization

lmaxit = ipar(1)
!lmaxit = params(2)

!restart_limit = 20 !default
restart_limit = ipar(2)

!if (length(params) >= 3) then
!    restart_limit = params(3)
!endif

if (lmeth == 1) then
  restart_limit = 0 
endif

!if length(params) == 3
!! default reorthogonalization
!     gmparms = [params(1), params(2), 1]
!     rgmparms = rpar
!     igmparms =(/ ipar(1), ipar(2), 1 /)
!elseif (length(params) == 4)
!! reorthogonalization method is params(4)
!     gmparms = [params(1), params(2), params(4)]
!else
!     gmparms = [params(1), params(2)]
!endif
     rgmparms = rpar
     igmparms =(/ ipar(1), 1, 0 /)

! lmaxit, default reorthogonalization method, initialize xinit
   
! linear iterative methods

if ((lmeth == 1) .or. (lmeth == 2)) then ! GMRES or GMRES(m) 

! compute the step using a GMRES routine especially designed
! for this purpose

!    [step, errstep, total_iters] = dgmres(f0, f, x, gmparms)
    xinit=0d0
    call dgmres(f0, fval, x, igmparms, rgmparms, xinit, step, errstep, &
total_iters,n, iupar, rupar)
    kinn = 0
    igmparms(3) =1 
!   restart at most restart_limit times
    do while ( (total_iters == lmaxit) .and.  &
          (errstep(total_iters) > rpar*norm(f0,n)) &
          .and. (kinn < restart_limit)) 
        kinn = kinn+1
        xinit=step
        call dgmres(f0, fval, x, igmparms, rgmparms, xinit, step, errstep, &
total_iters,n, iupar, rupar)
!        [step, errstep, total_iters] = dgmres(f0, f, x, gmparms,step)
    enddo
    total_iters = total_iters+kinn*lmaxit
    f_evals = total_iters+kinn

! Bi-CGSTAB

elseif (lmeth == 3) then
    call dcgstab(f0, fval, x, igmparms, rgmparms, xinit, step, errstep, &
total_iters,n, iupar, rupar)
!    [step, errstep, total_iters] = dcgstab(f0, f, x, gmparms)
    f_evals = 2*total_iters

! TFQMR

elseif (lmeth == 4) then 
    call dtfqmr(f0, fval, x, igmparms, rgmparms,xinit, step, errstep, &
total_iters,n, iupar, rupar)
!    [step, errstep, total_iters] = dtfqmr(f0, f, x, gmparms)
    f_evals = 2*total_iters
else
    write(6,*) "lmeth error in fdkrylov"
    stop
endif
return
end subroutine dkrylov

subroutine dirder(x,w,fval,f0,n,z, iupar, rupar)
!function z = dirder(x,w,f,f0)
! Finite difference directional derivative
! Approximate f'(x) w
! 
! C. T. Kelley, April 1, 2003
!
! This code comes with no guarantee or warranty of any kind.
!
! function z = dirder(x,w,f,f0)
!
! inputs:
!           x, w = point and direction
!           f = function
!           f0 = f(x), in nonlinear iterations
!                f(x) has usually been computed
!                before the call to dirder
implicit none
integer::n
real(kind=8):: x(n), w(n), f0(n)
integer:: iupar(*)
real(kind=8):: rupar(*)

external:: fval

real(kind=8):: z(n)
real(kind=8):: f1(n), del(n)
real(kind=8):: epsnew,xs
real(kind=8):: norm, sgn
! Use a hardwired difference increment.

epsnew = 1.d-8

!n = length(x)

! scale the step

if (norm(w,n) == 0) then
    z = 0d0
return
endif

! Now scale the difference increment.

!xs=(x'*w)/norm(w)
xs=sum(x*w)/norm(w,n)

if (xs .ne. 0.d0) then
     epsnew=epsnew*max(abs(xs),1.d0)*sgn(xs)
endif
epsnew=epsnew/norm(w,n)

! del and f1 could share the same space if storage
! is more important than clarity.

del = x+epsnew*w
call fval(del,f1, iupar, rupar) !feval
z = (f1 - f0)/epsnew

return
end subroutine dirder


subroutine dgmres(f0, fval, xc, ipar, rpar, xinit, x, error, total_iters,n, iupar, rupar)
use LinearSolvers
implicit none
!function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit)
! GMRES linear equation solver for use in Newton-GMRES solver
!
! C. T. Kelley, April 1, 2003
!
! This code comes with no guarantee or warranty of any kind.
!
! function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit)
!
!
! Input:  f0 = function at current point
!         f = nonlinear function
!              the format for f is  function fx = f(x)
!              Note that for Newton-GMRES we incorporate any 
!              preconditioning into the function routine.
!         xc = current point
!         params = two dimensional vector to control iteration
! rpar
!              params(1) = relative residual reduction factor
! ipar(1)
!              params(2) = max number of iterations 
! ipar(2)           params(3) (Optional) = reorthogonalization method
!                   1 -- Brown/Hindmarsh condition (default)
!                   2 -- Never reorthogonalize (not recommended)
!                   3 -- Always reorthogonalize (not cheap!)
! ipar(3) 0: xinit=0; 1: restart, xinit=step
!         xinit = initial iterate. xinit = 0 is the default. This
!              is a reasonable choice unless restarted GMRES
!              will be used as the linear solver.
!
! Output: x = solution
!         error = vector of residual norms for the history of
!                 the iteration
!         total_iters = number of iterations
!
! Requires givapp.m, dirder.m 

integer:: ipar(3), n
real(kind=8):: rpar
real(kind=8):: f0(n), xc(n), xinit(n)
external:: fval
integer:: iupar(*)
real(kind=8):: rupar(*)

real(kind=8):: norm !, giveapp
real(kind=8):: x(n), error(ipar(1))
integer:: total_iters

integer:: kmax, reorth, k, j
real(kind=8):: r(n), z(n)
real(kind=8):: errtol, b(n)
real(kind=8):: h(ipar(1),ipar(1)), ho(ipar(1),ipar(1))
real(kind=8):: v(n,ipar(1))
real(kind=8):: c(ipar(1)+1),  s(ipar(1)+1)
real(kind=8):: rho, beta, normav2, hr, normav, nu
real(kind=8):: g(ipar(1)+1), go(ipar(1)+1), y(ipar(1))

logical successful   ! Status of computations

! initialization


errtol = rpar !params(1)
kmax = ipar(1) !params(2)
reorth = ipar(2)
!if length(params) == 3
!    reorth = params(3)
!end

! The right side of the linear equation for the step is -f0. 

b = -f0
!n = length(b)

! Use zero vector as initial iterate for Newton step unless
! the calling routine has a better idea (useful for GMRES(m)).

x = 0d0
r = b

if (ipar(3) .gt. 0) then !nargin == 5
    x = xinit
    call dirder(xc, x, fval,f0,n,r, iupar, rupar)
    r = -r-f0
!    r = -dirder(xc, x, fval, f0)-f0
endif

!
h=0d0
v=0d0

! h = zeros(kmax), v = zeros(n,kmax), c = zeros(kmax+1,1), s = zeros(kmax+1,1)
g=0d0
rho = norm(r,n)
g(1) = rho !*eye(kmax+1,1)
errtol = errtol*norm(b,n)

! Test for termination on entry.

error(1) = rho ! [error,rho]
total_iters = 0

if(rho < errtol) then
write(6,*) "early termination"
return
endif

!
v(:,1) = r/rho
beta = rho
k = 0

! GMRES iteration

do while ((rho > errtol) .and. (k < kmax))
    k = k+1

!   Call directional derivative function.
    call dirder(xc, v(:,k), fval, f0, n, r, iupar, rupar)
    v(:,k+1) = r
!    v(:,k+1) = dirder(xc, v(:,k), f, f0)
    normav = norm(v(:,k+1),n)

!   Modified Gram-Schmidt

    do j = 1,k
        h(j,k) = DOT_PRODUCT(v(:,j),v(:,k+1))
        v(:,k+1) = v(:,k+1)-h(j,k)*v(:,j)
    enddo

    h(k+1,k) = dsqrt(DOT_PRODUCT(v(:,k+1),v(:,k+1)))
!    h(k+1,k) = norm(v(:,k+1))
    normav2 = h(k+1,k)

!   Reorthogonalize?

if  ( ((reorth == 1) .and. (normav + .001d0*normav2 == normav)) &
  .or. (reorth ==  3)) then

    do j = 1,k
        hr=DOT_PRODUCT(v(:,j),v(:,k+1))
!        hr = v(:,j)'*v(:,k+1)
	h(j,k) = h(j,k)+hr
        v(:,k+1) = v(:,k+1)-hr*v(:,j)
    enddo
    h(k+1,k) = norm(v(:,k+1),n)
endif

!   Watch out for happy breakdown.

    if(h(k+1,k) .ne. 0) then
    v(:,k+1) = v(:,k+1)/h(k+1,k)
    endif

!   Form and store the information for the new Givens rotation.

    if (k > 1) then
        call givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1,  ho(1:k,k))
         h(1:k,k) = ho(1:k,k)
!        h(1:k,k) = givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1)
    endif

!   Don't divide by zero if solution has  been found.

    nu = norm(h(k:k+1,k),2)
    if (nu .ne. 0) then
        c(k) = h(k,k)/nu
!        c(k) = conjg(h(k,k)/nu)

        s(k) = -h(k+1,k)/nu
        h(k,k) = c(k)*h(k,k)-s(k)*h(k+1,k)
        h(k+1,k) = 0d0
        call givapp(c(k),s(k),g(k:k+1),1,go(k:k+1))
         g(k:k+1) = go(k:k+1)
        !g(k:k+1) = givapp(c(k),s(k),g(k:k+1),1)
    endif

!   Update the residual norm.

    rho = abs(g(k+1))
    error(k) = rho ![error,rho]

!   end of the main while loop

enddo

! At this point either k > kmax or rho < errtol.
! It's time to compute x and leave.

!!! Solver
successful = gaussianElimination( h(1:k,1:k), g(1:k), y )
!y = h(1:k,1:k)\g(1:k)
if ( successful ) then
total_iters = k
x = x + matmul(v(1:n,1:k),y)
else
write(6,*) "Fail to solve the linear system"
endif
return
end subroutine dgmres 

!real(kind=8) function vrot = 
subroutine givapp(c,s,vin,k, vrot)
!  Apply a sequence of k Givens rotations, used within gmres codes.
! 
!  C. T. Kelley, April 1, 2003
!
! This code comes with no guarantee or warranty of any kind.
!
!  function vrot = givapp(c, s, vin, k)
!
implicit none
integer:: k
real(kind=8):: c(k), s(k), vin(k+1)
integer:: i
real(kind=8):: vrot(k+1)
real(kind=8):: w1,w2
vrot = vin
do i = 1,k
    w1 = c(i)*vrot(i)-s(i)*vrot(i+1)
!
!   Here's a modest change that makes the code work in complex
!   arithmetic. Thanks to Howard Elman for this.
!
    w2 = s(i)*vrot(i)+c(i)*vrot(i+1)
!    w2 = s(i)*vrot(i)+conj(c(i))*vrot(i+1)
    vrot(i:i+1) = (/w1, w2 /) ![w1,w2]
enddo
!
return
end subroutine givapp

subroutine dcgstab(f0, f, xc, ipar, rpar, xinit, x, error, total_iters,n)
!function [x, error, total_iters] = ...
!                     dcgstab(f0, f, xc, params, xinit)
! Forward difference Bi-CGSTAB solver for use in nsoli
!
! C. T. Kelley, April 1, 2003
!
! This code comes with no guarantee or warranty of any kind.
!
! function [x, error, total_iters]
!                    = dcgstab(f0, f, xc, params, xinit)
!
! Input:  f0 = function at current point
!         f = nonlinear function
!              the format for f is  function fx = f(x)
!              Note that for Newton-GMRES we incorporate any
!              preconditioning into the function routine.
!         xc = current point
!         params = two dimensional vector to control iteration
! rpar
!              params(1) = relative residual reduction factor
! ipar(1)
!              params(2) = max number of iterations
! ipar(2) =0: xinit = 0; 1: xinit = step 
!         xinit = initial iterate. xinit = 0 is the default. This
!              is a reasonable choice unless restarts are needed.
!
!
! Output: x = solution
!         error = vector of residual norms for the history of
!                 the iteration
!         total_iters = number of iterations
!
! Requires: dirder.m
!
implicit none
integer:: ipar(3), n
real(kind=8):: rpar
real(kind=8):: f0(n), xc(n) 
real(kind=8):: xinit(n)
external:: f

real(kind=8):: norm !, giveapp
real(kind=8):: x(n), error(ipar(1))
integer:: total_iters

integer:: kmax, reorth, k, j
real(kind=8):: r(n), z(n), hatr0(n)
real(kind=8):: errtol, b(n)
real(kind=8):: v(n), p(n), s(n), t(n)

real(kind=8):: alpha, beta, omega, zeta, tau

real(kind=8):: rho(ipar(1)+1) !, beta, normav2, hr, normav, nu



!
! initialization
!
b = -f0 
!n = length(b) 
errtol = rpar*norm(b,n) 
kmax = ipar(1)
rho = 0d0

! Use zero vector as initial iterate for Newton step unless
! the calling routine has a better idea (useful for GMRES(m)).

x = 0d0 
r = b
 if (ipar(3) .gt. 0) then
    !if nargin == 5
    x = xinit
    call dirder(xc, x, f, f0, n, r)
    r = -r-f0
!    r = -dirder(xc, x, f, f0)-f0
endif

!
hatr0 = r
k = 0 
rho(1) = 1 
alpha = 1 
omega = 1
v = 0d0
p = 0d0

rho(2) = DOT_PRODUCT(hatr0,r)

zeta = norm(r,n) 
error(1) = zeta ![error,zeta]

! Bi-CGSTAB iteration

do while( (zeta > errtol) .and.  (k < kmax))
    k = k+1
    if (omega == 0) then
       write(6,*) "Bi-CGSTAB breakdown, omega = 0"
       stop
    endif
    beta = (rho(k+1)/rho(k))*(alpha/omega)
    p = r+beta*(p - omega*v)

    call dirder(xc, p, f, f0, n, v)
!    v = dirder(xc,p,f,f0)

     tau=DOT_PRODUCT(hatr0,v)
!    tau = hatr0'*v

    if (tau == 0) then
        write(6,*) "Bi-CGSTAB breakdown, tau = 0"
        stop
    endif 
    alpha = rho(k+1)/tau
    s = r-alpha*v 

    call dirder(xc, s, f, f0, n, t)
!    t = dirder(xc,s,f,f0)
    tau = DOT_PRODUCT(t,t)

    if (tau == 0) then
       write(6,*) "Bi-CGSTAB breakdown, t = 0"
     stop
    endif
    omega = DOT_PRODUCT(t,s)/tau 
    rho(k+2) = -omega*DOT_PRODUCT(hatr0,t)
!    rho(k+2) = -omega*(hatr0'*t)
    x = x+alpha*p+omega*s
    r = s-omega*t
    zeta = norm(r,n)
    total_iters = k
    error(k) = zeta ![error, zeta]
enddo

return
end subroutine dcgstab

subroutine dtfqmr(f0, fval, xc, ipar, rpar, xinit, x, error, total_iters,n, iupar, rupar)
!function [x, error, total_iters] = ...
!                     dtfqmr(f0, f, xc, params, xinit)
! Forward difference TFQMR solver for use in nsoli
!
! C. T. Kelley, April 1, 2003
!
! This code comes with no guarantee or warranty of any kind.
!
! function [x, error, total_iters]
!                    = dtfqmr(f0, f, xc, params, xinit)
!
!
!
! Input:  f0 = function at current point
!         f = nonlinear function
!              the format for f is  function fx = f(x)
!              Note that for Newton-GMRES we incorporate any
!              preconditioning into the function routine.
!         xc = current point
!         params = two dimensional vector to control iteration
! rpar
!             params(1) = relative residual reduction factor
! ipar(1)
!             params(2) = max number of iterations
! ipar(3) = 0: xinit = 0; 1: xinit = step
!         xinit = initial iterate. xinit = 0 is the default. This
!              is a reasonable choice unless restarts are needed.
!
!
! Output: x = solution
!         error = vector of residual norms for the history of
!                 the iteration
!         total_iters = number of iterations
!
! Requires: dirder.m
!
implicit none
integer:: ipar(3), n
real(kind=8):: rpar
real(kind=8):: f0(n), xc(n)
real(kind=8):: xinit(n)
external:: fval
integer:: iupar(*)
real(kind=8):: rupar(*)

real(kind=8):: norm 
real(kind=8):: x(n), error(ipar(1))
integer:: total_iters

integer:: kmax, reorth, k, j, m
real(kind=8):: errtol, b(n)

real(kind=8):: u(n,2)
real(kind=8):: y(n,2)
real(kind=8):: w(n), r(n)
real(kind=8):: d(n), v(n)

real(kind=8):: theta, eta, tau, rho, rhon, sigma, c, alpha, beta


!
! initialization
!

b = -f0
errtol = rpar*norm(b,n) 
kmax = ipar(1)
x = 0 !zeros(n,1)
r = b
if (ipar(3) .gt. 0) then
   !if nargin == 5
    x = xinit
    call dirder(xc, x, fval, f0, n, r, iupar, rupar)
    r = -r-f0
!    r = -dirder(xc, x, f, f0)-f0
endif

!
u = 0d0 !zeros(n,2) 
y = 0d0 !zeros(n,2) 
w = r 
y(:,1) = r 
k = 0 

d = 0d0
    call dirder(xc, y(:,1), fval, f0, n, v, iupar, rupar)
!v = dirder(xc, y(:,1),f,f0)
u(:,1) = v
theta = 0 
eta = 0 
tau = norm(r,n) 
error(1) = tau ![error,tau]
rho = tau*tau

! TFQMR iteration

do while( k < kmax)
    k = k+1
    sigma =  DOT_PRODUCT(r,v)
!    sigma = r'*v

    if (sigma == 0) then
        write(6,*) "TFQMR breakdown, sigma = 0"
        stop
    endif
!
    alpha = rho/sigma
!
    do j = 1,2

!   Compute y2 and u2 only if you have to

        if (j == 2) then 
            y(:,2) = y(:,1)-alpha*v
            call dirder(xc, y(:,2), fval, f0, n, u(:,2), iupar, rupar)
!            u(:,2) = dirder(xc, y(:,2),f,f0)
        endif

        m = 2*k-2+j
        w = w-alpha*u(:,j)
        d = y(:,j)+(theta*theta*eta/alpha)*d
        theta = norm(w,n)/tau 
        c = 1d0/dsqrt(1d0+theta*theta)
        tau = tau*theta*c 
        eta = c*c*alpha
        x = x+eta*d
!
!   Try to terminate the iteration at each pass through the loop
!
        if (tau*sqrt(dble(m+1)) <=  errtol) then
            error(k) = tau ! [error, tau]
            total_iters = k
            return
        endif
    enddo !j
!
!
!
    if (rho == 0) then
        write(6,*) "TFQMR breakdown, rho = 0"
        stop
    endif
!
    rhon = DOT_PRODUCT(r,w)
!    rhon = r'*w 
    beta = rhon/rho 
    rho = rhon
    y(:,1) = w + beta*y(:,2)

    call dirder(xc, y(:,1), fval, f0, n, u(:,1), iupar, rupar)
!    u(:,1) = dirder(xc, y(:,1),f,f0)
    v = u(:,1)+beta*(u(:,2)+beta*v)
    error(k) = tau ![error, tau]
    total_iters = k
enddo
return
end subroutine dtfqmr
!
real(kind=8) function sgn(x)
implicit none
real(kind=8)::x
if (x .gt. 0d0) then
 sgn=1d0
elseif  (x .lt. 0d0) then
 sgn=-1d0
else
 sgn=0d0
endif
return
end function sgn

real(kind=8) function norm(x,n)
implicit none
integer:: n

real(kind=8),intent(in):: x(n)

!sum=dot_product(x,x)
!norm=dsqrt(sum)

norm=sqrt(dot_product(x,x))
return
end function norm

real(kind=8) function norm2(x,n)
implicit none
integer:: n
real(kind=8), intent(in):: x(n)
real(kind=8)::sum
!n=size(x)
sum=dot_product(x,x)
norm2=dsqrt(sum)

!norm=sqrt(dot_product(x,x))
return
end function norm2
