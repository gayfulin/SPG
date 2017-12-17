module set_generator_parameters

contains

Subroutine Calc_corr(a_mnl, time, corr)

!*********************************************************************************************************
! PURPOSE
!
! We calculate the temporal correlation in spectral space by explicit assymptotic formula
!
! Author: D Gayfulin, 2013
! Current code owner: D. Gayfulin
! Last update: December, 2013

Real, intent(in) :: a_mnl
Real, intent(in) :: time
Real, intent(out):: corr

!===========Exexcution==========================
corr=(exp(-a_mnl*time))*(1+a_mnl*time+(1./3.)*(a_mnl*time)**2)

End Subroutine Calc_corr

subroutine L05_true_definer(xhalfngrid, yhalfngrid,verthalfngrid, mesh_size, q, lambda, delta, L05_true)

!*********************************************************************************************************
! PURPOSE
!
! We calculate L05 for given lambda, using a theoretical formula
! 
!
! Authors:  I.Mamay, 2012
!           D Gayfulin, 2013 -2014
! Current code owner: D. Gayfulin
! Last update: March, 2014
!**********************************************************************************************************

Use FFT, Only: SET99, FFT991
Implicit None
!==========================External parameters=========
Integer, intent(in)  :: xhalfngrid, yhalfngrid     ! hor. resolution
Integer, intent(in)  :: verthalfngrid ! vert. resolution
Real   , intent(in)  :: lambda, delta       
Real   , intent(in)  :: q            ! smoothness parameter
!======================Local parameters================
Real    :: sigma
Real    :: gamma
Real    :: corr(0:2*xhalfngrid)
Real    :: Dist(0:xhalfngrid)
Real    :: F_sum(0:2*xhalfngrid)   !sum of b(k, l, m) for fixed k
Real   , intent(out) :: L05_true
Integer :: xngrid, weight2
Real :: mesh_size
Integer :: i, j, k, l,m  ! loop counters
Real :: mu               ! temporal pseudoscales
Real :: pi               ! number pi
Real :: factor           ! value ff3(1,1,1)
Real :: Calc_L           ! result of function Calc_L
Real :: P_d              ! function P_d
Real :: b_nml            ! variance of xi_nml
Real :: p                ! value of P_d
Integer  :: l_05         ! lower estimate of L05 value

Integer :: ifax(1:13)                                   ! array  for subroutine FFT991
Real, Allocatable :: trigs(:), work(:), xx(:)           ! arrays for subroutine FFT991
Real, Allocatable :: ff(:,:,:), ff2(:,:,:), ff3(:,:,:)  ! arrays using for Fourier transform
Real, Allocatable :: ans(:,:,:)                         ! array of desired correlation
pi       = Atan(1.)*4
xngrid     = 2 * xhalfngrid
mu       = 0.0006
sigma=1

Allocate(trigs(1:3*xhalfngrid+1))
Allocate(work(1:xngrid+1))
Allocate(xx(0:xngrid+1))
Allocate(ff (1:xngrid+2,1:xhalfngrid+1,1:xhalfngrid+1))
Allocate(ff2(1:xngrid+2,1:    xngrid+2,1:xhalfngrid+1))

gamma=(1.*xhalfngrid)/yhalfngrid
Call SET99 (trigs,ifax,xngrid)       ! calculation of the FFT coeff
F_sum(:)=0
Do k=0, xhalfngrid
   Do m=0, yhalfngrid
      Do l=0, verthalfngrid
         p=P_d(k,m,l,q,lambda,mu, gamma, delta)
         b_nml=(3*sigma**2)/(16*p**5)
         F_sum(k)=F_sum(k)+b_nml*weight2(m, l)
      End do
   End Do
End do
    Do j = 1, xhalfngrid+1
      xx(2*j-1) = 0
      xx(2*(j-1))=F_sum(j-1)
    End Do
    Dist(0)=0; l_05=0
    Call  FFT991(xx, work, trigs, ifax, 1, xngrid+2, xngrid, 1, 1)       ! inverse  FFT
    Do j = 0, xhalfngrid
      Dist(j)=mesh_size*j
      If (j>0) Then
         xx(j)= xx(j)/xx(0)
      End if
    End do
    xx(0)=1   
    j=0
  Do j = 1, xhalfngrid
     If (xx(j) .ge. 0.5) Then
        l_05=j
     End if
  End Do
L05_true=((l_05)*(0.5-xx(l_05+1))+(l_05+1)*(xx(l_05)-0.5))/(xx(l_05)-xx(l_05+1))! theoretucal L05 for given lambda
Deallocate(xx)

end subroutine L05_true_definer

subroutine T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu, delta, T05_true)

!*********************************************************************************************************
! PURPOSE
!
! We calculate T05 for given mu and lambda, using a theoretical formula
!
! Authors:  I.Mamay, 2012
!           D Gayfulin, 2013 -2014
! Current code owner: D. Gayfulin
! Last update: March, 2014
!
!**********************************************************************************************************

Use FFT, Only: SET99, FFT991
Implicit None

Real, intent(in)  :: q
Integer, intent(in)  :: xhalfngrid, yhalfngrid 
Integer, intent(in)  :: verthalfngrid
Real   , intent(in)  :: lambda, mu, delta
Real   , intent(out) :: T05_true
Real    :: sigma
Integer :: ngrid
Integer :: k, l, i, j, m            ! loop counters
Real :: pi                          ! number pi, time, normalize time
Real :: P_d                         ! function P_d
Integer :: weight3,step, maxnsteps
Real    :: dt                       ! value of time steps
Real    :: F_sum(-1:800)
Real    :: T_max, p, a_nml, t
Integer :: T05_low, T05_up
Real    :: gamma

pi = 4 * Atan(1.)

sigma=1
F_sum(:)=1;
i=0
step =300
maxnsteps=800

gamma=(1.*xhalfngrid)/yhalfngrid

Do While ((F_sum(i-1) .gt. 0.5) .AND. (i .le. maxnsteps))
!   Print *, mu, i, F_sum(i-1)
   t=i*step
   Do k=0, xhalfngrid
      Do m=0, yhalfngrid
         Do l=0, verthalfngrid
            p=P_d(k,m,l,q,lambda,mu, gamma, delta)
            a_nml=(3*sigma**2)/(16*p**5)*exp(-p*t)*(1+p*t+((p*t)**2)/3.)
            F_sum(i)=F_sum(i)+a_nml*weight3(k, m, l)
         End do
      End Do
  End Do
  If (i .ge. 1) Then
     F_sum(i)=F_sum(i)/F_sum(0)
  End if
  i=i+1
End do
T05_up=i-1
T05_low=T05_up-1
T05_true=step*(T05_up*(F_sum(T05_low)-0.5)+T05_low*(0.5-F_sum(T05_up)))/(F_sum(T05_low)-F_sum(T05_up))


end subroutine T05_true_definer


subroutine Calc_lambda(xngrid, yngrid, vertngrid, q, L05, Size_km, lam, delta, stat)

!=========================================================================================================
! PURPOSE
!
! We calculate spatial correltion parameter lambda
!
! ARGS
! ngrid - horizontal gridzize, 
! vertngrid - vertical gridzize
! q - smoothness
! L05 - spatio scales (km)
! Size_km - Size of the domain (km)
! lam - lambda parameter, we define
!
! Authors:  I.Mamay, 2012
!           D Gayfulin, 2013 -2014
! Current code owner: D. Gayfulin
! Last update: October, 2016
!
!=========================================================================================================

Implicit None

Integer, intent(in)  :: xngrid, yngrid  ! number of a points of the grid
Integer, intent(in)  :: vertngrid
Real, intent(in)  :: q           ! smoothness
Real,    intent(in)  :: L05         ! spatio scales (km)
Real,    intent(in)  :: Size_km     ! Size of domain (km)
Real,    intent(out) :: lam         ! spatio pseudoscales
Real, intent(in)  :: delta
Integer, intent(out) :: stat        ! status
Real               ::  L05arr(1:20)
Real    :: alpha                    ! coefficient  for linear interpolation
Real    :: mesh_size
Real    :: lambda_min, lambda_middle, lambda_max   ! min/max value of lambda
Real    :: lambda_low, lambda_up    ! up/low estimate of lambda
Real    :: L05precision, L05_point, L05_true, L05_true_old, L05_true_low, L05_true_middle, L05_true_up  
Integer :: m                        ! index number for linear interpolation
Integer :: i                        ! loop counter
Integer :: n1, n2, n3               ! ngrid = 100*n1 + 10n2 + n3
Integer :: Lmin, Lmax               ! lambda(L05):[Lmin; Lmax] -> R
Integer :: nmesh, mesh              ! parameters of fragmentation of interval [Lmin; Lmax]
Integer :: mesh_km, L_min, L_max    ! result of functions mesh_km, L_min, L_max
Integer :: ier, ie                  ! return code
Integer :: xhalfngrid, yhalfngrid, verthalfngrid

Real :: Pi  
!=========================
!===== Preparations ======
Pi=4*atan(1.)
stat = 0
xhalfngrid=xngrid*0.5
yhalfngrid=yngrid*0.5
verthalfngrid=vertngrid*0.5
!=========================
L05_point=(L05/Size_km)*xngrid
mesh_size=Size_km/xngrid

!==========================================
!=== Calculating of mesh, Lmin and Lmax ===

!==========================================
L05precision=1./xngrid
lambda_min=1./xngrid
lambda_up=Pi/2
lambda_low=lambda_min

Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_low, delta, L05_true_low)
Do while (L05_true_low .gt. L05_point)
   Print *, 'lower estimate not enough', lambda_low, L05_true_low
   lambda_low=lambda_low/2
   Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_low, delta, L05_true_low)
End Do
Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_up, delta, L05_true_up)
Do while (L05_true_up .lt. L05_point)
   Print *, 'upper estimate not enough', lambda_up, L05_true_up
   lambda_up=lambda_up*2
   Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_low, delta, L05_true_low)
End Do


lambda_middle=(lambda_low+lambda_up)/2
Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_middle, delta, L05_true_middle)
   Do while ((lambda_up-lambda_low) .gt. L05precision)
       If (L05_true_middle .gt. L05_point) Then
          lambda_up=lambda_middle
       Else
          lambda_low=lambda_middle
       End If
       lambda_middle=(lambda_low+lambda_up)/2
       Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_middle, delta, L05_true_middle)
   End do
Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_low, delta, L05_true_low)
Call L05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, mesh_size, q, lambda_up, delta, L05_true_up)

lam=(lambda_up*(L05_point-L05_true_low) + lambda_low*(L05_true_up-L05_point))/(L05_true_up-L05_true_low)
Print "(' ',a, f9.5)", 'lambda=', lam

End subroutine Calc_lambda


subroutine Calc_mu(xngrid, yngrid, vertngrid, q, lambda, T05, mu, delta, stat)

!=========================================================================================================
! PURPOSE
!
! We calculate temporal correlation parameter  mu
!
! ARGS
! xngrid, yngrid - horizontal gridzize,
! vertngrid - vertical gridzize
! q - smoothness
! L05 - spatio scales (km)
! T05 - temporal scales (sec)
! Size_km - Size of the domain (km)
! mu - parameter, we define
!
! Authors:  I.Mamay, 2012
!           D Gayfulin, 2013 -2014
! Current code owner: D. Gayfulin
! Last update: March, 2014
!
!=========================================================================================================

Implicit None

Integer, intent(in)  :: xngrid, yngrid  ! number of grid points
Integer, intent(in)  :: vertngrid
Real, intent(in)  :: q         ! smoothness
Real   , intent(in)  :: lambda    ! lambda parameter
Real   , intent(in)  :: T05       ! temporal scales (sec)
Real   , intent(in)  :: delta
Real   , intent(out) :: mu        ! temporal pseudoscales
Integer, intent(out) :: stat      ! status
Real    :: alpha                  ! coefficient  for linear interpolation
Real    :: mu_min, mu_middle, mu_max   ! min/max value of lambda
Real    :: mu_low, mu_up    ! up/low estimate of lambda
Real    :: T05precision, T05_true, T05_true_up, T05_true_low, T05_true_middle

Integer :: m                      ! index number for linear interpolation
Integer :: i                      ! loop counter
Integer :: n1, n2, n3             ! ngrid = 100*n1 + 10n2 + n3
Integer :: Lmin, Lmax             ! lambda(L05):[Lmin; Lmax] -> R
Integer :: nmesh, mesh            ! parameters of fragmentation of interval [Lmin; Lmax]
Integer :: mesh_km, L_min, L_max  ! result of functions mesh_km, L_min, L_max
Integer :: ier, ie                ! return code
Integer :: xhalfngrid, yhalfngrid, verthalfngrid
!=========================
!===== Preparations ======
stat = 0
verthalfngrid=vertngrid/2
xhalfngrid=xngrid/2
yhalfngrid=yngrid/2
!========================
T05precision=0.000001
mu_low=0.00008
mu_up=0.001
mu_middle=(mu_low+mu_up)/2

!======================================================================
!=========== Calculation of linear interpolation parameters ===========
Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_low, delta, T05_true_low)
Do while (T05_true_low .lt. T05)
   Print *, 'lower estimate not enough', mu_low, T05_true_low
   mu_low=mu_low/2
   Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_low, delta, T05_true_low)
End Do
Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_up, delta, T05_true_up)
Do while (T05_true_up .gt. T05)
   Print *, 'upper estimate not enough', mu_up, T05_true_up
   mu_up=mu_up*2
   Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_up, delta, T05_true_up)
End Do


Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_middle, delta, T05_true_middle)
   Do while ((mu_up-mu_low) .gt. T05precision)
       If (T05_true_middle .gt. T05) Then
          mu_low=mu_middle
       Else
          mu_up=mu_middle
       End If
       mu_middle=(mu_low+mu_up)/2
       Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_middle, delta, T05_true_middle)
   End do
Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_low, delta, T05_true_low)
Call T05_true_definer(xhalfngrid, yhalfngrid, verthalfngrid, q, lambda, mu_up, delta, T05_true_up)

mu=(mu_up*(T05_true_low-T05) + mu_low*(T05-T05_true_up))/(T05_true_low-T05_true_up)
Print "(' ',a, e11.4)", 'mu=', mu
End subroutine Calc_mu

Subroutine Calc_sigma(q, xhalfngrid, yhalfngrid, verthalfngrid, lambda, mu, delta, sigma, stat)

!=========================================================================================================
! PURPOSE
!
! This subroutine calculates the field sigma parameter using variance theoretical formula
!
! ARGS
! q - smoothness
! halfngrid = ngrid / 2
! ngrid - horizontal gridzize,
! vertngrid - vertical gridzize
! lambda, mu - spatio and temporal pseudoscales
! sigma - prameter, we define
!
! Authors: I. Mamay, 2012
!          D. Gayfulin, 2013-2014, 2016
! Current code owner: D.Gayfulin
! Last updated: Oct, 2016

!=========================================================================================================

Implicit None

Real, intent(in)  :: q               ! smoothness parameter
Integer, intent(in)  :: xhalfngrid, yhalfngrid        ! halfngrid = ngrid / 2
Integer, intent(in)  :: verthalfngrid      ! verthalfngrid = vertngrid / 2
Real   , intent(in)  :: lambda, mu, delta  ! spatio and temporal pseudoscales
Real   , intent(out) :: sigma              ! sigma parameter of the generator main equation 
Integer, intent(out) :: stat               ! status
Real    :: pi           ! number pi
Real    :: P_d          ! result of function P_d 
Integer :: k, l, j      ! loop counters
Real    :: B0           ! theoritical field variance
Real    :: gamma
!=========================
!===== Preparations ======
stat  = 0
B0    = 0.0
pi    = 4.0 * Atan(1.0)
gamma=(1.*xhalfngrid)/yhalfngrid
!=========================

!=================================
!======= Calculation of B0 =======
Do j = -xhalfngrid, xhalfngrid
  Do l = -yhalfngrid, yhalfngrid
    Do k = -verthalfngrid, verthalfngrid
      B0 = B0 +  3./(16.*P_d(j,l,k,q,lambda,mu,gamma,delta,stat)**5)
    End Do
  End Do
End Do
If (verthalfngrid .ne. 0) Then
   B0 = B0 / ((2*pi)**6)
Else
   B0 = B0 / ((2*pi)**4)
End If

sigma = 1./ sqrt(B0)
Print "(' ',a, e11.4)", 'sigma=', sigma
!=================================

end subroutine Calc_sigma


end module set_generator_parameters



Real Function P_d(j, l, k, q, lambda, mu, gamma, delta)
!========================================
! This Function calculates the polynomial
! of the Laplacian operator
!========================================
Implicit None
Integer, intent(in) :: j, l, k                          ! indexes of Fourier coefficient
Real, intent(in) :: q                                ! smoothness
Real, intent(in)    :: lambda, mu, gamma, delta          ! spatio and temporal pseudoscales
P_d = mu * ((1. + (j**2+(gamma*l)**2+(delta*k)**2)*lambda**2)**q)       ! polynomyal P(Delta)
End Function P_d

Integer Function weight2(l, m) ! number of different elements of the form (+-l, +-m)
Implicit None
Integer, intent(in) :: l, m
integer t1, t2
t1=2; t2=2
If (l==0) Then
  t1=1
End if
If (m==0) Then
  t2=1
End if
weight2=t1*t2
End Function weight2

Integer Function weight3(l, m, n) ! numbers of different elements of the form (+-l, +-m, +-n)
Implicit None
Integer, intent(in) :: l, m, n
integer t1, t2, t3
t1=2; t2=2; t3=2
If (l==0) Then
  t1=1
End if
If (m==0) Then
  t2=1
End if
If (n==0) Then
  t3=1
End if
weight3=t1*t2*t3
End Function weight3
