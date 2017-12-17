module generator

Use set_generator_parameters, Only: Calc_lambda, Calc_mu, Calc_sigma, Calc_corr

contains

subroutine generator_launch(xngrid, yngrid, vertngrid, q, std, L05, T05, delta, dseed, timestride, nstrides, nsample,&
                            a_dt_min, a_dt_max, mesh_size, n_min_coarse_grid, exp_factor_coarse_grid,&
                            intpl_acceleration, comput_resolution, rand_field, stat)

!=========================================================================================================
! PURPOSE
!
! We calculate the internal parameters lambda, mu and sigma.
! After that we generate Fourier coeffitients = rexi_nml + i * imxi_nml 
! such that Re( xi(-n,-m,-l) ) =  Re( xi(n,m,l) )
! and       Im( xi(-n,-m,-l) ) = -Im( xi(n,m,l) )
! Then we calculate the field  rand_field  using                  
! Fourier tranform of rexi_nml + i * imxi_nml
!
! Arguments
! halfngrid = ngrid / 2
! ngrid - number of a points of the grid, 
! q - smoothness parameter
! std - standard deviation
! L05 - spatio scales (km)
! T05 - temporal scales (sec)
! dseed - seed for random number generator
! timestride - length of time interval for one step
! nstrides - number of time steps
! nsample - number of samples generated
! rand_field - 3D random field
! domain_size - size of the domain (in km)
! a_dt_max - the relation (step of integration)/(time scale).
!
! Authors:  I.Mamay, 2012
!           D Gayfulin, 2013 -2014
! Current code owner: D. Gayfulin
! Last update: October, 2016
!=========================================================================================================

Implicit None

Integer          :: xngrid, yngrid ! number of horizontal points
Integer          :: vertngrid      ! number of vertical points

Real            , intent(in)   :: q                  ! smoothness parameter
Real            , intent(in)   :: std                ! standard deviation
Real            , intent(in)   :: L05, T05           ! spatio and temporal scales (km, sec)
Real            , intent(in)   :: delta              ! vertical scale parameter 
Double Precision, intent(inout):: dseed              ! seed for random number generator
Real            , intent(in)   :: timestride          ! length of time interval
Real            , intent(in)   :: a_dt_min, a_dt_max ! min and max of relation between timestep an timescale in sp. space
Real            , intent(in)   :: mesh_size          ! size of horizontal mesh (in km)
Integer         , intent(in)   :: nstrides              ! number of answers on the orbit
Integer         , intent(in)   :: nsample            ! number of random orbits 
Integer         , intent(in)   :: n_min_coarse_grid  !
Integer         , intent(in)   :: intpl_acceleration ! 1 -use acceleration -1 - otherwise
Integer         , intent(in)   :: comput_resolution  ! we recomend to use values from 5 to 20
Real            , intent(in)   :: exp_factor_coarse_grid ! 
Real,allocatable, intent(out)  :: rand_field(:,:,:,:,:)  ! 3D random field
Integer         , intent(out)  :: stat               ! status

Integer :: gen_dim            ! spatial dimension of the field (equals 2 or 3)
!Integer :: verthalfngrid   ! (number of points of the grid)/2
Integer :: xngrid_ext, xhalfngrid_ext, yngrid_ext, yhalfngrid_ext ! number of points of the extended grid
Integer :: vertngrid_ext, verthalfngrid_ext ! number of points of the extended grid
Integer :: n, m, l            ! 3-dimensional index
Integer :: retcod             ! return code
Real    :: mu, la             ! temporal and spatio pseudoscales
Real    :: sigma              ! standard deviation parameter 
Real    :: xdomain_size, ydomain_size        ! size of the domain (km)
Real, Allocatable :: rexi_nml(:,:,:,:,:)     ! Real  part of fourier coefficient of function rand_field with index n,m,l
Real, Allocatable :: imxi_nml(:,:,:,:,:)     ! Imaginary part of fourier coefficient of function rand_field with index n,m,l
Real,allocatable  :: rand_field_big(:,:,:,:,:),rand_field_full(:,:,:,:,:)   ! extended 3D random field
Character(Len=10) :: ctime0, ctime1, ctime2, ctime3, ctime4  ! CPU time
Character(Len=10) :: ctime_spec1, ctime_spec2
Character(Len=2)  :: c_hh0, c_hh1, c_hh2, c_hh3, c_hh4, c_mm0, c_mm1, c_mm2, c_mm3, c_mm4, char_m
Character(Len=6)  :: c_ss0, c_ss1, c_ss2, c_ss3, c_ss4
Character(Len=20) :: field_name
Real :: time1, time2, time3    ! CPU time in real format
Real :: hh0, hh1, hh2, hh3, hh4, mm0, mm1, mm2, mm3, mm4, ss0, ss1, ss2, ss3, ss4
Real :: define_par_time, integration_time, integration_time_no_interp, fourier_time, final_interp_time
Real :: rand_phase
Real :: intpl_coeff
Real :: percentage
Integer ::  xngrid_coarse, yngrid_coarse, xngrid_big, yngrid_big
Integer,  Allocatable :: ireg_grid(:,:,:)
!=========================
!===== Preparations ======
Call Date_and_Time(time=ctime0)  ! time=hhmmss.sss
stat  = 0
xdomain_size=mesh_size*xngrid
ydomain_size=mesh_size*yngrid
Print *, 'Calculation of the extended grid parameters started'
Call definegrid(xngrid, yngrid, vertngrid, xngrid_ext, yngrid_ext, vertngrid_ext, comput_resolution, L05,& 
                intpl_acceleration, mesh_size,xdomain_size,ydomain_size, xngrid_coarse, yngrid_coarse, intpl_coeff) 
Print "(' ',a, i4)", 'extended xngrid=  ', xngrid_ext
Print "(' ',a, i4)", 'extended yngrid=  ', yngrid_ext
Print "(' ',a, i4)", 'extended vertngrid=  ', vertngrid_ext
If (intpl_acceleration .eq. 1) Then
   xngrid_big=xngrid_ext
   yngrid_big=yngrid_ext
   xngrid_ext=xngrid_coarse
   yngrid_ext=yngrid_coarse
   Print "(' ',a, i4)", 'coarse xngrid=  ', xngrid_coarse
   Print "(' ',a, i4)", 'coarse yngrid=  ', yngrid_coarse
End If
Print "(' ',a, f8.2, a)", 'extended xdomain size=', xdomain_size, ' km'
Print "(' ',a, f8.2, a)", 'extended ydomain size=', ydomain_size, ' km'
xhalfngrid_ext=xngrid_ext/2
yhalfngrid_ext=yngrid_ext/2
verthalfngrid_ext=vertngrid_ext/2
!verthalfngrid = vertngrid/2
! seting generated field dimension (2 or 3)
If (vertngrid .eq. 0) Then
   gen_dim=2
   vertngrid=1
   Print *, 'The generator is launched in 2D mode'
Else
   gen_dim=3
   Print *, 'The generator is launched in 3D mode'
End If
Print *, 'Calculation of internal parameters started'
! =========================Calculation of internal parameters=====================================
!===================================================
!============= Calculation of lambda ===============
Call Calc_lambda(xngrid_ext, yngrid_ext, vertngrid_ext, q, L05, xdomain_size, la, delta, stat)
If (stat .NE. 0) Then
  Print"(/'Subroutine generator. Error in Calc_lambda. stat=',i4)", stat
  Return
End If
!===================================================

!===================================================
!=============== Calculation of mu =================
Call Calc_mu (xngrid_ext, yngrid_ext, vertngrid_ext, q, la, T05, mu, delta, stat)
If (stat .NE. 0) Then
  Print"(/'Subroutine generator. Error in Calc_mu. stat=',i4)", stat
  Return
End If
!===================================================

!===================================================
!========== Calculation of sigma ====== ============
Call Calc_sigma(q, xhalfngrid_ext, yhalfngrid_ext, verthalfngrid_ext, la, mu, delta, sigma, stat)
If (stat .NE. 0) Then
  Print"(/'Subroutine generator. Error in Calc_B0_theoretic. stat=',i4)", stat
  Return
End If
Print *, 'The internal parameters are defined'
!======================================================
!===================== Allocation =====================
Allocate(rand_field(0:xngrid-1,0:yngrid-1,0:vertngrid-1,1:nstrides,1:nsample))
Allocate(rexi_nml(-xhalfngrid_ext:xhalfngrid_ext, -yhalfngrid_ext:yhalfngrid_ext,&
                  -verthalfngrid_ext:verthalfngrid_ext, 1:nstrides, 1:nsample))
Allocate(imxi_nml(-xhalfngrid_ext:xhalfngrid_ext, -yhalfngrid_ext:yhalfngrid_ext,&
                  -verthalfngrid_ext:verthalfngrid_ext, 1:nstrides, 1:nsample))
If (vertngrid .eq. 1) Then
   Allocate(rand_field_big(0:xngrid_ext+1,0:yngrid_ext+1,0:0,1:nstrides,1:nsample))
Else
   Allocate(rand_field_big(0:xngrid_ext+1,0:yngrid_ext+1,0:vertngrid_ext-1,1:nstrides,1:nsample))
End If
Allocate(ireg_grid(0:xhalfngrid_ext, -yhalfngrid_ext:yhalfngrid_ext,&
                  -verthalfngrid_ext:verthalfngrid_ext))
!======================================================

Call Date_and_Time(time=ctime1)  ! time=hhmmss.sss
Print *, 'Integration in spectral space started'; Print *, ' '
!=========================================================================================
!=========== Generating of Fourier coeffitients xi = rexi_nml + i * imxi_nml such that ===
!=========== Re( xi(-n,-m,-l) ) =  Re( xi(n,m,l) ) =======================================
!=========== Im( xi(-n,-m,-l) ) = -Im( xi(n,m,l) ) =======================================
Call generate_coarse_grid(xhalfngrid_ext, yhalfngrid_ext, verthalfngrid_ext,&
                          n_min_coarse_grid, exp_factor_coarse_grid, ireg_grid)
Do l = -verthalfngrid_ext, verthalfngrid_ext       ! loop by set {n.NE.0}
  Do m = -yhalfngrid_ext, yhalfngrid_ext
    Call Date_and_Time(time=ctime_spec1)
    Do n = 1, xhalfngrid_ext
       If (ireg_grid(n, m, 0) .eq. 1) Then
          Call solver(xngrid_ext, yngrid_ext, vertngrid_ext,n,m,l,q,sigma,la,mu,delta,dseed,timestride,nstrides,nsample,&
                     rexi_nml(n,m,l,:,:),imxi_nml(n,m,l,:,:),a_dt_min,a_dt_max,stat)
          If (stat .NE. 0) Then
            Print"(/'Subroutine generator. Error in solver. stat=',i4)", stat
            Return
          End If
          rexi_nml(-n,-m,-l,:,:) =  rexi_nml(n,m,l,:,:)
          imxi_nml(-n,-m,-l,:,:) = -imxi_nml(n,m,l,:,:)
       End If
    End Do
  End Do
End Do
Do l = -verthalfngrid_ext, verthalfngrid_ext       ! loop by set {n=0 .and. m.NE.0}
  Do m = 1, yhalfngrid_ext
    n = 0
    Call solver(xngrid_ext,yngrid_ext,vertngrid_ext,n,m,l,q,sigma,la,mu,delta,dseed,timestride,nstrides,nsample,&
               rexi_nml(n,m,l,:,:),imxi_nml(n,m,l,:,:),a_dt_min,a_dt_max, stat)
      If (stat .NE. 0) Then
      Print"(/'Subroutine generator. Error in solver. stat=',i4)", stat
      Return
    End If
    rexi_nml(-n,-m,-l,:,:) =  rexi_nml(n,m,l,:,:)
    imxi_nml(-n,-m,-l,:,:) = -imxi_nml(n,m,l,:,:)
  End Do
End Do
If (gen_dim .eq. 3) then
Do l = 1, verthalfngrid_ext                ! loop by set {n=0 .and. m=0 .and. l.NE.0}
   m = 0; n = 0
   Call solver(xngrid_ext,yngrid_ext,vertngrid_ext,n,m,l,q, sigma,la,mu,delta,dseed, timestride,nstrides,nsample,&
               rexi_nml(n,m,l,:,:),imxi_nml(n,m,l,:,:),a_dt_min,a_dt_max,stat)
   If (stat .NE. 0) Then
      Print"(/'Subroutine generator. Error in solver. stat=',i4)", stat
      Return
   End If
   rexi_nml(-n,-m,-l,:,:) =  rexi_nml(n,m,l,:,:)
   imxi_nml(-n,-m,-l,:,:) = -imxi_nml(n,m,l,:,:)
  End Do
End if                                   ! call for point {n=0 .and. m=0 .and. l=0}
Call solver(xngrid_ext,yngrid_ext,vertngrid_ext,0,0,0,q,sigma,la,mu,delta,dseed,timestride,nstrides,nsample,&
            rexi_nml(0,0,0,:,:),imxi_nml(0,0,0,:,:),a_dt_min,a_dt_max,stat)
If (stat .NE. 0) Then
  Print"(/'Subroutine generator. Error in solver. stat=',i4)", stat
  Return
End If

!=========================================================================================

Call Date_and_Time(time=ctime_spec2)
c_hh1=ctime1(1:2); c_hh2=ctime_spec2(1:2)
c_mm1=ctime1(3:4); c_mm2=ctime_spec2(3:4)
c_ss1=ctime1(5:10); c_ss2=ctime_spec2(5:10)
Read(c_hh1, *) hh1; Read(c_hh2, *) hh2
Read(c_mm1, *) mm1; Read(c_mm2, *) mm2
Read(c_ss1, *) ss1; Read(c_ss2, *) ss2
Integration_time_no_interp=3600*(hh2-hh1)+60*(mm2-mm1)+ss2-ss1

Print *, 'Integration in spectral space ended'
Print *, 'Interpolation in spectral space started'
Call spectral_interpolation(dseed, xhalfngrid_ext, yhalfngrid_ext, verthalfngrid_ext, ireg_grid, la, mu, delta, sigma, q,&
                                  rexi_nml, imxi_nml, nstrides, nsample)
Call Date_and_Time(time=ctime2)  ! time=hhmmss.sss
Print *, 'Interpolation in spectral space ended'
!==================================================================
!=========== Calculation of function rand_field by ====================
!= Fourier trasnsform of the coefficients rexi_nml + i * imxi_nml =
Print *, ' '
Print *, 'Fourier trasnsform started'
Do n = 1, nsample
  Do m = 1, nstrides
     Percentage=100.*((n-1.)*nstrides+m)/(nstrides*nsample)
     Print "(' ',i3,' % of the Fourier transform done')", Nint(Percentage)
     If (gen_dim .eq. 3) Then
      Call fft3d(xhalfngrid_ext, yhalfngrid_ext, verthalfngrid_ext,&
                 rexi_nml(:,:,:,m,n), imxi_nml(:,:,:,m,n), rand_field_big(:,:,:,m,n), stat)
     Else 
      Call fft2d(xhalfngrid_ext, yhalfngrid_ext, rexi_nml(:,:,:,m,n), imxi_nml(:,:,:,m,n), rand_field_big(:,:,0,m,n), stat)
     End If
    If (stat .NE. 0) Then
      Print"(/'Subroutine generator. Error in fft3d. stat=',i4)", stat
      Return
    End If
  End Do
End Do
Print *, 'Fourier trasnsform ended'
Call Date_and_Time(time=ctime3)
If ((intpl_acceleration .eq. 1) .AND. (intpl_coeff .gt. 1)) Then
   If (vertngrid .eq. 1) Then 
      Allocate(rand_field_full(0:xngrid_big-1, 0:yngrid_big-1, 0:0, 1:nstrides, 1:nsample))
      Call bilinear_interpolation(rand_field_big,rand_field_full,xngrid_ext, yngrid_ext, vertngrid,&
                               xngrid_big, yngrid_big, nstrides, nsample)
   Else
      Allocate(rand_field_full(0:xngrid_big-1, 0:yngrid_big-1, 0:vertngrid_ext-1, 1:nstrides, 1:nsample))
      Call bilinear_interpolation(rand_field_big,rand_field_full,xngrid_ext, yngrid_ext, vertngrid_ext,&
                               xngrid_big, yngrid_big, nstrides, nsample)
   End If
   ! seting field variance and cut
   Do l=0, xngrid-1
      Do m=0, yngrid-1
         Do n=0, vertngrid-1 
            rand_field(l,m,n,:,:)=rand_field_full(l,m,n,:,:)*std
         End Do
      End Do
   End Do
Else
   Do l=0, xngrid-1
      Do m=0, yngrid-1
         Do n=0, vertngrid-1
            rand_field(l,m,n,:,:)=rand_field_big(l,m,n,:,:)*std
         End Do
      End Do
   End Do
End If
Deallocate(rand_field_big)
If ((intpl_acceleration .eq. 1) .AND. ((intpl_coeff .gt. 1))) Deallocate(rand_field_full)
Deallocate(rexi_nml)
Deallocate(imxi_nml)
!=================================================================

Call Date_and_Time(time=ctime4)  ! time=hhmmss.sss
Read(ctime1, *) time1; Read(ctime2, *) time2; Read(ctime3, *) time3
c_hh0=ctime0(1:2); c_hh1=ctime1(1:2); c_hh2=ctime2(1:2);  c_hh3=ctime3(1:2); c_hh4=ctime4(1:2)
c_mm0=ctime0(3:4); c_mm1=ctime1(3:4); c_mm2=ctime2(3:4);  c_mm3=ctime3(3:4); c_mm4=ctime4(3:4)
c_ss0=ctime0(5:10); c_ss1=ctime1(5:10); c_ss2=ctime2(5:10); c_ss3=ctime3(5:10); c_ss4=ctime4(5:10)
Read(c_hh0, *) hh0; Read(c_hh1, *) hh1; Read(c_hh2, *) hh2;  Read(c_hh3, *) hh3;  Read(c_hh4, *) hh4;
Read(c_mm0, *) mm0; Read(c_mm1, *) mm1; Read(c_mm2, *) mm2;  Read(c_mm3, *) mm3;  Read(c_mm4, *) mm4;
Read(c_ss0, *) ss0; Read(c_ss1, *) ss1; Read(c_ss2, *) ss2;  Read(c_ss3, *) ss3;  Read(c_ss4, *) ss4;
define_par_time=3600*(hh1-hh0)+60*(mm1-mm0)+ss1-ss0
If (define_par_time .lt. 0) define_par_time=define_par_time+86400
Integration_time=3600*(hh2-hh1)+60*(mm2-mm1)+ss2-ss1
If (Integration_time .lt. 0) Integration_time=Integration_time+86400
fourier_time=3600*(hh3-hh2)+60*(mm3-mm2)+ss3-ss2
If (fourier_time .lt. 0) fourier_time=fourier_time+86400
final_interp_time=3600*(hh4-hh3)+60*(mm4-mm3)+ss4-ss3
If (final_interp_time .lt. 0) final_interp_time=final_interp_time+86400
Print *, ' '
Print "(' ',a,f8.3,a)", 'Parameters definition CPU time=', define_par_time, ' sec.'
Print "(' ',a,f8.3,a)", 'Integration CPU time without interpolation=', Integration_time_no_interp, ' sec.'
!Print "(' ',a,f8.3,a)", 'Integration CPU time=', Integration_time, ' sec.'
Print "(' ',a,f8.3,a)", 'Interpolation CPU time=', Integration_time-Integration_time_no_interp, ' sec.'
Print "(' ',a,f8.3,a)", 'Fourier transform CPU time=', fourier_time, ' sec.'
Print "(' ',a,f8.3,a)", 'Final interpolation CPU time=', final_interp_time, ' sec.'
end subroutine generator_launch

subroutine solver(xngrid, yngrid, vertngrid, n, m, l, q, sigma, lambda, mu, delta, dseed, timestride, nstrides,&
                  nsample, rexi_nml, imxi_nml,a_dt_min, a_dt_max, stat)

!=========================================================================================================
! PURPOSE
!
! We generate random process xi_nml = rexi_nml + i * imxi_nml
!
! ARGS
! ngrid - number of a points of the grid, 
! n,m,l - 3-dimensional index
! q - smoothness parameter
! sigma - standard deviation 
! lambda, mu - parameters for resp. spatial and temporal correlation
! dseed - seed for random numbers generator
! timestride - length of time step
! nstrides - number of time steps
! nsample - number of samples
! rexi_nml, imxi_nml - real and imaginary part of random process xi_nml                
!
! Authors:  I.Mamay, 2012
!           D Gayfulin, 2013 -2014
! Current code owner: D. Gayfulin
! Last update: March, 2014
!
!=========================================================================================================

Implicit None
Integer         , intent(in)   :: xngrid, yngrid     ! hor. resolution        
Integer         , intent(in)   :: vertngrid          ! vert. resoultion
Integer         , intent(in)   :: n, m, l            ! 3-dimensional index
Real            , intent(in)   :: q                  ! smoothness
Real            , intent(in)   :: sigma              ! standart deviation 
Real            , intent(in)   :: lambda, mu, delta  ! spatio and temporal pseudoscales
Double Precision, intent(inout):: dseed              ! seed for random number generator
Real            , intent(in)   :: timestride          ! length of time interval
Real            , intent(in)   :: a_dt_min, a_dt_max      
Integer         , intent(in)   :: nstrides              ! number of answers on the orbit
Integer         , intent(in)   :: nsample            ! number of random orbits 
Real            , intent(out)  :: rexi_nml(1:nstrides, 1:nsample) ! Real  part of random process xi_nml                  
Real            , intent(out)  :: imxi_nml(1:nstrides, 1:nsample) ! Imaginary part of random process xi_nml                   
Integer         , intent(out)  :: stat               ! status
Integer :: i, j, k           ! loop counters
Integer :: nsubsteps         ! number of steps between two answers
Integer :: nsteps            ! nsteps = nstrides * nsubsteps + 1
Integer :: ier               ! return code
Real    :: delta_t           ! delta_t = timestride / nsubsteps 
Real    :: a_nml             ! a_nml = P_dn(n,m,l,...)
Real    :: b_nml             ! Var(xi_mnl)
Real    :: P_dn              ! value of the function P_dn
Real    :: discr_corr        ! we multiply our field by this factor because the variance in continious case
                             ! is different from variance of the field, obtained from our discrete integration scheme
Real, Allocatable :: alpha(:)     ! array of random uncorrelated normally distributed values
Real, Allocatable :: rexi(:,:)    ! Real  part of random process xi
Real, Allocatable :: imxi(:,:)    ! Imaginary part of random process xi
Real :: xi_first(1:3*nsample)     ! Since our differential equation has order three, we should define 
                                  ! three initial steps to start the integration
Real :: kappa                     ! parameter of the inegrated equation (1+a_nml*dt)
Real :: tau                       ! 1/a_nml - temporal scale of integrated process
Real :: percentage                ! percentage of integration done
Real :: c1, c2             ! correlation for 1 step and for 2 steps of integration
!
Real :: a_dt_curr
Real :: max_wavenum, curr_wavenum, start_var_wavenum
Real :: gamma

!=========================
!===== Preparations ======
stat = 0
gamma=(1.*xngrid)/yngrid
max_wavenum=Sqrt((xngrid/2)**2+(yngrid/2)**2+(delta*vertngrid/2)**2)
curr_wavenum=Sqrt(n**2+(gamma*m)**2+(delta*l)**2)
a_dt_curr=a_dt_min+((curr_wavenum)/(max_wavenum))**2*(a_dt_max-a_dt_min)
!=========================
If ((n == xngrid*0.5) .AND. (m==0) .AND. (vertngrid .ne. 0)) Then
    Percentage=100*(l+vertngrid*0.5)/(1.*vertngrid)
    Print "(' ',i3,' % of the time integration done')", Nint(Percentage)
End if
!Print *, n, m
!============================================
!=========== Calculation of a_nml ===========
a_nml = P_dn(n, m, l, q, lambda, mu, gamma, delta, stat)
If (stat .NE. 0) Then
  Print"(/'solver. Error in Function P_dn. stat=',i4)", stat
  Return
End If
!============================================

!====================================================
!=== Calculation of nsubsteps, delta_t and nsteps ===
If (timestride/4 < a_dt_curr / a_nml) Then
  nsubsteps = 4
Else
  nsubsteps = 1 + NINT(timestride / (a_dt_curr / a_nml))
End If
tau=1./a_nml
delta_t = timestride / (nsubsteps * 1.) 
nsteps  = nstrides * nsubsteps + 1
!====================================================

!================================================================
!========================= Allocations ==========================
Allocate(alpha(1:nsteps), Stat=ier)
If(ier.ne.0) Then
  Print "(/'solver error. Allocate alpha failed. ier=', i5)", ier 
  stat=ier ;  Return
End If
Allocate(rexi(1:nsteps,1:nsample))
If(ier.ne.0) Then
  Print "(/'solver error. Allocate rexi failed. ier=', i5)", ier 
  stat=ier ;  Return
End If
Allocate(imxi(1:nsteps,1:nsample))
If(ier.ne.0) Then
  Print "(/'solver error. Allocate imxi failed. ier=', i5)", ier 
  stat=ier ;  Return
End If
!================================================================
discr_corr=sqrt(1. + 0.5 * a_nml * delta_t)**5/Sqrt((a_nml * delta_t)**4*1./6.+&
(a_nml * delta_t)**3*2./3.+(a_nml * delta_t)**2*5./3.+(a_nml * delta_t)*2.+1)

!==============================================================================
!======= Generating of random field rexi_0 and calculation of rexi(1,:) ======= 
kappa=1+a_nml*delta_t
Call GGNML(dseed, 3*nsample, xi_first) 
!Call Calc_corr(a_nml, delta_t, c1)
!Call Calc_corr(a_nml, 2*delta_t, c2)
c1=(3*kappa*(kappa**2+1))/(kappa**4+4*kappa**2+1)
c2=(6*kappa**2)/(kappa**4+4*kappa**2+1)
Do k=1, nsample
   rexi(1,k)=xi_first((k-1)*3+1)
   rexi(2,k)=c1*xi_first((k-1)*3+1)+Sqrt(1-c1**2)*xi_first((k-1)*3+2)
   rexi(3,k)=c2*xi_first((k-1)*3+1)+((c1*(1-c2))/Sqrt(1-c1**2))*xi_first((k-1)*3+2)+&
             Sqrt(1-c2**2-((c1*(1-c2))/Sqrt(1-c1**2))**2)*xi_first((k-1)*3+3)
End do
If ((n .ne. 0) .OR. (m .ne. 0) .OR. (l .ne. 0)) Then
   rexi(1:3,:)=rexi(1:3,:)*sigma*Sqrt(3./(16.*a_nml**5))* (1. / sqrt(2.)) ! seting variance of imxi
Else
  rexi(1:3,:)=rexi(1:3,:) * ( sigma / sqrt((16./3.)*a_nml**5) )! seting variance of rexi
End If

!==============================================================================

!===================================================================
!===== Generating of white noise alpha and calculation of rexi ===== 
Do j = 1, nsample
  Call GGNML(dseed, nsteps, alpha(:))
  If(n == 0 .and. m == 0 .and. l == 0) Then
    alpha = alpha * sqrt(delta_t)
  Else
    alpha = alpha * sqrt(delta_t) / Sqrt(2.)  
  End If
  Do i = 4, nsteps
    ! integraton step
    rexi(i,j) = (3*kappa**2*rexi(i-1,j)-3*kappa*rexi(i-2,j)+rexi(i-3,j)+sigma*delta_t**2*alpha(i-1))/(kappa**3)
  End Do 
End Do
rexi(4:nsteps,:)=rexi(4:nsteps,:)*discr_corr       

!===================================================================

!==============================================================================
!======= Generating of random field imxi_0 and calculation of imxi(1,:) ======= 
If(n == 0 .and. m == 0 .and. l == 0) Then
  imxi(:,:) = 0
Else  
  Call GGNML(dseed, 3*nsample, xi_first) 
Do k=1, nsample
   imxi(1,k)=xi_first((k-1)*3+1)
   imxi(2,k)=c1*xi_first((k-1)*3+1)+Sqrt(1-c1**2)*xi_first((k-1)*3+2)
   imxi(3,k)=c2*xi_first((k-1)*3+1)+((c1*(1-c2))/Sqrt(1-c1**2))*xi_first((k-1)*3+2)+&
             Sqrt(1-c2**2-((c1*(1-c2))/Sqrt(1-c1**2))**2)*xi_first((k-1)*3+3)
End Do
   imxi(1:3,:)=imxi(1:3,:)*sigma*Sqrt(3./(16.*a_nml**5))* (1. / sqrt(2.)) ! seting variance of imxi
End If
!==============================================================================

!===================================================================
!===== Generating of white noise alpha and calculation of imxi ===== 
If(n == 0 .and. m == 0 .and. l == 0) Then
  imxi(:,:) = 0
Else  
  Do j = 1, nsample
    Call GGNML(dseed, nsteps, alpha(:))
    alpha = alpha * sqrt(delta_t) / Sqrt(2.)
    ! integraton step
    Do i = 4, nsteps
       imxi(i,j) = (3*kappa**2*imxi(i-1,j)-3*kappa*imxi(i-2,j)+imxi(i-3,j)+sigma*delta_t**2*alpha(i-1))/(kappa**3)
    End Do 
  End Do  
End If
imxi(4:nsteps,:)=imxi(4:nsteps,:)*discr_corr
!==========================================
!=========== Writing the answer ===========
Do j = 1, nsample  
  Do i = 1, nstrides 
    rexi_nml(i, j) = rexi(i * nsubsteps + 1, j)
    imxi_nml(i, j) = imxi(i * nsubsteps + 1, j)
  End Do
End Do
!=================
!= Deallocations =
Deallocate(rexi)
Deallocate(imxi)
Deallocate(alpha)
!=================

end subroutine solver

subroutine fft2d(xhalfngrid, yhalfngrid, rexi_nml, imxi_nml, rand_field, stat)

!=========================================================================================================
! PURPOSE
!
! This subroutine makes the inverse 2-dimensional fast fourier transformation
!
! Foe detailed description see FFT manual
!
! ARGS
! halfngrid = ngrid / 2
! ngrid - number of a points of the grid,
! rexi_nml, imxi_nml - real and imaginary part of Fourier coeffitients
! rand_field = Inverse_FFT( xi_nml )
!
! Author: D. Gayfulin, 2014
! Current code owner: D.Gayfulin
! Last updated: Jan, 2014
!=========================================================================================================

Use FFT, Only: SET99, FFT991
Implicit None

Integer, intent(in) :: xhalfngrid, yhalfngrid ! halfngrid = ngrid / 2
Real   , intent(in) :: rexi_nml(-xhalfngrid:xhalfngrid, -yhalfngrid:yhalfngrid)
                                                       ! Real  part of coeffitients
Real   , intent(in) :: imxi_nml(-xhalfngrid:xhalfngrid, -yhalfngrid:yhalfngrid)
                                                     ! Imaginary part of coeffitients
Real   , intent(out):: rand_field  (0:2*xhalfngrid + 1, 0:2*yhalfngrid + 1)       ! 2D random field
Integer, intent(out):: stat       ! status
Real                :: pi         ! number pi
Integer             :: xngrid, yngrid, vertngrid       ! grid resolution
Integer             :: i, j, k, m,p ,l    ! 3-dimensional indexes
Integer             :: x_ifax(1:13)                ! array  for subroutine FFT991
Integer             :: y_ifax(1:13)                ! array  for subroutine FFT991
Real, Allocatable   :: x_trigs(:), x_work(:), x_xx(:)  ! arrays for subroutine FFT991
Real, Allocatable   :: y_trigs(:), y_work(:), y_xx(:)  ! arrays for subroutine FFT991
Real, Allocatable   :: Reel(:,:), Imel(:,:)  ! arrays using for Fourier transform
Real, Allocatable   :: phi3(:,:), psi3(:,:)  ! arrays using for Fourier transform

!=========================
!===== Preparations ======
stat = 0
xngrid = 2 * xhalfngrid
yngrid = 2 * yhalfngrid
pi   = 4.0 * Atan(1.)
!============================================================
!======================= Allocations ========================
Allocate(x_trigs(1:3*xhalfngrid+1)); Allocate(y_trigs(1:3*yhalfngrid+1))
Allocate(x_work(1:xngrid+1)); Allocate(y_work(1:yngrid+1))
Allocate(x_xx(0:xngrid+1)); Allocate(y_xx(0:yngrid+1))
Allocate(Reel(-xhalfngrid:xhalfngrid,0:yngrid-1))
Allocate(Imel(-xhalfngrid:xhalfngrid,0:yngrid-1))

!============================================================
Call SET99 (x_trigs,x_ifax,xngrid)       ! calculation of the FFT coeff
Call SET99 (y_trigs,y_ifax,yngrid)       ! calculation of the FFT coeff
!=======================First stage==========================
Do l = -xhalfngrid, xhalfngrid

    Do m = 0, yhalfngrid
      y_xx(2*m) = 0.5 * (rexi_nml(l,m) + rexi_nml(-l,m))
    End Do

    Do m = 0, yhalfngrid-1
      y_xx(2*m+1) = 0.5 * (imxi_nml(l,m) + imxi_nml(-l,m))
    End Do

    y_xx(1)    = 0
    y_xx(yngrid)=2*y_xx(yngrid) ! this coeffitent should be multiplied by 2, see FFT manual
    y_xx(yngrid+1)=0
    Call  FFT991(y_xx, y_work, y_trigs, y_ifax, 1, yngrid+2, yngrid, 1, 1)       ! inverse  FFT
    Reel(l,0:yngrid-1) = y_xx(0:yngrid-1)

    Do m = 0, yhalfngrid
       y_xx(2*m) = 0.5 * (imxi_nml(l,m) - imxi_nml(-l,m))
    End Do

    Do m = 0, yhalfngrid-1
      y_xx(2*m+1) =  0.5 * (rexi_nml(-l,m) - rexi_nml(l,m))
    End Do
    y_xx(yngrid)=2*y_xx(yngrid) ! this coeffitent should be multiplied by 2, see FFT manual
    y_xx(1)      = 0
    y_xx(yngrid+1) = 0

    Call  FFT991(y_xx, y_work, y_trigs, y_ifax, 1, yngrid+2, yngrid, 1, 1)       ! inverse  FFT
    Imel(l,0:yngrid-1) = y_xx(0:yngrid-1)

End Do
Allocate( phi3(0:xngrid-1,0:yngrid-1))
Allocate( psi3(0:xngrid-1,0:yngrid-1))
!======================================Second stage=========================
Do j = 0, yngrid-1
    Do l = 0, xhalfngrid
       x_xx(2*l) = Reel(l,j)
    End Do

    Do l = 0, xhalfngrid-1
        x_xx(2*l+1)=Imel(l,j)
    End Do
    x_xx(xngrid)=2*x_xx(xngrid)    ! this coeffitent should be multiplied by 2, see FFT manual
    x_xx(1)      = 0
    x_xx(xngrid+1) = 0

    Call  FFT991(x_xx, x_work, x_trigs, x_ifax, 1, xngrid+2, xngrid, 1, 1)       ! inverse  FFT
    phi3(0:xngrid-1,j) = x_xx(0:xngrid-1) / (4*pi**2)

    Do l = 0, xhalfngrid
      x_xx(2*l) = 0.5 * (Imel(l,j) - Imel(-l,j))
    End Do

    Do l = 0, xhalfngrid-1
      x_xx(2*l+1) = 0.5 * (Reel(l,j) - Reel(-l,j))
    End Do
    x_xx(xngrid)=2*x_xx(xngrid)    ! this coeffitent should be multiplied by 2, see FFT manual
    x_xx(1)      = 0
    x_xx(xngrid+1) = 0

    Call  FFT991(x_xx, x_work, x_trigs, x_ifax, 1, xngrid+2, xngrid, 1, 1)       ! inverse  FFT
    psi3(0:xngrid-1,j) = x_xx(0:xngrid-1) / (4*pi**2)

End Do
!======================Deallocation=================================
Deallocate(Reel)
Deallocate(Imel)
Deallocate(x_trigs, y_trigs)
Deallocate(x_work, y_work)
Deallocate(x_xx, y_xx)
!======= Output of the result =======
Do i = 0, xngrid - 1
  Do j = 0, yngrid - 1
      rand_field(i,j) = phi3(i,j)
  End Do
End Do
!=====================================

Deallocate(phi3)
Deallocate(psi3)

End subroutine fft2d

subroutine fft3d(xhalfngrid, yhalfngrid, verthalfngrid, rexi_nml, imxi_nml, rand_field, stat)

!=========================================================================================================
! PURPOSE
!
! This subroutine makes the inverse 3-dimensional fast fourier transformation
!
! Foe detailed description see FFT manual
!
! ARGS
! halfngrid = ngrid / 2, where ngrid is the horizontal gridsize,
! verthalfngrid= vertngrid/ 2, where vertngrid is the vertical gridsize
! rexi_nml, imxi_nml - real and imaginary part of Fourier coeffitients
! rand_field = Inverse_FFT( xi_nml )
!
! Authors: I. Mamay, 2012
!          D. Gayfulin, 2013-2014
! Current code owner: D.Gayfulin
! Last updated: Jan, 2014
!=========================================================================================================


Use FFT, Only: SET99, FFT991
Implicit None

Integer, intent(in) :: xhalfngrid, yhalfngrid   ! halfngrid = ngrid / 2     
Integer, intent(in) :: verthalfngrid
Real   , intent(in) :: rexi_nml(-xhalfngrid:xhalfngrid, -yhalfngrid:yhalfngrid, -verthalfngrid:verthalfngrid)
                                                       ! Real  part of random process xi
Real   , intent(in) :: imxi_nml(-xhalfngrid:xhalfngrid, -yhalfngrid:yhalfngrid, -verthalfngrid:verthalfngrid) 
                                                     ! Imaginary part of random process xi
Real   , intent(out):: rand_field  (0:2*xhalfngrid + 1, 0:2*yhalfngrid + 1, 0:2*verthalfngrid + 1)       ! 3D random field
Integer, intent(out):: stat       ! status
Real                :: pi         ! number pi
Integer             :: xngrid, yngrid, vertngrid       ! number of a points of the grid
Integer             :: i, j, k, m,p ,l    ! 3-dimensional index
Integer             :: x_ifax(1:13)                ! array  for subroutine FFT991
Integer             :: y_ifax(1:13)                ! array  for subroutine FFT991
Integer             :: v_ifax(1:13)                ! array  for subroutine FFT991
Real, Allocatable   :: x_trigs(:), x_work(:), x_xx(:)  ! arrays for subroutine FFT991
Real, Allocatable   :: y_trigs(:), y_work(:), y_xx(:)  ! arrays for subroutine FFT991
Real, Allocatable   :: v_trigs(:), v_work(:), v_xx(:)  ! arrays for subroutine FFT991
Real, Allocatable   :: Redlm (:,:,:), Imdlm (:,:,:)  ! arrays using for Fourier transform
Real, Allocatable   :: Reel(:,:,:), Imel(:,:,:)  ! arrays using for Fourier transform
Real, Allocatable   :: phi3(:,:,:), psi3(:,:,:)  ! arrays using for Fourier transform

!=========================
!===== Preparations ======
stat = 0
xngrid = 2 * xhalfngrid
yngrid = 2 * yhalfngrid
vertngrid=2*verthalfngrid
pi   = 4.0 * Atan(1.)
!=========================

!============================================================
!======================= Allocations ========================
Allocate(x_trigs(1:3*xhalfngrid+1)); Allocate(y_trigs(1:3*yhalfngrid+1)); Allocate(v_trigs(1:3*verthalfngrid+1))
Allocate(x_work(1:xngrid+1)); Allocate(y_work(1:yngrid+1)); Allocate(v_work(1:vertngrid+1))
Allocate(x_xx(0:xngrid+1)); Allocate(y_xx(0:yngrid+1)); Allocate(v_xx(0:vertngrid+1))
Allocate(Redlm(-xhalfngrid:xhalfngrid,-yhalfngrid:yhalfngrid, 0:vertngrid-1))
Allocate(Imdlm(-xhalfngrid:xhalfngrid,-yhalfngrid:yhalfngrid, 0:vertngrid-1))
Allocate(Reel(-xhalfngrid:xhalfngrid,0:yngrid-1,0:vertngrid-1))
Allocate(Imel(-xhalfngrid:xhalfngrid,0:yngrid-1,0:vertngrid-1))
!============================================================

!================================================================
!=============== First stage=====================================
Call SET99 (x_trigs,x_ifax,xngrid)       ! calculation of the FFT coeff
Call SET99 (y_trigs,y_ifax,yngrid)       ! calculation of the FFT coeff
Call SET99 (v_trigs,v_ifax, vertngrid)
Do l = -xhalfngrid, xhalfngrid
  Do m = -yhalfngrid, yhalfngrid

    Do p = 0, verthalfngrid
      v_xx(2*p) = 0.5 * (rexi_nml(l,m,p) + rexi_nml(-l,-m,p))
    End Do

    Do p = 0, verthalfngrid-1
      v_xx(2*p+1) = 0.5 * (imxi_nml(l,m,p) + imxi_nml(-l,-m,p))
    End Do

    v_xx(1)    = 0
    v_xx(vertngrid)=2*v_xx(vertngrid)
    v_xx(vertngrid+1)=0
    Call  FFT991(v_xx, v_work, v_trigs, v_ifax, 1, vertngrid+2, vertngrid, 1, 1)       ! inverse  FFT
    Redlm(l,m,0:vertngrid-1) = v_xx(0:vertngrid-1)
   
    Do p = 0, verthalfngrid
       v_xx(2*p) = 0.5 * (imxi_nml(l,m,p) - imxi_nml(-l,-m,p))
    End Do

    Do p = 0, verthalfngrid-1
      v_xx(2*p+1) =  0.5 * (rexi_nml(-l,-m,p) - rexi_nml(l,m,p))
    End Do
    v_xx(vertngrid)=2*v_xx(vertngrid)
    v_xx(1)      = 0
    v_xx(vertngrid+1) = 0

    Call  FFT991(v_xx, v_work, v_trigs, v_ifax, 1, vertngrid+2, vertngrid, 1, 1)       ! inverse  FFT
    Imdlm(l,m,0:vertngrid-1) = v_xx(0:vertngrid-1)
    End Do
End Do
Deallocate(v_xx); Deallocate(v_trigs); Deallocate(v_work)
!================================================================

!================================================================
!=============== Second stage ===================================
Do l = -xhalfngrid, xhalfngrid
  Do k = 0, vertngrid-1

    Do m = 0, yhalfngrid
      y_xx(2*m) = 0.5 * ( Redlm(l,m,k) + Redlm(-l,m,k))
    End Do

    Do m = 0, yhalfngrid-1
      y_xx(2*m+1) = 0.5 * (Imdlm(l,m,k) + Imdlm(-l,m,k))
    End Do

    y_xx(1)      = 0
    y_xx(yngrid+1) = 0
    y_xx(yngrid)=2*y_xx(yngrid)
    Call  FFT991(y_xx, y_work, y_trigs, y_ifax, 1, yngrid+2, yngrid, 1, 1)       ! inverse  FFT
    Reel(l,0:yngrid-1,k) = y_xx(0:yngrid-1)

    Do m = 0, yhalfngrid
      y_xx(2*m) = 0.5 * (Imdlm(l,m,k) - Imdlm(-l,m,k))
    End Do

    Do m = 0, yhalfngrid-1
       y_xx(2*m+1) = 0.5 * (Redlm(-l,m,k) - Redlm(l,m,k))
    End Do

    y_xx(1)      = 0
    y_xx(yngrid+1) = 0
    y_xx(yngrid)=2*y_xx(yngrid)
    Call  FFT991(y_xx, y_work, y_trigs, y_ifax, 1, yngrid+2, yngrid, 1, 1)       ! inverse  FFT
    Imel(l,0:yngrid-1,k) = y_xx(0:yngrid-1)

  End Do
End Do
!================================================================

Deallocate(Redlm)
Deallocate(Imdlm)
Allocate(phi3(0:xngrid-1,0:yngrid-1,0:vertngrid-1))
Allocate(psi3(0:xngrid-1,0:yngrid-1,0:vertngrid-1))

!================================================================
!=============== Third stage ===================================
Do j = 0, yngrid-1
  Do k = 0, vertngrid-1

    Do l = 0, xhalfngrid
       x_xx(2*l) =  Reel(l,j,k)
    End Do

    Do l = 0, xhalfngrid-1
        x_xx(2*l+1)=Imel(l,j,k)
    End Do
    x_xx(xngrid)=2*x_xx(xngrid)
    x_xx(1)      = 0
    x_xx(xngrid+1) = 0

    Call  FFT991(x_xx, x_work, x_trigs, x_ifax, 1, xngrid+2, xngrid, 1, 1)       ! inverse  FFT
    phi3(0:xngrid-1,j,k) = x_xx(0:xngrid-1) / (8*pi**3)

    Do l = 0, xhalfngrid
      x_xx(2*l) = 0.5 * (Imel(l,j,k) - Imel(-l,j,k))
    End Do

    Do l = 0, xhalfngrid-1
      x_xx(2*l+1) = 0.5 * (Reel(l,j,k) - Reel(-l,j,k))
    End Do
    x_xx(xngrid)=2*x_xx(xngrid)
    x_xx(1)      = 0
    x_xx(xngrid+1) = 0

    Call  FFT991(x_xx, x_work, x_trigs, x_ifax, 1, xngrid+2, xngrid, 1, 1)       ! inverse  FFT
    psi3(0:xngrid-1,j,k) = x_xx(0:xngrid-1) / (8*pi**3) 

  End Do
End Do
!================================================================

Deallocate(Reel)
Deallocate(Imel)
Deallocate(x_trigs, y_trigs)
Deallocate(x_work, y_work)
Deallocate(x_xx, y_xx)

!=====================================
!======= Calculation of rand_field =======
Do i = 0, xngrid - 1
  Do j = 0, yngrid - 1
    Do k = 0, vertngrid - 1
      rand_field(i,j,k) = phi3(i,j,k)
    End Do
  End Do
End Do
!=====================================

Deallocate(phi3)
Deallocate(psi3)

end subroutine fft3d 

Subroutine generate_coarse_grid(xhalfngrid, yhalfngrid, verthalfngrid, n_min_coarse_grid, exp_factor_coarse_grid, ireg_grid)
!=========================================================================================================
! PURPOSE
!
! We perform integration in spectral space not for all wavevectors, but for coarser grid and then interpolate
! Author: D. Gayfulin, 2016
! Current code owner: D.Gayfulin
! Last updated: Feb, 2016
!
!==========================================================
Implicit None
!
Integer, intent(in)  :: xhalfngrid, yhalfngrid, verthalfngrid
Integer, intent(in)  :: n_min_coarse_grid
Real, Intent(in)     :: exp_factor_coarse_grid
Integer, intent(Out) :: ireg_grid(0:xhalfngrid,-yhalfngrid:yhalfngrid,-verthalfngrid:verthalfngrid)
!
Integer :: i, j, m, n, l
Integer :: coarse_points(0:xhalfngrid+yhalfngrid)
!===================Execituon========================
ireg_grid=0
Do i=0, n_min_coarse_grid
   coarse_points(i)=i
End Do
n=n_min_coarse_grid; i=n_min_coarse_grid
m=n; j=i
Do While (n .lt. xhalfngrid)
   n=Nint(exp_factor_coarse_grid*n); i=i+1
   If (n .gt. xhalfngrid) n=xhalfngrid
   coarse_points(i)=n
End do
Do While (m .lt. yhalfngrid)
   m=Nint(exp_factor_coarse_grid*m); j=j+1
   If (m .gt. yhalfngrid) m=yhalfngrid
   coarse_points(j)=m
End do
Do n=0, i
   Do m= 0, j
      ireg_grid(coarse_points(n), coarse_points(m), :)=1
      ireg_grid(coarse_points(n), -coarse_points(m), :)=1
   End Do
End Do
End Subroutine generate_coarse_grid

Subroutine spectral_interpolation(dseed, xhalfngrid, yhalfngrid, verthalfngrid, ireg_grid, lambda, mu, sigma, delta, q,&
                                  rexi_nml, imxi_nml, nstrides, nsample)
!=========================================================================================================
! PURPOSE
!
! Interpolation of random field in spectral space
! Author: D. Gayfulin, 2016
! Current code owner: D.Gayfulin
! Last updated: Feb, 2016
!

!==========================================================
Implicit None
!
Double Precision, intent(inout):: dseed

Integer, intent(in) :: xhalfngrid, yhalfngrid, verthalfngrid, nstrides, nsample
Integer, intent(in) :: ireg_grid(0:xhalfngrid,-yhalfngrid:yhalfngrid,-verthalfngrid:verthalfngrid)
Real, intent(in)    :: lambda, mu, sigma, delta, q
Real, intent(inout) :: rexi_nml(-xhalfngrid:xhalfngrid,-yhalfngrid:yhalfngrid,&
                                   -verthalfngrid:verthalfngrid, 1:nstrides, 1:nsample)
Real, intent(inout) :: imxi_nml(-xhalfngrid:xhalfngrid,-yhalfngrid:yhalfngrid,&
                                   -verthalfngrid:verthalfngrid, 1:nstrides, 1:nsample)
Integer :: i, j, m, n, l, curr_answ, curr_sample, stat
Integer :: x1, x2, y1, y2, x,y
Real :: P_dn
Real :: fx1y1, fx1y2, fx2y2, fx2y1
Real :: varnml, varx1y1, varx1y2, varx2y2, varx2y1
Real :: a_nml, a_nmlx1y1, a_nmlx1y2, a_nmlx2y2, a_nmlx2y1
Real :: wx1y1, wx1y2, wx2y1, wx2y2
Real :: rand_phase
Real :: gamma
Real :: Pi
!===================Execituon========================
Pi=4*atan(1.)
gamma=(1.*xhalfngrid)/yhalfngrid
Do m = -yhalfngrid, yhalfngrid
   Do n = 1, xhalfngrid
      If (ireg_grid(n, m, 0) .eq. 0) Then
         x1=n
         Do While (ireg_grid(x1, yhalfngrid, 0) .eq. 0)
            x1=x1-1
         End Do
         x2=n
         Do While (ireg_grid(x2, yhalfngrid, 0) .eq. 0)
            x2=x2+1
         End Do
         y2=m
         Do While (ireg_grid(xhalfngrid, Nint(Abs(1.*y2)), 0) .eq. 0)
            y2=y2+1
         End Do
         y1=m
         Do While (ireg_grid(xhalfngrid, Nint(Abs(1.*y1)), 0) .eq. 0)
            y1=y1-1
         End Do
         Do l=-verthalfngrid, verthalfngrid
            a_nml=P_dn(n, m, l, q, lambda, mu, gamma, delta, stat)
            varnml=sigma*Sqrt(3./(16.*a_nml**5))*(1./sqrt(2.))

            a_nmlx1y1=P_dn(x1, y1, l, q, lambda, mu, gamma, delta, stat)
            varx1y1=sigma*Sqrt(3./(16.*a_nmlx1y1**5))*(1./sqrt(2.))

            a_nmlx1y2=P_dn(x1, y2, l, q, lambda, mu, gamma, delta, stat)
            varx1y2=sigma*Sqrt(3./(16.*a_nmlx1y2**5))*(1./sqrt(2.))

            a_nmlx2y2=P_dn(x2, y2, l, q, lambda, mu, gamma, delta, stat)
            varx2y2=sigma*Sqrt(3./(16.*a_nmlx2y2**5))*(1./sqrt(2.))

            a_nmlx2y1=P_dn(x2, y1, l, q, lambda, mu, gamma, delta, stat)
            varx2y1=sigma*Sqrt(3./(16.*a_nmlx2y1**5))*(1./sqrt(2.))
            Do curr_answ=1, nstrides
               Do curr_sample=1, nsample
                  fx1y1=rexi_nml(x1,y1,l,curr_answ,curr_sample)
                  fx1y2=rexi_nml(x1,y2,l,curr_answ,curr_sample)
                  fx2y2=rexi_nml(x2,y2,l,curr_answ,curr_sample)
                  fx2y1=rexi_nml(x2,y1,l,curr_answ,curr_sample)
                  x=(n-x1)/(1.*x2-x1)
                  y=(m-y1)/(1.*y2-y1)
                  wx1y1=1-x-y+x*y
                  wx1y2=x-x*y
                  wx2y1=y-x*y
                  wx2y2=x*y
                  rexi_nml(n,m,l,curr_answ,curr_sample)=wx1y1*fx1y1+wx1y2*fx1y2+wx2y1*fx2y1+wx2y2*fx2y2
                  rexi_nml(n,m,l,curr_answ,curr_sample)=rexi_nml(n,m,l,curr_answ,curr_sample)*Sqrt(varnml)/&
                                Sqrt(wx1y1**2*varx1y1+wx1y2**2*varx1y2+wx2y1**2*varx2y1+wx2y2**2*varx2y2)

                  fx1y1=imxi_nml(x1,y1,l,curr_answ,curr_sample)
                  fx1y2=imxi_nml(x1,y2,l,curr_answ,curr_sample)
                  fx2y2=imxi_nml(x2,y2,l,curr_answ,curr_sample)
                  fx2y1=imxi_nml(x2,y1,l,curr_answ,curr_sample)
                  imxi_nml(n,m,l,curr_answ,curr_sample)=wx1y1*fx1y1+wx1y2*fx1y2+wx2y1*fx2y1+wx2y2*fx2y2
                  imxi_nml(n,m,l,curr_answ,curr_sample)=imxi_nml(n,m,l,curr_answ,curr_sample)*Sqrt(varnml)/&
                                Sqrt(wx1y1**2*varx1y1+wx1y2**2*varx1y2+wx2y1**2*varx2y1+wx2y2**2*varx2y2)
              End Do
            End Do
            Call GGUBS(dseed, 1, rand_phase)
            rand_phase=2*Pi*rand_phase
            rexi_nml(n,m,l,:,:)=rexi_nml(n,m,l,:,:)*Cos(rand_phase)-imxi_nml(n,m,l,:,:)*Sin(rand_phase)
            imxi_nml(n,m,l,:,:)=rexi_nml(n,m,l,:,:)*Sin(rand_phase)+imxi_nml(n,m,l,:,:)*Cos(rand_phase)
            rexi_nml(-n,-m,-l,:,:) =  rexi_nml(n,m,l,:,:)
            imxi_nml(-n,-m,-l,:,:) = -imxi_nml(n,m,l,:,:)
         End Do
      End If
   End Do
End Do

End Subroutine spectral_interpolation

Subroutine definegrid(xngrid_old, yngrid_old, vertngrid_old, xngrid_new, yngrid_new, vertngrid_new, comput_resolution,&
           L05, intpl_acceleration,  mesh_size, xdomain_size, ydomain_size, xngrid_coarse, yngrid_coarse, intpl_coeff)
!=========================================================================================================
! PURPOSE
!
! This subroutine calculates the prameters of extended area. 
! We expand our domain size to prevent correlations of the borders
!
! Author: D. Gayfulin, 2016
! Current code owner: D.Gayfulin
! Last updated: Oct, 2016
!
!==========================================================
Implicit None
!
Integer, intent(in)    :: xngrid_old, yngrid_old, vertngrid_old
Integer, intent(in)    :: intpl_acceleration
Integer, intent(in)    :: comput_resolution
Integer, intent(out)   :: xngrid_new, yngrid_new, vertngrid_new, xngrid_coarse, yngrid_coarse
Real   , intent(in)    :: L05
Real   , intent(in)    :: mesh_size
Real   , intent(inout) :: xdomain_size, ydomain_size
Real   , intent(out)   :: intpl_coeff
!Local parameters
Real :: L02
Integer :: xngrid_coarse_tmp, yngrid_coarse_tmp       
Real :: xcut_partion, ycut_partion
Integer :: xngrid_est, yngrid_est, vertngrid_est
!===================Execituon========================
L02=2.5*L05
xcut_partion=xdomain_size/(xdomain_size+L02)
Print "(' ',a, f7.4)", 'xcut_partition=', xcut_partion
ycut_partion=ydomain_size/(ydomain_size+L02)
Print "(' ',a, f7.4)", 'ycut_partition=', ycut_partion
xngrid_est=nint(xngrid_old/xcut_partion)
yngrid_est=nint(yngrid_old/ycut_partion)
If (intpl_acceleration .eq. 1) Then
   intpl_coeff=max((L05/mesh_size)/comput_resolution, 1.)
   Print *, 'intpl_coeff=', intpl_coeff
   xngrid_coarse_tmp=nint(xngrid_est/intpl_coeff)
   yngrid_coarse_tmp=nint(yngrid_est/intpl_coeff)
End If
Call nearest235(xngrid_est, xngrid_new)
Call nearest235(yngrid_est, yngrid_new)
If (intpl_acceleration .eq. 1) Then
   Call nearest235(xngrid_coarse_tmp, xngrid_coarse)
   Call nearest235(yngrid_coarse_tmp, yngrid_coarse)
Else
   xngrid_coarse=xngrid_est
   yngrid_coarse=yngrid_est
End If
If (vertngrid_old .ne. 0) Then
   vertngrid_est=nint(vertngrid_old/xcut_partion)
   Call nearest235(vertngrid_est, vertngrid_new)
Else
   vertngrid_new=0
End If
xdomain_size=mesh_size*xngrid_new
ydomain_size=mesh_size*yngrid_new

End subroutine definegrid

Subroutine bilinear_interpolation(rand_field_coarse, rand_field_full,&
                                  xngrid_ext, yngrid_ext, vertngrid_ext, xngrid_big, yngrid_big, nstrides, nsample)
!=========================================================================================================
! PURPOSE
!
! We interpolate our field from the coarse grid to the full grid
!
! Author: D. Gayfulin, 2016
! Current code owner: D.Gayfulin
! Last updated: Nov, 2016
Implicit None
Integer, intent(in)  ::  xngrid_ext, yngrid_ext, vertngrid_ext, xngrid_big, yngrid_big, nstrides, nsample
Real, intent(in)  :: rand_field_coarse(0:xngrid_ext+1, 0:yngrid_ext+1, 0:vertngrid_ext-1, 1:nstrides, 1:nsample)
Real,  intent(Out) :: rand_field_full(0:xngrid_big-1, 0:yngrid_big-1, 0:vertngrid_ext-1, 1:nstrides, 1:nsample)
! local vars
Integer :: i, j, k, m, t
Integer :: x1, x2, y1, y2, x, y
Real :: wx1y1, wx1y2, wx2y1, wx2y2
!Real :: fx1y1, fx1y2, fx2y2, fx2y1
Real :: x_coeff, y_coeff
!===================Execituon========================
x_coeff=(xngrid_big-1.)/(xngrid_ext-1)
y_coeff=(yngrid_big-1.)/(yngrid_ext-1)
Print *, 'bilinear interpolation started'
If (vertngrid_ext .eq. 1) Then
   Do i=0, xngrid_big-1
      Do j=0, yngrid_big-1
         x1=floor(i/x_coeff); x2=x1+1;
         y1=floor(j/y_coeff); y2=y1+1;
         x=i/x_coeff-x1
         y=j/y_coeff-y1
         wx1y1=1-x-y+x*y
         wx1y2=x-x*y
         wx2y1=y-x*y
         wx2y2=x*y
          Do m=1, nsample
            Do t=1, nstrides
               rand_field_full(i, j, 0, t, m)=wx1y1*rand_field_coarse(x1,y1,0,t,m)+&
                                              wx1y2*rand_field_coarse(x1,y2,0,t,m)+&
                                              wx2y1*rand_field_coarse(x2,y1,0,t,m)+&
                                              wx2y2*rand_field_coarse(x2,y2,0,t,m)
            End Do  
         End Do
      End Do
   End Do
Else
   Do i=0, xngrid_big-1
      Do j=0, yngrid_big-1
         x1=floor(i/x_coeff); x2=x1+1;
         y1=floor(j/y_coeff); y2=y1+1;
         x=i/x_coeff-x1
         y=j/y_coeff-y1
         wx1y1=1-x-y+x*y
         wx1y2=x-x*y
         wx2y1=y-x*y
         wx2y2=x*y
        Do k=0, vertngrid_ext-1
            Do t=1, nstrides
               Do m=1, nsample
                  rand_field_full(i, j, k, t, m)=wx1y1*rand_field_coarse(x1,y1,k,t,m)+&
                                                 wx1y2*rand_field_coarse(x1,y2,k,t,m)+&
                                                 wx2y1*rand_field_coarse(x2,y1,k,t,m)+&
                                                 wx2y2*rand_field_coarse(x2,y2,k,t,m)
               End Do
            End Do
        End Do
      End Do
   End Do
End If
Print *, 'bilinear interpolation ended'
End Subroutine bilinear_interpolation

Subroutine nearest235(ngrid, ngridnew)
!=========================================================================================================
! PURPOSE
!
! we find the smallest number ngridnew, greater than ngrid such that ngridnew=2^n*3^m*5^l
!
! Author: D. Gayfulin, 2016
! Current code owner: D.Gayfulin
! Last updated: Feb, 2016
!
Implicit None
Integer, intent(in)  :: ngrid
Integer, intent(out) :: ngridnew
!
integer :: answ
!===================Execituon========================
ngridnew=ngrid
Call is235(ngridnew, answ)
Do While ((answ .ne. 1) .OR. (MOD(ngridnew, 2) .eq. 1))
  ngridnew=ngridnew+1
  Call is235(ngridnew, answ)
End do
End Subroutine nearest235

Subroutine is235(input, answer)
!=========================================================================================================
! PURPOSE
!
! we check, if input nomber has the form 2^n*3^m*5^l
!
! Author: D. Gayfulin, 2016
! Current code owner: D.Gayfulin
! Last updated: Feb, 2016
!
! 1 for true, otherwise for false
Implicit None
Integer, intent(in)  :: input
Integer, intent(out) :: answer
!
Integer :: curr_num
!===================Execituon========================
curr_num=input
Do While (MOD(curr_num, 2) .eq. 0)
  curr_num=curr_num/2
End Do
Do While (MOD(curr_num, 3) .eq. 0)
  curr_num=curr_num/3
End Do
Do While (MOD(curr_num, 5) .eq. 0)
  curr_num=curr_num/5
End Do
answer=curr_num
!===================Execituon========================
End Subroutine is235
end module generator
 

Real Function P_dn(n, m, l, q, lambda, mu, gamma, delta, stat)
!========================================
! This Function calculates the polynomial
! of the Laplacian operator
!========================================
Implicit None
Integer, intent(in)  :: n, m, l                      ! indexes of Fourier coefficient
Real,    intent(in)  :: q                            ! smoothness
Real   , intent(in)  :: lambda, mu, gamma, delta     ! spatio and temporal pseudoscales
Integer, intent(out) :: stat                         ! status
stat = 0
P_dn = mu * ((1. + (n**2+(gamma*m)**2+(delta*l)**2)*lambda**2)**q)   ! polynomial 
End Function P_dn
