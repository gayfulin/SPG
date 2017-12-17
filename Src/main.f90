Program main

Use generator, Only: generator_launch   
Implicit None
!====================================
!  This is the main program of the SPG. Here we read the external parameters from the config file,
!  call the subroutine that computes the pseudo-random spatio-temporal field, and calculate the diagnostic statistics
!  of the computed field
! 
!  The external parametrs are:
!  1)ngrid - number of gridpoints in each horizontal dimension.
!  2)vertngrid - number of gridpoints in the vertical dimension. 
!    !!!Note that for the two-dimensional version vertngrid should be equal to zero or one!!!
!  3)timeslice - time step in seconds, i.e. the time interval between the subsequent output spatial fields.
!  4)nstrides - number of time steps.
!  5)nsample - number of samples of fields generated. The fields from different samples are completely independent.
!  6)Std - field standard deviation.
!  7)L05 - The spatial correlation function of the field equals 0.5 at range L05.
!  8)T05 - The temporal correlation function of the field equals 0.2 at range T05.
!  9)a_dt_min - mininmal over the spectrum ratio step of integration/time scale. Should be less than 0.2.
!  10)a_dt_max - maximal over the spectrum ratio step of integration/time scale. Should be less than 3.
!  11)n_min_coarse_grid  - parameter n_0 of the coarse grid. The recommended value is 20 or more.
!  12)exp_factor_coarse_grid - parameter epsilon of the coarse grid. The recommended value is 0.1-0.2.
!  13)mesh_size - horizontal grid spacing (km).
!  14)dseed - seed for random numbers generation. Real number, should be greater than 1.
!
!  Authors: I.Mamay - 2012
!           D. Gayfulin - 2013-2016
!  Last update: October, 2016
!  Current code owner - D. Gayfulin, HMC of Russia, gayfulin@rambler.ru
!====================================
! External parameters:
! See the description above
!
Integer :: xngrid, yngrid, vertngrid
Integer :: n_min_coarse_grid
Integer :: nstrides, nsample
Integer :: calc_stat

Integer :: intpl_acceleration
Integer :: comput_resolution
!--------------------------------
Real :: timeslice
Real :: q
Real :: delta
Real :: std
Real :: L05, T05
Real :: a_dt_max, a_dt_min
Real :: mesh_size
Real :: exp_factor_coarse_grid
!-------------------------------
! Iternal parameters
!-------------------------------
! generated field
Real, Allocatable :: rand_field(:,:,:,:,:) ! dimensions: 1,2 - horysontal, 3 - vertical, 4 -time, 5-samples
!-------------------------------
! field parameters
Integer :: xhalfngrid, yhalfngrid, verthalfngrid ! 1/2 of resp. ngrid and vertngrid
!-------------------------------
! file access parameters
CHARACTER(Len=100):: filename ! config file
Integer, Parameter :: fileid=57              ! config file ID
!-------------------------------
! statistics parameters
Real :: avr ! variance of the whole field
Integer :: curr_range ! current distance between points (in km)
Real, Allocatable :: Corr_arr_x(:) ! array of spatial correlations 
Real, Allocatable :: Corr_arr_y(:) ! array of spatial correlations
Real, Allocatable :: Corr_arr_vert(:) ! array of spatial correlations
Real, Allocatable :: Corr_t(:)   ! array of temporal correlations 
Real, Allocatable :: var_arr(:)  ! array of variances for different time moments
Real, Allocatable :: num_l(:)    ! number of points, among which 
                                 ! the spatial correlation is being calculated for given range in space
Real, Allocatable :: num_l_y(:)    !
Real, Allocatable :: num_vert_l(:) ! number of vertical points, among which
                                 ! the spatial correlation is being calculated for given range in space
Real, Allocatable :: num_t(:)    ! number of points, among which
                                 ! the spatial correlation is being calculated for given time range
!-------------------------------
! other parameters
Integer :: n, m, l, i, k, j, t       ! loop counters
Integer :: ret_code, rtc, rtcc, stat ! return codes for subroutines
Double Precision  :: dseed           ! parameter for generating random values

Character(Len=10) :: ctime1, ctime2  ! CPU time

Real :: Computation_time
Character(Len=2)  :: c_hh1, c_hh2, c_mm1, c_mm2
Character(Len=6)  :: c_ss1, c_ss2
Real :: hh1, hh2, mm1, mm2, ss1, ss2

!INTEGER iargc,  getarg               ! input arguments parameters
!=========================================EXECUTION=============================================
!
!Opening config file
Print *, 'Stochastic Pattern Generator is launched'
Call Date_and_Time(time=ctime1)  ! time=hhmmss.sss
Print *, ''
IF( iargc() < 1) THEN
   Print *, 'Required parameters are:'
   Print *, ' 1. Config file name'
   Call Exit
END IF
Call getarg(1, filename)
Print *, 'Config file is ', Trim(filename)

Open(fileid, file = filename, form = 'formatted', &
Action='Read', Iostat=ret_code)
If (ret_code .eq. 0) Then
  Print *, 'Config file successfully opened'
Else
  Print *, 'Config file open failed, ret_code=', ret_code
  Call Exit(1)
End If
! reading config file
rtcc=0
Read(fileid, *, Iostat=rtc); rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) xngrid; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) yngrid; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) vertngrid; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) timeslice; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) nstrides; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) nsample; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) std; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) L05; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) T05; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) delta; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) a_dt_min; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) a_dt_max; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) n_min_coarse_grid; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) exp_factor_coarse_grid; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) intpl_acceleration; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) comput_resolution; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) calc_stat; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) mesh_size; rtcc=rtcc+rtc
Read(fileid, *, Iostat=rtc) dseed; rtcc=rtcc+rtc

If (xngrid .lt. 1) Then
   Print *, 'xngrid should be >=1, we set the default value 256'
   xngrid=256
End If

If (yngrid .lt. 1) Then
   Print *, 'yngrid should be >=1, we set the default value 256'
   yngrid=256
End If

If (vertngrid .lt. 0) Then
   Print *, 'vertngrid should be >=0, we set the default value 32'
   vertngrid=32
End If

If (timeslice .le. 0) Then
   Print *, 'timeslice should be >=0, we set the default value 3600'
   timeslice=3600
End If

If (nstrides .lt. 1) Then
   Print *, 'nstrides should be >=1,  we set the default value 12'
   nstrides=12
End If

If (nsample .lt. 1) Then
   Print *, 'nsample should be >=1,  we set the default value 4'
   nsample=4
End If

If (std .le. 0) Then
   Print *, 'std should be >0,  we set the default value 1'
   std=1
End If

If (L05 .le. 0) Then
   Print *, 'L05 should be >0,  we set the default value 100'
   L05=100
End If

If (T05 .le. 0) Then
   Print *, 'T02 should be >0,  we set the default value 7200'
   T05=7200
End If

If ((a_dt_min .le. 0) .OR. (a_dt_min .gt. 0.25)) Then
   Print *, 'a_dt_min should be between 0 and 0.25, we set the default value 0.05'
   a_dt_min=0.05
End If

If ((a_dt_max .lt. a_dt_min) .OR. (a_dt_max .gt. 5)) Then
   Print *, 'a_dt_max should be between a_dt_min and 5, we set the default value 3'
   a_dt_max=3
End If

If (n_min_coarse_grid .lt. 20) Then
    Print *, 'n_min_coarse_grid  should be >=20, we set the default value 20'
    n_min_coarse_grid=20
End If

If ((exp_factor_coarse_grid .le. 0.01) .OR. (exp_factor_coarse_grid .gt. 0.2)) Then
   Print *, 'exp_factor_coarse_grid should be between 0.01 and 0.2, we set the default value 0.1'
   exp_factor_coarse_grid=0.1
End If

If (mesh_size .le. 0) Then
   Print *, 'mesh_size should be >0,  we set the default value 7'
   mesh_size=7
End If

If (dseed .lt. 1) Then
   Print *, 'dseed should be >1,  we set the default value 405.7'
   dseed=405.7
End If
If (rtcc .eq. 0) Then
  Print *, 'Config file successfully read. The internal parameters are:'
  Print "(' ',a, i4, a)", 'x-dimension gridsize=', xngrid, ' points'
  Print "(' ',a, i4, a)", 'y-dimension gridsize=', yngrid, ' points'
  Print "(' ',a, i4, a)", 'vertical gridsize=', vertngrid, ' points'
  Print "(' ',a, f7.2, a)", 'time step=', timeslice, ' sec'
  Print "(' ',a, i4)", 'number of time steps=', nstrides
  Print "(' ',a, i4)", 'number of samples=', nsample
  Print "(' ',a, f6.3)", 'STD=', std
  Print "(' ',a, f8.2, a)", 'L05=', L05, ' km'
  Print "(' ',a, f9.2, a)", 'T05=', T05, ' sec'
  Print "(' ',a, f6.3)", 'delta=', delta
  Print "(' ',a, f6.3)", 'a_dt_min=', a_dt_min
  Print "(' ',a, f6.3)", 'a_dt_max=', a_dt_max
  Print "(' ',a, i4)", 'paramemter n_0 of the coarse grid=', n_min_coarse_grid
  Print "(' ',a, f6.3)", 'parameter epsilon of the coarse grid=', exp_factor_coarse_grid
  Print "(' ',a, i4)", 'interpolation acceleration, 1-yes, -1 -no=', intpl_acceleration
  Print "(' ',a, i4)", 'computational resolution=', comput_resolution
  Print "(' ',a, i4)", 'Calculate field statistics, 1-yes, -1 -no=', calc_stat
  Print "(' ',a, f6.3, a)", 'mesh_size=', mesh_size, ' km'
  Print "(' ',a, f16.8)", 'dseed=', dseed
  Print *, ''
Else
  Print *, 'Config file read failed, ret_code=', rtcc
  Call Exit(1)
End if
Close(fileid)
! defining parameters
q=0.5
exp_factor_coarse_grid=1+exp_factor_coarse_grid
xhalfngrid= Floor(xngrid*0.5)
yhalfngrid= Floor(yngrid*0.5)
verthalfngrid=Floor(vertngrid*0.5)
!Print *, 'verthalfngrid', verthalfngrid
!   main Call
Call Date_and_Time(time=ctime1)
Call generator_launch(xngrid, yngrid, vertngrid, q, std, L05, T05, delta, dseed, timeslice, nstrides, nsample,&
                      a_dt_min, a_dt_max, mesh_size, n_min_coarse_grid, exp_factor_coarse_grid,&
                      intpl_acceleration, comput_resolution, rand_field, stat)
If (stat .NE. 0) Then
  Call Exit(stat)
End If
Call Date_and_Time(time=ctime2)  ! time=hhmmss.sss
!*******************************************************************************************!
! The field rand_field is already calculated at this point and can be written to an output file !
!*******************************************************************************************!

! Now we calculate some statistics of the generated field
If (calc_stat .eq. 1) Then

!Call Date_and_Time(time=ctime2)
!Print *, ctime1, ctime2
! Allocating statistics arrays
avr=0
Allocate(Corr_arr_x(0:xhalfngrid))
Allocate(Corr_arr_y(0:yhalfngrid))
Allocate(Corr_arr_vert(0:vertngrid-2))
Allocate(Corr_t(0:nstrides-1))
Allocate(var_arr(1:nstrides))
Allocate(num_l(0:xhalfngrid))
Allocate(num_l_y(0:yhalfngrid))
Allocate(num_vert_l(0:vertngrid-2))
Allocate(num_t(0:xngrid))
num_t(:)=0; Corr_arr_x(:)=0; Corr_arr_y(:)=0; Corr_arr_vert(:)=0; Corr_t(:)=0; var_arr(:)=0; num_l(:)=0; num_vert_l(:)=0
num_l_y(:)=0
Print "(' ',a, 2f11.5)", 'Field minimum and maximum', Minval(rand_field(:,:,:,:,:)), Maxval(rand_field(:,:,:,:,:))
! Calculating variance for different time moments
If (vertngrid .eq. 0) vertngrid=1
Do t=1, nstrides
Do m=1, nsample
   Do i = 0, xngrid-1
      Do j = 0, yngrid-1
        Do k = 0, vertngrid-1
            var_arr(t)= var_arr(t)+rand_field(i,j,k,t,m)*rand_field(i,j,k,t,m)
         End Do
      End Do
    End Do
End do
End do
var_arr(:)=var_arr(:)/((nsample*xngrid*yngrid)*(vertngrid))
! Variance output
!Do i=1, nstrides
!   Print *, 'Std, time=', i*timeslice,'sec',Sqrt(var_arr(i))
!End Do
! Calculating the total variance
Do t=1, nstrides
Do m=1, nsample
   Do i = 0, xngrid-1
      Do j = 0, yngrid-1
        Do k = 0, vertngrid-1
            avr=avr+rand_field(i,j,k,t,m)*rand_field(i,j,k,t,m)
         End Do
      End Do
    End Do
End do
End do
! Variance output
avr=avr/(nsample*xngrid*yngrid*nstrides*vertngrid)
Print "(' ',a, f9.5)", 'Standard deviation of the field=',Sqrt(avr)
! Calculating the spatial correlation function
Print *, 'Estimating correlation functions'
Do l=0, xhalfngrid
    num_l(l)=0
Do m=1, nsample
   Do t=1, nstrides
      Do j = 0, yngrid-1
         Do k = 0, vertngrid-1
            i=0
            Do While (i+l+1 .le. xngrid-1) 
               Corr_arr_x(l)=Corr_arr_x(l)+rand_field(i+l+1,j,k,t,m)*rand_field(i,j,k,t,m)
               num_l(l)=num_l(l)+1
               i=i+l+1
            End do
         End do
      End Do
   End Do
   End do
End Do
Do l=0, yhalfngrid
    num_l_y(l)=0
Do m=1, nsample
   Do t=1, nstrides
      Do i = 0, xngrid-1
         Do k = 0, vertngrid-1
            j=0
            Do While (j+l+1 .le. yngrid-1)
               Corr_arr_y(l)=Corr_arr_y(l)+rand_field(i,j+l+1,k,t,m)*rand_field(i,j,k,t,m)
               num_l_y(l)=num_l_y(l)+1
               j=j+l+1
            End do
         End do
      End Do
   End Do
   End do
End Do
Do l=0, vertngrid-2
    num_vert_l(l)=0
Do m=1, nsample
   Do t=1, nstrides
      Do i = 0, xngrid-1
         Do j = 0, yngrid-1
            k=0
            Do While (k+l+1 .le. vertngrid-1)
               Corr_arr_vert(l)=Corr_arr_vert(l)+rand_field(i,j,k+l+1,t,m)*rand_field(i,j,k,t,m)
               num_vert_l(l)=num_vert_l(l)+1
               k=k+l+1
            End do
         End do
      End Do
   End Do
   End do
End Do
! Calculating the temporal correlation function
Do l=0, nstrides-2
num_t(l)=0
Do m=1, nsample
  Do j = 0, yngrid-1
     Do k = 0, vertngrid-1
        Do i = 0, xngrid-1
            t=1
            Do While (t+l+1 .le. nstrides)
                Corr_t(l)=Corr_t(l)+rand_field(i,j,k,t,m)*rand_field(i,j,k,t+l+1,m)
                num_t(l)=num_t(l)+1
                t=t+l+1
            End do
         End do
      End do
   End Do
End Do
End do


Print *, ''
Print *, ' Spatial correlation x'
Print *, '---------------------------------'
Print *, '|  Distance(km) | Correlation   |'
Print *, '---------------------------------'
Print "(' |',i5, '          |    ',f7.4,'    |')", 0, 1.0
Do l=0, nint(1.*xhalfngrid/2)
   curr_range=nint((mesh_size)*(l+1))  ! range in km in the domain, 
                                                     ! which corresponds to l points distance on torus
   Print "(' |',i5, '          |    ',f7.4,'    |')", curr_range, Corr_arr_x(l)/(avr*num_l(l))
End do
Print *, '---------------------------------'
Print *, ' '
Print *, ' Spatial correlation y'
Print *, '---------------------------------'
Print *, '|  Distance(km) | Correlation   |'
Print *, '---------------------------------'
Print "(' |',i5, '          |    ',f7.4,'    |')", 0, 1.0
Do l=0, nint(1.*yhalfngrid/2)
   curr_range=nint((mesh_size)*(l+1))  ! range in km in the domain,
                                                     ! which corresponds to l points distance on torus
   Print "(' |',i5, '          |    ',f7.4,'    |')", curr_range, Corr_arr_y(l)/(avr*num_l_y(l))
End do
Print *, '---------------------------------'
Print *, ' '

If (vertngrid .gt. 1) Then
Print *, ' Vertical correlation'
Print *, '---------------------------------'
Print *, '|   Distance,   | Correlation   |'
Print *, '|   gridpoints  |               |'
Print *, '---------------------------------'
Print "(' |',i5, '          |    ',f7.4,'    |')", 0, 1.0
Do l=0, nint(1.*vertngrid/2-1)
   curr_range=nint((mesh_size*xngrid/(1.*vertngrid))*(l+1))  ! range in km in the domain,
                                                            ! which corresponds to l points distance on torus
   Print "(' |',i5, '          |    ',f7.4,'    |')", l+1, Corr_arr_vert(l)/(avr*num_vert_l(l))
End do
Print *, '---------------------------------'
End If
Print *, ' '

Print *, 'Temporal correlation'

Print *, '---------------------------------'
Print *, '|  Time(hours)  | Correlation   |'
Print *, '---------------------------------'
Print "(' |',i5, '          |    ',f7.4,'    |')", 0, 1.0
Do l=0, nstrides-2
   Print "(' |',f7.2, '        |    ',f7.4,'    |')", (Nint((l+1)*timeslice))/3600., Corr_t(l)/(avr*num_t(l))
End do
Print *, '---------------------------------'
Print *, ' '

Deallocate(var_arr)
Deallocate(Corr_arr_x, Corr_arr_y, Corr_arr_vert)
Deallocate(Corr_t)
Deallocate(num_l, num_l_y)
Deallocate(num_t)
Deallocate(num_vert_l)
End If
Deallocate(rand_field)

!Call Date_and_Time(time=ctime2)  ! time=hhmmss.sss
c_hh1=ctime1(1:2); c_hh2=ctime2(1:2)
c_mm1=ctime1(3:4); c_mm2=ctime2(3:4)
c_ss1=ctime1(5:10); c_ss2=ctime2(5:10)
Read(c_hh1, *) hh1; Read(c_hh2, *) hh2
Read(c_mm1, *) mm1; Read(c_mm2, *) mm2
Read(c_ss1, *) ss1; Read(c_ss2, *) ss2
Computation_time=3600*(hh2-hh1)+60*(mm2-mm1)+ss2-ss1
Print "(' ',a,f8.3,a)", 'Total CPU time=', Computation_time, ' sec'
Print *, ' '
Print*, 'The program is successfully completed!'

End program main
