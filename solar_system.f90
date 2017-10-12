	!> @mainpage
	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>program to solve for the orbits of the planets
	!>Solar System Model (SSM): 
	!>Solves equations of motion in 3-D:
	!> <br><br>
	!> This program uses Newton's law of gravitation 
	!> which states that the force between two masses, 
	!> \f$m_i\f$ and \f$m_j\f$, separated by vector \f$\vec{r_{i,j}}\f$ is <br>
	!> \f{equation*}{ \vec{F}_{i,j}=-G\frac{m_i m_j}{r_{i,j}^3}\vec{r_{i,j}} \f}
	!><br><br>
	!> The coordinate system is cartesian based (i.e. x,y,z) 
	!> so Newton's Law of gravitation is written using this system. 
	!> The vector \f$\vec{r_{i,j}}\f$ can be written down in column form as: <br><br>
	!> \f{equation*}{ \vec{r_{i,j}}=\left(\begin{array}{c} x_j-x_i 
	!> \\y_j-y_i \\z_j-z_i \end{array} \right) \f}
	!> <br><br>
	!> Hence, the components of the force in the x,y and z direction are: <br>
	!>\f{equation*}{ \left(\begin{array}{c}
    !>  Fx_{i,j}\\
    !>  Fy_{i,j}\\
    !>  Fz_{i,j}
    !>  \end{array}\right)
    !> =-G\frac{m_i m_j}{\left({\left(x_j-x_i\right)^2 +
    !> \left(y_j-y_i\right)^2 +\left(z_j-z_i\right)^2}\right)^{3/2}}
    !>  \left(\begin{array}{c} x_j-x_i \\y_j-y_i \\z_j-z_i \end{array} \right) \f} <br>
    !> Newton's 2nd Law states for net force, \f$ \vec{F}\f$ , 
    !> applied to object mass, m causing the 
    !> object to accelerate with acceleration, a: <br>
    !> \f{equation*}{ \vec{F}=m\vec{a}\f} <br>
    !> Hence, to solve for the motion of each object (i.e. find the \f$x_i,y_i,z_i\f$ 
    !> positions 
    !> with time) we must solve three differential equations for each dimension: <br>
    !> \f{eqnarray*}{
    !> Fx_{i,j} &=& m_i\frac{d^2 x_i}{dt^2}\\
    !> Fy_{i,j} &=& m_i\frac{d^2 y_i}{dt^2}\\
    !> Fz_{i,j} &=& m_i\frac{d^2 z_i}{dt^2}
    !> \f} <br> <br>
    !> The code solves 1st order differential equations 
    !> (those having a first order derivative). Clearly the equations to solve are 
    !> 2nd order differential equations, but we turn them into 1st order differential 
    !> equations by specifying another 3 equations for the 3 components of velocity, 
    !> \f$vx_i,vy_i,vz_i\f$, which give the following 6 equations to solve for each 
    !> body in the solar system: <br><br>
    !> \f{eqnarray*}{
    !> Fx_{i,j} &=& m_i\frac{d vx_i}{dt}\\
    !> Fy_{i,j} &=& m_i\frac{d vy_i}{dt}\\
    !> Fz_{i,j} &=& m_i\frac{d vz_i}{dt}\\
    !> vx_i &=& \frac{d x_i}{dt}\\
    !> vy_i &=& \frac{d y_i}{dt}\\
    !> vz_i &=& \frac{d z_i}{dt}
    !> \f}
	!> compile using the Makefile (note requires netcdf) and then run using: <br>
	!> ./main.exe namelist.in
	!> <br><br>
	!> (namelist used for initialisation).
	!> <br><br>




	!>@author
	!>Paul Connolly, The University of Manchester
	!>@brief
	!>variables for the solar system model
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Module so variables can be easily passed to derivative function        !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	module consts_and_vars
		use, intrinsic :: iso_fortran_env
		
		implicit none
		integer, parameter :: dp = real64
		integer, parameter :: i4b = selected_int_kind(9)
		integer(i4b), parameter :: n_bodies=18
		integer(i4b) :: neq=n_bodies*6,nn_interact=n_bodies
		real(dp), parameter :: c=2.99792458e8_dp ! speed of light
		real(dp), dimension(n_bodies) :: gm, &
						 x,y,z, ux,uy,uz, meandist, lx,ly,lz,lt
		real(dp), dimension(n_bodies*6) :: yinit1
		integer(i4b), dimension(n_bodies) :: inds1, inds, &
													 interactions,interactions2

		logical, dimension(n_bodies) ::  interact=.true.
		logical :: general_relativity=.true.
		real(dp), parameter :: G=6.6719e-11_dp, m=1.98855e30_dp      
		real(dp) :: tt=0.e0_dp
	end module consts_and_vars
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>main program of Solar System Model (SSM): 
	program solar_system

	use netcdf
	use consts_and_vars, only : n_bodies, neq,gm,x,y,z,ux,uy,uz, yinit1, &
						 meandist, inds1,inds,G,m,interactions,tt, &
						 interact,nn_interact,interactions,interactions2, &
						 general_relativity, lx,ly,lz,lt, dp, i4b

	implicit none
	integer(i4b) :: i, j, allocatestatus
	real(dp), dimension(6*n_bodies) :: yinit, ydot, ysol, atol,rtol1
	real(dp), dimension(6*n_bodies) :: rpar
	integer(i4b), dimension(6*n_bodies) :: ipar

	real(dp) :: tout, dt=5e0_dp*86400.e0_dp, rtol, tfinal, interval_io, &
				tt_last
	integer(i4b) :: itol, itask, istate, iopt, ng, & 
								   lrw, liw, mf, mflag
	real(dp), allocatable, dimension(:) :: rwork 
	integer(i4b), allocatable, dimension(:) :: iwork

	! netcdf stuff
	character (len = 100) :: outputfile01='output02.nc'
	integer(i4b), parameter :: nsav=10000
	real(dp), dimension(3,2,n_bodies) :: pos
	real(dp), dimension(3,2,n_bodies,nsav) :: pos_save
	real(dp), dimension(nsav) :: tt_save
	integer(i4b) :: ncid, varid, x_dimid, i_dimid, j_dimid, & 
								 a_dimid, icur,isav
	logical :: run_forward_in_time = .true.

	! namelist
	character (len=200) :: nmlfile = ' '
	namelist /initial_state/ x,y,z,ux,uy,uz, meandist,gm
	namelist /run_vars/ outputfile01,run_forward_in_time, general_relativity, &
					  dt,tfinal, interval_io, interact

	external  solar01, jsolar01




	interactions=(/(i,i=1,n_bodies)/)

	tfinal=365.25e0_dp*86400.e0_dp*5.e2_dp






	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! namelist                                                               !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call getarg(1,nmlfile)      
	open(8,file=nmlfile,status='old', recl=80, delim='apostrophe')
	read(8,nml=initial_state)
	x=x*1e3_dp
	y=y*1e3_dp
	z=z*1e3_dp
	ux=ux*1e3_dp
	uy=uy*1e3_dp
	uz=uz*1e3_dp
	gm=gm*1e9_dp
	read(8,nml=run_vars)
	dt=dt*86400.e0_dp
	tfinal=tfinal*365.25e0_dp*86400.e0_dp
	close(8)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






	nn_interact=count(interact) ! number of true elements
	j=1
	do i=1,n_bodies 
	  if(interact(i)) then
		 interactions2(j)=interactions(i)
		 j=j+1
	  end if
	end do


	if(.not.run_forward_in_time) then
	 ux=-ux
	 uy=-uy
	 uz=-uz      
	endif

	! mean distance from the sun
!	meandist=(/7e8_dp, 5.79e10_dp, 1.082e11_dp,  &
!		  1.496e11_dp, 2.279e11_dp, 7.783e11_dp, &
!		 1.426e12_dp, 2.871e12_dp, 4.497e12_dp, 5.914e12_dp, &
!		 1.0159e13_dp, 7.0792e12_dp, 8.0725e12_dp, 5.8516e12_dp, &
!		  5.7248e12_dp, 5.4604e12_dp, 5.4460e12_dp, 5.3797e12_dp, &
!	   	5.5168e12_dp, 5.3324e12_dp, 1.2391e13_dp, 1.2932e13_dp/) 

	inds1=1

	! angular momenta per unit mass
	lx=uy*z-uz*y
	ly=uz*x-ux*z
	lz=ux*y-uy*x
	lt=sqrt(lx*lx+ly*ly+lz*lz)








	! open netcdf file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_create(outputfile01, nf90_clobber, ncid) )
	call check( nf90_def_dim(ncid, "times", nf90_unlimited, x_dimid) )
	call check( nf90_def_dim(ncid, "n_bodies", n_bodies, i_dimid) )
	call check( nf90_def_dim(ncid, "n_dims", 3 , j_dimid) )
	! close and free up any buffers
	call check( nf90_close(ncid) )

	! open file for writing:
	call check( nf90_open(outputfile01, nf90_write, ncid) )    
	! define mode
	call check( nf90_redef(ncid) )     



	! define a variable TIME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_def_var(ncid, "time", nf90_float, & ! was nf90_double
					  (/x_dimid/), varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(ncid, "time", a_dimid) )
	! units
	call check( nf90_put_att(ncid, a_dimid, "units", "s") )        


	! define a variable POS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_def_var(ncid, "pos", nf90_float, &
					  (/j_dimid, i_dimid,x_dimid/), varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(ncid, "pos", a_dimid) )
	! units
	call check( nf90_put_att(ncid, a_dimid, "units", "m") )        


	! define a variable VEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call check( nf90_def_var(ncid, "vel", nf90_float, &
					  (/j_dimid, i_dimid,x_dimid/), varid) )
	! get id to a_dimid
	call check( nf90_inq_varid(ncid, "vel", a_dimid) )
	! units
	call check( nf90_put_att(ncid, a_dimid, "units", "m s-1") )        

	! END DEFINE MODE. 
	! TELL netCDF WE ARE DONE DEFINING METADATA.
	call check( nf90_enddef(ncid) )
	call check( nf90_close(ncid) )
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Set-up variables for the ode solver                                    !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,n_bodies
	 yinit( 1+(i-1)*6: 6*i )=(/x(i), y(i), z(i), ux(i), uy(i), uz(i)/)
	end do
	yinit1=yinit


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! set up variables to pass to dvode                                      !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	tt=0.e0_dp
	ysol=yinit
	tout=tt+dt
	!      tfinal=365.25d0*86400.d0*1.d6
	!itol=2 ! both relative and absolute error convergence
	itol=4 ! both relative and absolute error convergence
	!      rtol=1.d-8
	rtol1=1.e-8_dp
	rtol1(7:12)=1e-10_dp
	rtol1=1d-12
	atol(:)=1.e-6_dp ! nearest micron and micron / second (i.e. negligible)
	itask=1
	istate=1
	iopt=1 ! optional input
	mf=10 !22 ! 
	lrw=22+9*neq+2*neq**2
	liw=30*neq

	allocate( rwork(lrw), stat = allocatestatus)
	if (allocatestatus /= 0) stop "*** not enough memory ***"

	allocate( iwork(liw), stat = allocatestatus)
	if (allocatestatus /= 0) stop "*** not enough memory ***" 

	iwork(6)=10000 ! max steps
	iwork(7)=10  ! max messages printed per problem
	iwork(5)=12 ! 5   ! order
	rwork(5)=0.e0_dp ! initial time-step
	rwork(6)=dt   ! max time-step
	rwork(7)=0.e0_dp ! min time-step allowed
	rwork(14)=2.e0_dp ! tolerance scale factor
	mflag=0
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!CALL XSETF(MFLAG) 





	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! solve system of equations calling dvode                                !
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	icur=1
	isav=1
	istate=1
	interval_io=interval_io*365.25_dp*86400._dp ! put in seconds
	tt_last=-interval_io
	do while (tt.lt.tfinal)
	 do while(tt.lt.tout)
	   call dvode(solar01,neq,ysol,tt,tout,itol,rtol1,atol,itask,istate, &
			  iopt,rwork,lrw,iwork,liw,jsolar01,mf,rpar,ipar)
	   istate=2
	 enddo
	 tout=tt+dt


	 tt_save(isav)=tt
	 pos_save(1:3,1:2,1:n_bodies,isav:isav)= &
		 reshape(ysol(1:neq),(/3,2,n_bodies,1/))

	 if(tt-tt_last >= interval_io) then
	 	 ! position in run:
	 	 print *,'Model is up to here: ',tt/(365.25*86400), &
	 	 		 ' done of ',tfinal/(365.25*86400)
	 	 tt_last=tt
	 endif

	 if (isav.eq.nsav.or.tt.ge.tfinal) then
	 	
		! open the netcdf file (do not overwrite if it exists)
		call check( nf90_open(outputfile01, nf90_write, ncid) )    

		! write netcdf variables
		call check( nf90_inq_varid(ncid, "time", varid ) )
		call check( nf90_put_var(ncid, varid, tt_save(1:isav), start = (/ICUR/) ) ) 


		call check( nf90_inq_varid(ncid, "pos", varid ) )
		call check( nf90_put_var(ncid, varid, pos_save(1:3,1,1:n_bodies,1:isav), &
		  start = (/1,1,ICUR/) ) ) 

		call check( nf90_inq_varid(ncid, "vel", varid ) )
		call check( nf90_put_var(ncid, varid, pos_save(1:3,2,1:n_bodies,1:isav), &
		 start = (/1,1,ICUR/) ) ) 

		! close the netCDF file:
		call check( nf90_close(ncid) )

		icur=icur+isav
		isav=0
	 end if

	 isav=isav+1

	end do
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	end program solar_system
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>provides derivatives for the RHS of the equations
	!>@param[inout] neq: number of odes to solve
	!>@param[inout] tt: time variable
	!>@param[inout] y, ydot: solution and derivative
	!>@param[inout] rpar, ipar: real and integer work arrays
	subroutine solar01(neq, tt, y, ydot, rpar, ipar)

	use consts_and_vars, only : gm, n_bodies,inds,interactions, yinit1, &
							  interact,interactions2,nn_interact,c, &
							  general_relativity,lx,ly,lz,lt, dp, i4b

	implicit none

	real(dp), intent(inout) :: tt
	real(dp), intent(inout), dimension(neq) :: y, ydot
	integer(i4b), intent(inout) :: neq
	real(dp), intent(inout),dimension(6*n_bodies) :: rpar
	integer(i4b), intent(inout),dimension(6*n_bodies) :: ipar

	! locals
	integer(i4b) :: i,n,k,j
	real(dp), dimension(n_bodies-1) :: r2
	logical, dimension(n_bodies) :: interact2
	real(dp) :: lx1,ly1,lz1, lt1

	ydot=0.d0


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Loop over all bodies                                                               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=1,n_bodies   
	 interact2=interact
 
	 n=count(interactions2(1:nn_interact).ne.i)
	 if (n.eq.0) cycle
 
 	 ! Do not interact with self:
	 if (i.eq.interactions(i)) interact2(i)=.false. 
	 
	 ! inds is the body that i interacts with:
	 inds(1:n)=pack(interactions,interact2)

	 ! inds indexes the body that body i interacts with.
	 ! however, the x,y,z's are separated by 6
	 
	 
	 ! square of distance between this planet and the other objects:
	 r2(1:n)=(y( (i-1)*6+1 )-y( (inds(1:n)-1)*6+1 ))**2 + &
		(y( (i-1)*6+2 )-y( (inds(1:n)-1)*6+2 ))**2 + &
		(y( (i-1)*6+3 )-y( (inds(1:n)-1)*6+3 ))**2





	 ! inverse square law between body i and the rest of them
	 ! dux/dt:
	 ydot((i-1)*6+4)=ydot((i-1)*6+4)- &
	  sum(gm( inds(1:n) )/ &
	  r2(1:n)*(y((i-1)*6+1)-y( (inds(1:n)-1)*6+1 ))/sqrt(r2(1:n)) )
	 ! duy/dt:
	 ydot((i-1)*6+5)=ydot((i-1)*6+5)- &
	  sum(gm( inds(1:n) )/ &
	  r2(1:n)*(y((i-1)*6+2)-y( (inds(1:n)-1)*6+2 ))/sqrt(r2(1:n)) )
	 ! duz/dt:
	 ydot((i-1)*6+6)=ydot((i-1)*6+6)- &
	  sum(gm( inds(1:n) )/ &
	  r2(1:n)*(y((i-1)*6+3)-y( (inds(1:n)-1)*6+3 ))/sqrt(r2(1:n)) )




	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! extra term for general relativity - do not apply to sun's motion                  !
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 if (general_relativity .and. (i.gt.1)) then
		! terms for General Relativity - 
		! http://farside.ph.utexas.edu/teaching/336k/Newtonhtml/node116.html
		n=1	
		! note that l^2=G*M0*r
		! dux/dt:
		ydot((i-1)*6+4)=ydot((i-1)*6+4)- &
		sum(3.e0_dp*gm( inds(1:n) )**2/c**2 / &
		r2(1:n)**2*(y((i-1)*6+1)-y( (inds(1:n)-1)*6+1 )) )
		! duy/dt:
		ydot((i-1)*6+5)=ydot((i-1)*6+5)- &
		sum(3.e0_dp*gm( inds(1:n) )**2/c**2 / &
		r2(1:n)**2*(y((i-1)*6+2)-y( (inds(1:n)-1)*6+2 )) )
		! duz/dt:
		ydot((i-1)*6+6)=ydot((i-1)*6+6)- &
		sum(3.e0_dp*gm( inds(1:n) )**2/c**2 / &
		r2(1:n)**2*(y((i-1)*6+3)-y( (inds(1:n)-1)*6+3 )) )
	 end if
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





	 ! set the rate of change of position to the y-values in the 2nd order diff eq
	 ydot( (i-1)*6+1 ) = y( (i-1)*6+4 )
	 ydot( (i-1)*6+2 ) = y( (i-1)*6+5 )
	 ydot( (i-1)*6+3 ) = y( (i-1)*6+6 )
	end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
#ifdef heliocentric
	! adjust positions and accelerations, so the sun remains at (0,0,0)
	ydot(7:n_bodies*6) = ydot(7:n_bodies*6)- &
		  reshape(spread( ydot(1:6),2,n_bodies-1),(/6*(n_bodies-1)/))
	ydot(1:6) = 0.d0
#endif
	end subroutine solar01
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>jacobian needs to be defined - does not need to do anything
	!>@param[inout] neq: number of odes to solve
	!>@param[inout] tt: time variable
	!>@param[inout] y, ydot: solution and derivative
	!>@param[inout] rpar, ipar: real and integer work arrays
	subroutine jsolar01(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

		use consts_and_vars, only : n_bodies, dp, i4b

		implicit none
		real(dp), intent(in) :: t
		real(dp), dimension(neq), intent(inout) :: y
		real(dp), dimension(nrpd, neq), intent(inout) :: pd
		integer(i4b), intent(inout) :: neq, ml, mu, nrpd
		real(dp), intent(inout),dimension(6*n_bodies) :: rpar
		integer(i4b), intent(inout),dimension(6*n_bodies) :: ipar



	end subroutine jsolar01
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!>@author
	!>Paul J. Connolly, The University of Manchester
	!>@brief
	!>helper routine for netcdf   
	!>@param[in] status of netcdf file pointer
	subroutine check(status)
		use netcdf
		use consts_and_vars, only : dp, i4b
		integer(i4b), intent (in) :: status

		if(status /= nf90_noerr) then 
		   print *, trim(nf90_strerror(status))
		   stop "stopped"
		end if
	end subroutine check  



