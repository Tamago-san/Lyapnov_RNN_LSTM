!-----------------------------------------------------------------------------
!ローレンツモデルのリアプノフ指数を求める。
!
!モデルのタイムスケールの0.01で進む。
!-----------------------------------------------------------------------------
module lyapnov
    implicit none
    integer, parameter :: par_node = 3
!    real(8), parameter :: par_a = 16.d0
!    real(8), parameter :: par_b = 40.d0
!    real(8), parameter :: par_c = 4.d0
    real(8), parameter :: par_a = 10.d0
    real(8), parameter :: par_b = 28.d0
    real(8), parameter :: par_c = 8.d0/3.d0
    real(8), parameter  :: dt=1.d-5
    real(8), parameter :: dt_Runge = 1.d-7
    real(8), parameter :: dx0 = 1.d0
    real(8), parameter :: dy0 = 4.d0
    real(8), parameter :: dz0 = 2.d0
    real(8), parameter :: x0 = 2.d0
    real(8), parameter :: y0 = 10.d0
    real(8), parameter :: z0 = 5.d0
    integer, parameter :: itrmax = 500000
    integer, parameter :: tau = 100
    integer, parameter :: measur_step = itrmax * tau
    integer, parameter :: skip_time = 10
    integer, parameter :: NOW = 2
    integer, parameter :: BEFORE = 1
  contains
    subroutine Jacobi(DF,x,y,z)
        real(8) DF(:,:)
        real(8),intent(in) :: x,y,z

        DF(1,1) = -par_a
        DF(2,1) = par_b - z
        DF(3,1) = y
        DF(1,2) = par_a
        DF(2,2) = -1.d0
        DF(3,2) = x
        DF(1,3) = 0.d0
        DF(2,3) = -x
        DF(3,3) = -par_c
    end subroutine Jacobi
    
!-----------------------------------------------------------------------------
!■オイラー法での近似
!-----------------------------------------------------------------------------
    subroutine Euler_method(r_now,x,y,z,step)
        real(8) r_now(1,1:3)
        real(8) x,y,z
        integer step
        
        r_now(1,1) = x + (-par_a*(x - y)) * dt
        r_now(1,2) = y + ((par_b - z)*x - y) *dt
        r_now(1,3) = z + (x*y -par_c*z) *dt
        open(20,file='output_001_Euler_x2y10z5.dat',position='append')
        write(20,*) step ,r_now(1,1:3)
        close(20)
    end subroutine Euler_method
!======================================================
!-----------------------------------------------------------------------------
!■ルンゲクッタ法での近似
!-----------------------------------------------------------------------------
	subroutine Runge_Kutta_method(r_now,x,y,z,step)
		real(8) x,y,z
		real(8) r_now(1,1:3)
		real(8) k1(1:3),k2(1:3),k3(1:3),k4(1:3)
		real(8) k(1:3)
		integer i,j,step
		
		
        do j=1, int(dt/dt_Runge)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k1(1)=f1(x,y)
		    k1(2)=f2(x,y,z)
		    k1(3)=f3(x,y,z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k2(1)=f1(x+k1(1)*0.5d0 , y+k1(2)*0.5d0)
		    k2(2)=f2(x+k1(1)*0.5d0, y+k1(2)*0.5d0 , z+k1(3)*0.5d0)
		    k2(3)=f3(x+k1(1)*0.5d0, y+k1(2)*0.5d0 , z+k1(3)*0.5d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k3(1)=f1(x+k2(1)*0.5d0, y+k2(2)*0.5d0)
		    k3(2)=f2(x+k2(1)*0.5d0, y+k2(2)*0.5d0 , z+k2(3)*0.5d0)
		    k3(3)=f3(x+k2(1)*0.5d0, y+k2(2)*0.5d0 , z+k2(3)*0.5d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k4(1)=f1(x+k3(1), y+k3(2))
		    k4(2)=f2(x+k3(1), y+k3(2) , z+k3(3) )
		    k4(3)=f3(x+k3(1), y+k3(2) , z+k3(3))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		    do i=1,3
		    	k(i)=( k1(i) + 2.d0*k2(i) + 2.d0*k3(i) + k4(i) )/6.d0
		    enddo
		    x = x+k(1)
		    y = y+k(2)
		    z = z+k(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
		r_now(1,1) = x
		r_now(1,2) = y
		r_now(1,3) = z
!        open(20,file='output_001_Euler_x2y10z5.dat',position='append')
!        write(20,*) step ,r_now(1,1:3)
!        close(20)
    contains
        function f1(x,y)
	    	real(8) x,y
	    	real(8) f1
	    	f1=dt_Runge*(-par_a*(x-y))
	    end function f1
	    function f2(x,y,z)
	    	real(8) x,y,z
	    	real(8) f2
	    	f2=dt_Runge*((par_b-z)*x - y)
	    end function f2
	    function f3(x,y,z)
	    	real(8) x,y,z
	    	real(8) f3
	    	f3=dt_Runge*(x*y -par_c*z)
	    end function f3
	end subroutine Runge_Kutta_method
!======================================================
!======================================================

    subroutine march(dr,DF)
        real(8) dr(:),DF(:,:)
        real(8) tmp
        real(8) dr_befor(par_node)
        integer i,j,k
        dr_befor(1:par_node) = dr(1:par_node)
!        write(*,*) dr_befor(1:par_node)

        do i=1,par_node
            tmp = 0.d0
            do j=1,par_node
                tmp = tmp + DF(i,j)*dr_befor(j)
            enddo
            dr(i) = dr_befor(i) +dt*tmp
        enddo
    end subroutine march

    subroutine output(dr,lyapnov,step)
        real(8),intent(in) :: dr(:)
        real(8) lyapnov
        real(8) dr_tau
        integer step
        integer step_now
            step_now=step*tau
            dr_tau = (dr(1)**2 +dr(2)**2+dr(3)**2 )**(0.5d0)
!            dr_tau = (dr(3)**2 )**(0.5d0)
            lyapnov = ( lyapnov*dble(step-1) +(log(dr_tau))/(dble(tau)*dt) )/dble(step)
!            lyapnov = log(dr_tau)/(tau*dt)
            open(20,file='lyapnov.dat',position='append')
            write(20,*) step_now,lyapnov
            close(20)
    end subroutine output
        
end module lyapnov

program main
use lyapnov
implicit none
    real(8) dr(par_node)
    real(8) dr0(3), lyap,abs_dr
    real(8) DF(par_node, par_node)
    real(8) r(3,par_node),b_syoki
    real(8) skip_step
    integer step
    integer istep,j,k,i
    skip_step=dble(skip_time)/dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!■初期化
!-------------------------------------------------
    open(20,file='lyapnov.dat',status='replace')
    open(21,file='lyapnov_end.dat',status='replace')
    open(10,file='output_001_Euler_x2y10z5.dat',status='replace')
    open(30,file="a.dat",status='replace')
    close(10)
    close(20)
    close(30)
    lyap = 0.d0
    dr(1) = dx0
    dr(2) = dy0
    dr(3) = dz0
    r(BEFORE,1) = x0
    r(BEFORE,2) = y0
    r(BEFORE,3) = z0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!s
    open(30,file="a.dat")
    do istep = 1, itrmax
        abs_dr = (dr(1)**2 +dr(2)**2 +dr(3)**2)**(0.5d0)
        dr(1:3) = dr(1:3)/ abs_dr
!        write(*,*) dr(1:3)
!        dr=1.d0/(3.d0)**(0.5d0)
        do j = 1,tau
            step = (istep-1)*tau +j
            !----------------------------------------------------------------------------
!            call       Euler_method(r(step,1:3),r(step-1,1),r(step-1,2),r(step-1,3),step)
            call Runge_Kutta_method(r(NOW,1:3),r(BEFORE,1),r(BEFORE,2),r(BEFORE,3),step)
            !----------------------------------------------------------------------------
            if (step > skip_step) then
                call Jacobi(DF,r(NOW,1),r(NOW,2),r(NOW,3))
                call march(dr,DF)
            endif
            abs_dr = (dr(1)**2 +dr(2)**2 +dr(3)**2)**(0.5d0)
            write(30,*) step, log((abs_dr)**(0.5d0)) /(tau* dt),abs_dr,r(NOW,1)
            r(BEFORE,1:3) = r(NOW,1:3)
        enddo
        call output(dr,lyap,istep)
        if(mod(istep,100)==0) write(*,*) istep*tau,lyap,int(istep*tau*100/measur_step),'%'
    enddo
    close(30)
end program