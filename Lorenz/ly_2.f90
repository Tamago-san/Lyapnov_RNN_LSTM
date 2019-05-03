!-----------------------------------------------------------------------------
!ローレンツモデルのリアプノフ指数を求める。
!
!モデルのタイムスケールの0.01で進む。
!-----------------------------------------------------------------------------
module lyapnov
    implicit none
    integer, parameter :: par_node = 3
    real(8), parameter :: par_a = 16.d0
    real(8) :: par_b = 40.d0
    real(8), parameter :: par_c = 4.d0
!    real(8), parameter :: par_a = 35.d0
!    real(8) :: par_b = 193.d0
!    real(8), parameter :: par_c = 8.d0/3.d0
    real(8), parameter  :: dt=1.d-4
    real(8), parameter :: dt_Runge = 1.d-6
    real(8), parameter :: dx0 = 1.d0
    real(8), parameter :: dy0 = 4.d0
    real(8), parameter :: dz0 = 2.d0
    real(8), parameter :: x0 = 0.d0
    real(8), parameter :: y0 = 0.1d0
    real(8), parameter :: z0 = 0.d0
    integer, parameter :: itrmax = 100000 !リヤプノフ指数計測回数
    integer, parameter :: tau = 100
    integer, parameter :: measur_step = itrmax * tau
    integer, parameter :: skip_time = 200
    integer, parameter :: Ly_skip_time = 500
    integer, parameter :: NOW = 2
    integer, parameter :: BEFORE = 1
  contains
    subroutine Jacobi(DF,x,y,z)
        real(8) DF(:,:)
        real(8),intent(in) :: x,y,z

        DF(1,1) = -par_a
        DF(2,1) = par_b - z
        DF(3,1) = y
        DF(1,2) = par_a +z
        DF(2,2) = -1.d0
        DF(3,2) = x
        DF(1,3) = y
        DF(2,3) = -x
        DF(3,3) = -par_c
    end subroutine Jacobi
    
!-----------------------------------------------------------------------------
!■オイラー法での近似
!-----------------------------------------------------------------------------
!    subroutine Euler_method(r_now,x,y,z,step)
!        real(8) r_now(1,1:3)
!        real(8) x,y,z
!        integer step
!
!        r_now(1,1) = x + (-par_a*(x - y) + y*z) * dt
!        r_now(1,2) = y + ((par_b - z)*x - y) *dt
!        r_now(1,3) = z + (x*y -par_c*z) *dt
!        open(20,file='output_001_Euler_x2y10z5.dat',position='append')
!        write(20,*) step ,r_now(1,1:3)
!        close(20)
!    end subroutine Euler_method
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
		    k1(1)=f1(x,y,z)
		    k1(2)=f2(x,y,z)
		    k1(3)=f3(x,y,z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k2(1)=f1(x+k1(1)*0.5d0 ,y+k1(2)*0.5d0 , z+k1(3)*0.5d0)
		    k2(2)=f2(x+k1(1)*0.5d0, y+k1(2)*0.5d0 , z+k1(3)*0.5d0)
		    k2(3)=f3(x+k1(1)*0.5d0, y+k1(2)*0.5d0 , z+k1(3)*0.5d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k3(1)=f1(x+k2(1)*0.5d0, y+k2(2)*0.5d0 , z+k2(3)*0.5d0)
		    k3(2)=f2(x+k2(1)*0.5d0, y+k2(2)*0.5d0 , z+k2(3)*0.5d0)
		    k3(3)=f3(x+k2(1)*0.5d0, y+k2(2)*0.5d0 , z+k2(3)*0.5d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		    k4(1)=f1(x+k3(1), y+k3(2) , z+k3(3))
		    k4(2)=f2(x+k3(1), y+k3(2) , z+k3(3))
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
        function f1(x,y,z)
	    	real(8) x,y,z
	    	real(8) f1
	    	f1=dt_Runge*(-par_a*(x-y)+ y*z)
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

    subroutine Ly_Calculation(dr,lyapnov,step)
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
!            open(20,file='lyapnov.dat',position='append')
!            write(20,*) step_now,lyapnov
!            close(20)
    end subroutine Ly_Calculation
        
end module lyapnov

program main
use lyapnov
implicit none
    real(8) dr(par_node)
    real(8) dr0(3), lyap,abs_dr
    real(8) DF(par_node, par_node)
    real(8) r(3,par_node),b_syoki
    real(8) skip_step,Ly_skip_step
    integer step
    integer istep,j,k,i
    skip_step=dble(skip_time)/dt
    Ly_skip_step = dble(Ly_skip_time)/dt
    b_syoki=par_b

    open(21,file='./data_out/lyapnov_end.dat',status='replace')
    close(21)
    
    
    do i=1,1
 !   par_b = b_syoki + dble(i)*0.01
    par_b = b_syoki
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!■初期化
!-------------------------------------------------
    open(20,file='./data_out/lyapnov.dat',status='replace')
    open(10,file='./data_out/output_Runge_Lorenz.dat',status='replace')
    open(30,file="./data_out/Ly_cal.dat",status='replace')
    open(11,file="./data_out/Ly_skip.dat",status='replace')
    close(10)
    close(11)
    close(20)
    close(30)
    lyap = 0.d0
    dr(1) = dx0
    dr(2) = dy0
    dr(3) = dz0
    r(BEFORE,1) = x0
    r(BEFORE,2) = y0
    r(BEFORE,3) = z0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(10,file='./data_out/output_Runge_Lorenz.dat')
    open(11,file='./data_out/Ly_skip.dat')
    open(30,file="./data_out/Ly_cal.dat")
    open(21,file='./data_out/lyapnov_end.dat',position='append')
write(*,*) "=========================================="
write(*,*) "START   SKIP STEP"
write(*,*) "=========================================="
    do istep = 1, int(skip_step)
        call Runge_Kutta_method(r(NOW,1:3),r(BEFORE,1),r(BEFORE,2),r(BEFORE,3),istep)
        if(mod(istep,10) == 0) write(10,*) dble(step)*dt*1.d2 ,r(NOW,1:3)
        r(BEFORE,1:3) = r(NOW,1:3)
    enddo
    close(10)
write(*,*) "=========================================="
write(*,*) "START   LY SKIP STEP"
write(*,*) "=========================================="
    abs_dr = (dr(1)**2 +dr(2)**2 +dr(3)**2)**(0.5d0)
    dr(1:3) = dr(1:3)/ abs_dr
    do istep = 1, int(Ly_skip_step)
        call Runge_Kutta_method(r(NOW,1:3),r(BEFORE,1),r(BEFORE,2),r(BEFORE,3),istep)
        call Jacobi(DF,r(NOW,1),r(NOW,2),r(NOW,3))
        call march(dr,DF)
        abs_dr = (dr(1)**2 +dr(2)**2 +dr(3)**2)**(0.5d0)
        if(mod(istep,100) == 0) write(11,*) istep,abs_dr,r(NOW,1),r(NOW,2),dr(1)
        r(BEFORE,1:3) = r(NOW,1:3)
    enddo
    close(11)
write(*,*) "=========================================="
write(*,*) "START   CALUCULATE LY"
write(*,*) "=========================================="
    do istep = 1, itrmax
        abs_dr = (dr(1)**2 +dr(2)**2 +dr(3)**2)**(0.5d0)
        dr(1:3) = dr(1:3)/ abs_dr
!        write(*,*) dr(1:3)
        do j = 1,tau
            step = (istep-1)*tau +j
            !----------------------------------------------------------------------------
!            call       Euler_method(r(step,1:3),r(step-1,1),r(step-1,2),r(step-1,3),step)
            call Runge_Kutta_method(r(NOW,1:3),r(BEFORE,1),r(BEFORE,2),r(BEFORE,3),step)
            !----------------------------------------------------------------------------
            call Jacobi(DF,r(NOW,1),r(NOW,2),r(NOW,3))
            call march(dr,DF)
            abs_dr = (dr(1)**2 +dr(2)**2 +dr(3)**2)**(0.5d0)
            r(BEFORE,1:3) = r(NOW,1:3)
        enddo
        write(30,*) istep,lyap,abs_dr,r(NOW,2),r(NOW,3),dr(1)
        call Ly_Calculation(dr,lyap,istep)
        if(mod(istep,100)==0) write(*,*)  par_c,istep*tau,lyap ,int(istep*tau*100/measur_step),'%'
    enddo
    close(30)
    write(21,*) par_b ,lyap
    close(21)
    enddo
end program