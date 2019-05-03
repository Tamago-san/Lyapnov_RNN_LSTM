!-----------------------------------------------------------------------------
!ローレンツモデルのリアプノフ指数を求める。
!
!モデルのタイムスケールの0.01で進む。
!gfortran How-Lyapnov.f90 rc_tanh.f90 -o Lyapnov.out -llapack -lblas
!-----------------------------------------------------------------------------
module lyapnov
    implicit none
    integer, parameter :: rc_node = 5
    integer, parameter :: par_node = 3
!    real(8), parameter :: par_a = 16.d0
!    real(8), parameter :: par_b = 40.d0
!    real(8), parameter :: par_c = 4.d0
    real(8), parameter :: par_a = 10.d0
    real(8), parameter :: par_b = 28.d0
    real(8), parameter :: par_c = 8.d0/3.d0
    real(8), parameter  :: dt=1.d-2
    real(8), parameter :: dt_Runge = 1.d-4
    real(8), parameter :: dx0 = 1.d0
    real(8), parameter :: dy0 = 4.d0
    real(8), parameter :: dz0 = 2.d0
    real(8), parameter :: x0 = 0.d0
    real(8), parameter :: y0 = 0.1d0
    real(8), parameter :: z0 = 0.d0
    integer, parameter :: lyapnov_step = 2000 !リヤプノフ指数計測回数
    integer, parameter :: tau = 20
    integer, parameter :: rc_step = lyapnov_step * tau
    integer, parameter :: skip_step = 50000
    integer, parameter :: Ly_skip_step = 0
    integer, parameter :: traning_step = 5000
    integer, parameter :: god_step = 10
    integer, parameter :: NOW = 2
    integer, parameter :: BEFORE = 1
    integer, parameter :: in_node = 1
    integer, parameter :: out_node = 1
    real(8), parameter :: G0 = 1.d0
    real(8), parameter :: NU0 = 1.d0
    real(8), parameter :: gusai0 = 0.d0
    real(8) w_out(rc_node,out_node)
    real(8) w_in (in_node,rc_node)
    real(8) W_rnn(rc_node,rc_node)

  contains
    subroutine march(r,dr,z,U_rc,S_rc,abs_dr,fstep,ftau,gusai,alpha,G,allstep)
        real(8) U_rc(:,:)
        real(8) S_rc(:,:)
        real(8) r(0:tau,1:rc_node)
        real(8) dr(0:tau,1:rc_node)
        real(8) z(0:tau,1:rc_node)
        real(8) u_tmp(1:in_node)
        real(8) r_tmp(1:rc_node)
        real(8) DF(rc_node, rc_node)
        real(8) wu(rc_node),wr(rc_node)
        real(8) wu2,wr2,abs_dr
        real(8) gusai,alpha,G
        integer i,j,k,f,fstep,ftau,step,allstep
!        write(*,*) dr_befor(1:par_node)
        step=(fstep-1)*tau +ftau
        wu=0.d0
        wr=0.d0
        do i=1,rc_node
        do j=1,in_node
            wu(i)=wu(i)+u_rc(step,j)*W_in(j,i)
        enddo
        enddo
        do i=1,rc_node
        do j=1,rc_node
            wr(i)=wr(i)+r(ftau-1,j)*W_rnn(j,i)
        enddo
        enddo
        do i=1,rc_node
            z(ftau,i) =G*( wr(i)+wu(i)+gusai )
            r(ftau,i)=tanh(z(ftau,i))
        enddo
        do i = 1,out_node
		    s_rc(ftau,i)  = 0.d0
		    do k = 1,rc_node
		        s_rc(ftau,i) = s_rc(ftau,i) + r(ftau,k) * w_out(k,i)
		    enddo
		enddo
        do j=1,rc_node
        do i=1,rc_node
            DF(j,i)=(4.d0 /(( exp(z(ftau-1,j))+exp(-z(ftau-1,j) ) )**2) )*G*W_rnn(i,j)
!            DF(j,i)=(1.d0 /(cosh(z(ftau-1,j)) )**2 ) *G*W_rnn(i,j)
!            write(*,*) ftau,DF(i,j), dr(ftau-1,j),abs_dr
        enddo
        enddo
        do i=1,rc_node
            dr(ftau,i) = 0.d0
            do k=1,rc_node
!                write(*,*) ftau, dr(ftau-1,k),DF(k,i),dr(ftau-1,k)*DF(k,i)
                dr(ftau,i) = dr(ftau,i) + dr(ftau-1,k)*DF(k,i)
!            write(*,*) dr(ftau,i)
            enddo
        enddo
    end subroutine march
    
    subroutine calu_ERR(S_rc_out,S_rc_data,ERR)
        real(8) S_rc_data(1:tau,1:out_node)
        real(8) S_rc_out(1:tau,1:out_node)
        real(8) ERR
        integer fi,fj,fk
        ERR=0.d0
        do fi=1,tau
            do fj=1,out_node
                ERR=ERR+(S_rc_data(fi,fj)-S_rc_out(fi,fj))**2
            enddo
        enddo
        ERR=(ERR/dble(out_node*tau))**0.5d0
    end subroutine calu_ERR

    subroutine create_dataset(dataset,trmax)
        real(8) dataset(:,:)
        real(8) X(3,1:3)
        integer istep,trmax
        X(BEFORE,:)=1.d0
!        open(10,file='./data/output_Runge_Lorenz.dat')
        write(*,*) "=========================================="
        write(*,*) "START   CREATE TRANING SET"
        write(*,*) "=========================================="
        do istep = 1, skip_step
            call Runge_Kutta_method(X(NOW,1:3),X(BEFORE,1),X(BEFORE,2),X(BEFORE,3),istep)
!            if(mod(istep,10) == 0) write(10,*) istep ,x(NOW,1),x(NOW,2),x(NOW,3)
!            write(10,*) istep ,x(NOW,1),x(NOW,2),x(NOW,3)
            x(BEFORE,1:3) = x(NOW,1:3)
        enddo
!        close(10)
        do istep=1,trmax
            call Runge_Kutta_method(X(NOW,1:3),x(BEFORE,1),x(BEFORE,2),x(BEFORE,3),istep)
            x(BEFORE,1:3) = X(NOW,1:3)
            dataset(istep,1:3)= X(NOW,1:3)
        enddo
        call standard(dataset,3,trmax)
    end subroutine create_dataset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine standard(a,data_node,step)
		real(8) a(:,:)
		integer data_node,step
		real(8) a_mean(data_node),a_var(data_node)
		integer i,j
		do i = 1,data_node
			a_mean(i)=mean(a(1:step,i),step)
		enddo
		do i = 1,data_node
			a_var(i)=variance(a(1:step,i),a_mean(i),step)
		enddo
		do j = 1,data_node
			do i=1,step
				a(i,j)=(a(i,j)-a_mean(j))/a_var(j)
			enddo
		enddo
	end subroutine standard
	function mean(a,step) result(out)
		real(8), intent(in) :: a(:)
		real(8) :: out
		integer i,step
		out=0.d0
		do i=1,step
			out = out + a(i)
		enddo
		out= out/dble(step)
	end function mean
	function variance(a,a_mean,step) result(out)
		real(8), intent(in) :: a(:),a_mean
		real(8) :: out
		integer i,step
		out=0.d0
		do i=1,step
			out = out + (a(i)-a_mean)**2
		enddo
		out= (out/dble(step))**0.5
	end function variance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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

    subroutine Initialization_file
        open(20,file='./data_out/lyapnov.dat',status='replace')
        open(21,file="tmp.dat",status='replace')
        open(10,file='./data/output_Runge_Lorenz.dat',status='replace')
        open(22,file="./data_out/rc_out.dat",status='replace')
        close(10)
        close(20)
        close(21)
        close(22)
    end subroutine Initialization_file
    subroutine Ly_Calculation(abs_dr,lyapnov,step)
        real(8),intent(inout) :: lyapnov
        real(8),intent(in   ) ::  abs_dr
        real(8) ly_be,ly_now
        integer step
        integer fi
!        write(*,*) step,"Ly_Calculation",lyapnov,abs_dr
        ly_be  =lyapnov
        ly_now =log(abs_dr)/dble(tau)
!        write(*,100) step,"  Ly> before,now",tmp1,tmp2
        if(step>ly_skip_step) lyapnov =( ly_now  +ly_be*dble(step-1))/dble(step)
!            open(20,file='lyapnov.dat',position='append')
!            write(20,*) step_now,lyapnov
!            close(20)
100 format(i5,a,f15.10,f15.10)
    end subroutine Ly_Calculation
    function sigmoid(x)
        real(8) x, sigmoid   ! 順番は関係ありません
        sigmoid = 1.d0 / (1.d0 + exp(-x) )
    end function sigmoid
    
    function d_sigmoid(x)
        real(8) x, d_sigmoid   ! 順番は関係ありません
        d_sigmoid = sigmoid(x)*(1-sigmoid(x))
    end function d_sigmoid
    
    function d_tanh(x)
        real(8) x, d_tanh   ! 順番は関係ありません
        d_tanh = 4.d0 /(( exp(x)+exp(-x) )**2)
    end function d_tanh
    
    subroutine calu_absdr(abs_dr,dr,ftau)
        real(8) abs_dr
        real(8),intent(in):: dr(0:tau,1:rc_node)
        integer fi,ftau
        abs_dr=0.d0
        do fi=1,rc_node
!            write(*,*) fi,dr(ftau,fi)
            abs_dr = abs_dr+dr(ftau,fi)**2
        enddo
        abs_dr=abs_dr**0.5d0
    end subroutine calu_absdr
end module lyapnov

program main
use lyapnov
implicit none
    real(8) U_tr(traning_step,in_node) !今は一次元、列サイズはトレーニング時間
    real(8) S_tr(traning_step,out_node)  !出力次元数、列サイズはトレーニング時間
    real(8) U_rc(rc_step,in_node) !今は一次元、列サイズはトレーニング時間
    real(8) S_rc(rc_step,out_node)  !出力次元数、列サイズはトレーニング時間
    real(8) S_rc_out(tau,out_node)  !出力次元数、列サイズはトレーニング時間
    real(8) U_skip(Ly_skip_step,in_node) !今は一次元、列サイズはトレーニング時間
    real(8) S_skip(Ly_skip_step,out_node)  !出力次元数、列サイズはトレーニング時間
    real(8) dr0(rc_node), lyap,abs_dr,abs_dr2
    real(8) DF(rc_node, rc_node)
    real(8) dataset(1:rc_step+traning_step+god_step+100,par_node)
    real(8) r(0:tau,rc_node)
    real(8) dr(0:tau,rc_node)
    real(8) z(0:tau,rc_node)
    real(8) G_tmp,NU_tmp,ERR,ERR_tmp,lyapnov_now
    real(8) ERR_min(3)
    integer step,all_step
    integer istep,itau,iG,iNU,j,k,i,n
    character(6) :: cist1,cist2
    
    
    open(42,file="./data_out/lyapnov_end_ly.dat",status='replace')
    open(43,file="./data_out/lyapnov_end_err.dat",status='replace')
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!■create_dataset
!-------------------------------------------------
    all_step=ly_skip_step+rc_step+traning_step+god_step
    call create_dataset(dataset,all_step+100)
    U_tr(1:traning_step,1 )=dataset(1:traning_step,1 )
    s_tr(1:traning_step,1)=dataset(1+god_step:traning_step+god_step,1)
    U_rc(1:rc_step,1 )=dataset(1+traning_step:traning_step+rc_step,1)
    S_rc(1:rc_step,1)=dataset(1+traning_step+god_step:traning_step+rc_step+god_step,1)
    !open(10,file='./data/output_Runge_Lorenz.dat')
    !do i=1,3000
    !write(10,*) i,U_tr(i,1),S_tr(i,1)
    !enddo
    !close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ERR_min=1.d0
    do iNU=0,20
    do iG=1,20
    G_tmp =dble(iG)*0.1d0
    NU_tmp=dble(iNU)*0.1d0
!    G_tmp =G0
!    NU_tmp=NU0
    
    write(cist1,'(i6.6)') nint(G_tmp*100000.d0)+nint(NU_tmp*100.d0)
!    write(cist2,'(i3.3)') nint(NU_tmp*100)
!    write(*,*) int(G_tmp*1000)
    open(41,file='./data_renban/lyapnov.'//cist1,status='replace')
    open(51,file='./data_renban2/rc_out.'//cist1,status='replace')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!■初期化
!-------------------------------------------------
    dr=0.d0
    lyap=0.d0
    abs_dr=0.d0
    call Initialization_file
    call random_number(dr(0,:))
!    write(*,*) dr(0,:)
!    call random_seed
!    call random_number(dr(0,:))
!    write(*,*) dr(0,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(40,file="./data_out/lyapnov.dat")
    open(30,file="./data_out/rc_out.dat")
    open(31,file="tmp.dat")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!■Wout calu
!-------------------------------------------------
!    call calu_wout(U_tr,S_tr)
    call rc_traning_own_fortran(in_node,out_node,rc_node,traning_step,tau,gusai0,1.d0,G_tmp,NU_tmp,&
                    u_tr,s_tr,U_rc,w_in,w_rnn,w_out)
!    do i=1,rc_node
!        write(*,*) w_rnn(10,i)
!    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!■Lyaapnov calu
!-------------------------------------------------
!write(*,*) "=========================================="
!write(*,*) "START   LY SKIP STEP"
!write(*,*) "=========================================="
!    abs_dr=0.d0
!    do i=1,rc_node
!        abs_dr = abs_dr + dr(0,i)**2
!    enddo
!    abs_dr=abs_dr**0.5d0
!    write(*,*) "abs 0 === ",abs_dr
!    r(0,:) =1.d0
!    dr(0,:) = dr(0,:)/ abs_dr
!    do istep = 1,Ly_skip_step
!        call march(r,dr,U_skip,S_skip,abs_dr,1,istep,gusai0,1.d0,G_tmp,Ly_skip_step)
!        call calu_absdr(abs_dr,dr,istep)
!    enddo
!    call system("sleep 2")


    write(*,*) "=========================================="
    write(*,*) "START   CALUCULATE LY"
    write(*,*) "=========================================="
!    dr(0,:) =1.d0
    r(0,:) =0.d0
    z(0,:) =0.d0
    call calu_absdr(abs_dr,dr,0)
    dr(0,:) = dr(0,:)/ abs_dr
!    write(*,*) "dr 0 >.... ",dr(0,:)
!    call calu_absdr(abs_dr,dr,0)
    do istep = 1, lyapnov_step
        dr(0,:) = dr(0,:)/ abs_dr
!        write(*,*) dr(0,:)
!        write(*,*) dr(1:3)
        do itau = 1,tau
            step = (istep-1)*tau +itau
            call march(r,dr,z,U_rc,S_rc_out,abs_dr,istep,itau,gusai0,1.d0,G_tmp,rc_step)
            call calu_absdr(abs_dr,dr,itau)
            !abs_dr=0.d0
            !do i=1,rc_node
            !    abs_dr = abs_dr + dr(itau,i)**2
            !enddo
            !abs_dr=abs_dr**0.5d0
            write(31,'(3e14.6)') abs_dr ,dr(itau,5)
        enddo
        call Ly_Calculation(abs_dr,lyap,istep)
        !write(*,*) 'aaaaa',lyap,lyapnov_now
        call calu_ERR(S_rc_out(1:tau,1:out_node),S_rc((istep-1)*tau+1:istep*tau,1:out_node),ERR_tmp)
        ERR=(ERR*dble(istep-1)+ERR_tmp)/dble(istep)
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*) "---------------------------------------------------------"
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*)  istep ,int(istep*100/lyapnov_step),'%',int(iG),'%'
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*)  "G          ====", G_tmp
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*)  "NU         ====", NU_tmp
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*)  "abs_dr     ====", abs_dr
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*)  "lyapnov    ====", lyap
        if(mod(istep,nint(lyapnov_step*1.d0))==0) write(*,*)  "ERROR      ====", ERR
        r(0,1:rc_node)  = r(tau,1:rc_node)
        dr(0,1:rc_node) = dr(tau,1:rc_node)
        z(0,1:rc_node) = z(tau,1:rc_node)
        write(40,*) G_tmp,log(lyap),abs_dr,ERR
        write(41,*) istep,lyap,ERR
        do i=1,tau
            if(istep>=lyapnov_step-7) write(51,*) S_rc_out(i,1), S_rc((istep-1)*tau+i,1),U_rc((istep-1)*tau+i,1)
!            if(istep==lyapnov_step/2) write(51,*) dataset(1+traning_step+i,1:3)
        enddo
    enddo
    write(42,*) G_tmp,NU_tmp,lyap
    write(43,*) G_tmp,NU_tmp,ERR
    close(30)
    close(31)
    close(40)
    close(41)
    close(51)
    enddo
    write(42,*) " "
    write(43,*) " "
    enddo
    
200 format(a,f15.10,f15.10)
    close(42)
    close(43)
end program