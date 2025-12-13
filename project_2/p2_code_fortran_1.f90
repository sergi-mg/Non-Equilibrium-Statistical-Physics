!Author: Sergi Martínez Galindo
!2nd Practice Non-Equilibrium Statistical Physics (2,3)
program practica2
    implicit none
    integer :: N, Np,i,i_2,j,k,m,k2
    double precision :: L,D,dt,tau,r2,genrand_real3,gamma,x_0(1:2),tau_vec(1:3),&
    t_max,C
    double precision, dimension(:), allocatable :: vector_var,val_t,MSD,delta&
    ,vector_PBC
    double precision, dimension(:,:), allocatable :: correlacio, correlacio_j
    external :: zero, pas_euler_PBC_2
    common/parameters/gamma,D,L,tau
    common/euler/C
    call init_genrand(12345)
    !parameters
    L=50.d0 
    D=15.d0 
    tau_vec=[1.d0,10.d0,100.d0]
    N=1000
    gamma=1.d0
    allocate(delta(1),MSD(1),correlacio(1,1),correlacio_j(1,1))
    allocate(val_t(1),vector_var(2*(1)),vector_PBC(2*(1)))

    !2
    m=25
    open(5,file="delta_t.dat")
    open(7,file="correlacio.dat")
    do j=1,3
        deallocate(delta,MSD,val_t,vector_var,vector_PBC,correlacio,correlacio_j)
        tau=tau_vec(j)
        dt=tau/20.d0
        t_max=tau*100.d0
        Np=int(t_max/dt)
        C=sqrt(2.d0*D*dt)/tau
        write(*,*)tau,Np
        allocate(delta(Np+1),MSD(Np+1),correlacio(1:Np+1,1:Np+1),&
            correlacio_j(1:Np+1,1:Np+1))
        allocate(val_t(Np+1),vector_var(2*(Np+1)),vector_PBC(2*(Np+1)))
        do k=1,Np+1
            delta(k)=0.d0
            do k2=1,Np+1
                correlacio(k,k2)=0.d0
            enddo
        enddo
        !average over particles
        do i=1,N
            do k=1,Np+1
                MSD(k)=0.d0
                do k2=1,Np+1
                    correlacio_j(k,k2)=0.d0
                enddo
            enddo
            !expected value of MSD for each particle
            do i_2=1,m
                x_0=[genrand_real3()*L-L/2.d0,0.d0]
                call recurrencia_v2_PBC(Np,2,1,x_0,0.d0,t_max,pas_euler_PBC_2,&
                    val_t,vector_var,vector_PBC,zero)
                do k=1,2*(Np+1)-1,2
                    r2=(vector_var(k)-x_0(1))**2.d0
                    MSD((k-1)/2+1)=MSD((k-1)/2+1)+r2/(m*1.d0)
                    do k2=1,k,2
                        correlacio_j((k-1)/2+1,(k2-1)/2+1)=&
                        correlacio_j((k-1)/2+1,(k2-1)/2+1)+&
                        vector_var(k+1)*vector_var(k2+1)/dble(m)
                    enddo
                enddo
            enddo
            do i_2=1,Np+1
                delta(i_2)=delta(i_2)+MSD(i_2)/dble(N)
                do k=1,i_2
                    correlacio(i_2,k)=correlacio(i_2,k)+correlacio_j(i_2,k)/dble(N)
                enddo
            enddo
        enddo
        write(5,*)
        write(5,*)
        write(5,*)"#Tau=",tau
        write(5,"(2(A16,3X))")"#t","Delta^2"
        write(7,*)
        write(7,*)
        write(7,*)"#Tau=",tau
        write(7,"(2(A16,3X))")"#t-t'","Correlacio"
        do k=1,Np+1
            write(5,"(2(e16.8,3X))")val_t(k),delta(k)
            do k2=1,k 
                write(7,"(2(e16.8,3X))")abs(val_t(k)-val_t(k2)),correlacio(k,k2)
            enddo
        enddo
    enddo
    close(5)
    close(7)
    


    write(*,*)"S'ha arribat al final del programa."
end program

double precision function zero(x)
    implicit none
    double precision :: x 
    zero=0.d0
    return
end function

subroutine boxmuller_2(sigma,mu,g1,g2)
!retorna 2 numeros aleatoris segons una dist. gaussiana
! amb valor esperat mu i desvacio estandard sigma
    implicit none
    integer :: comptar,i
    double precision :: g1,g2,sigma,mu
    double precision :: x,e1,e2,pi,y, genrand_real3
    pi=4.d0*atan(1.d0)
    e1=genrand_real3()
    e2=genrand_real3()
    g1=(-2.d0*log(1.d0-e1))**0.5d0*cos(2.d0*pi*e2)
    g2=(-2.d0*log(1.d0-e1))**0.5d0*sin(2.d0*pi*e2)
    g1=mu+sigma*g1
    g2=mu+sigma*g2
    return
end subroutine


subroutine pas_euler_PBC_2(x,y_in,y_in_PBC,nvar,h,y_out,y_out_PBC,F)
    implicit none
    integer :: nvar,i
    double precision :: y_in(nvar),y_out(nvar),h,x,f_calc(nvar),g,&
    gamma,D,tau,y_in_PBC(nvar),y_out_PBC(nvar),L,F,g2,C
    common/parameters/gamma,D,L,tau
    common/euler/C
    call boxmuller_2(1.d0,0.d0,g,g2)
    !position
    y_out(1)=y_in(1)+(y_in(2)+F(y_in(1)))*h/gamma
    y_out_PBC(1)=y_in_PBC(1)+(y_in_PBC(2)+F(y_in(1)))*h/gamma
    !random
    y_out(2)=y_in(2)*(1-h/tau)+C*g
    y_out_PBC(2)=y_in_PBC(2)*(1-h/tau)+C*g 
    !PBC conditions
    if (y_out_PBC(1).lt.(-L/2.d0)) then
        y_out_PBC(1)=L/2.d0+mod(y_out_PBC(1),L)
    else if (y_out_PBC(1).gt.L/2.d0) then
        y_out_PBC(1)=-L/2.d0+mod(y_out_PBC(1),L)
    endif
    return
end subroutine


subroutine recurrencia_v2_PBC(np,nv,grau,y0,xmin,xmax,metode,valx,valy,valy_PBC,F)
    implicit none
    integer :: np,nv,grau,i,j,k
    double precision :: xmax,xmin,y0(nv*grau),y_pas_in(nv*grau),&
    y_pas_out(nv),valx(np+1),valy(nv*(np+1)),h,xin(grau),y_pas_out_PBC(nv),&
    y_pas_in_PBC(nv*grau),valy_PBC(nv*(np+1)),F
    !np: numero de passos
    !nv: numero de variables
    !grau: grau del metode a utilitzar &
    !(ex: euler grau 1, Adams-Bashfoth grau 2)
    !y0 (vector amb els valors necessaris per incialitzar un pas)
    !y0=[a0,b0,...,a1,b1,...]
    !valy=[a1,b1,c1,...,a2,b2,c2,...,anp,bnp,cnp,...]
    !f--> subrutina que retorna les derivades del vector y 
    !valx=[x1,x2,...,xn]---temps
    h=(xmax-xmin)/np
    valx(1)=xmin
    valx((np+1))=xmax

    !generar una llista de numeros gaussians

    do j=1,grau
        do k=1,nv 
            y_pas_in(k+(j-1)*nv)=y0(nv*(grau-j)+k) 
            y_pas_in_PBC(k+(j-1)*nv)=y0(nv*(grau-j)+k)  
        enddo
    enddo
    do i=1,nv*grau
        valy(i)=y0(i)
        valy_PBC(i)=y0(i)
    enddo

    do i=2,(np+1)
        valx(i)=xmin+h*(i-1)
    enddo
    do i=1,grau
        xin(i)=valx(grau+1-i)
    enddo

    do i=grau+1,(np+1) 
        call metode(xin,y_pas_in,y_pas_in_PBC,nv,h,y_pas_out,y_pas_out_PBC,F)
        do j=1,nv
            valy((i-1)*nv+j)=y_pas_out(j)
            valy_PBC((i-1)*nv+j)=y_pas_out_PBC(j)
        enddo
        do j=1,grau
            xin(j)=valx(i+1-j)
            do k=1,nv
                y_pas_in(k+(j-1)*nv)=valy((i-j)*nv+k)
                y_pas_in_PBC(k+(j-1)*nv)=valy_PBC((i-j)*nv+k)
            enddo
        enddo
    enddo
    
    return
end subroutine

!------------------------------------------------------------------------

!mt19937ar

subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
!
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*&
               ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
!
      return
      end
!-----------------------------------------------------------------------
!     initialize by an array with array-length
!     init_key is the array for initializing keys
!     key_length is its length
!-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
!
      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)&
                +init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
!
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,0xffffffff]-interval
!-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
!
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
!
      y=mt(mti)
      mti=mti+1
!
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
!
      genrand_int32=y
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,0x7fffffff]-interval
!-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1]-real-interval
!-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on (0,1)-real-interval
!-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1) with 53-bit resolution
!-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
!-----------------------------------------------------------------------
!     initialize large number (over 32-bit constant number)
!-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!    TOPBIT_MASK = Z'80000000'
!    ALLBIT_MASK = Z'ffffffff'
!    UPPER_MASK  = Z'80000000'
!    LOWER_MASK  = Z'7fffffff'
!    MATRIX_A    = Z'9908b0df'
!    T1_MASK     = Z'9d2c5680'
!    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
