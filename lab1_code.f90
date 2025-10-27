!Author: Sergi Martínez Galindo
!1st Practice Non-Equilibrium Statistical Physics
program practica1
    implicit none
    integer :: N, i_min, i_max, Nb, i,index_min, index_max,Np,N_part,j,k,i_2,m
    double precision :: mean,s2,s,bs,tau,dt, gamma1, gamma2, gamma1_vec(1:5),L,&
    x_0(2),genrand_real3, r2, D(5), log_D(5), log_gamma_1(5), a,da,b,db,r_reg
    double precision, dimension(:), allocatable :: numbers_list,xhis,yhis,dyhis,&
    val_t, val_x, val_x_PBC,delta,MSD
    double precision, dimension(:,:,:,:), allocatable :: sim, sim_PBC
    common/parameters/gamma1,gamma2,L
    external :: pas_euler,pas_euler_PBC
    call init_genrand(12345)
    !parameters
    gamma1=30.d0 !capital gamma (variable)
    gamma2=1.d0 !damping coefficient
    N=10000 !gaussian distribution
    Nb=30 !# boxes histogram

    !2. Gaussian distribution: Box-Muller
    !number generation
    allocate(numbers_list(N),xhis(Nb), yhis(Nb), dyhis(Nb))
    call boxmuller(N,numbers_list,1.d0,0.d0)
    !parameters calculation
    call par_stats(N,numbers_list,mean,s2,s)

    write(*,*)"Mean=",mean,"Sigma=",s

    i_max=index_max(numbers_list, N)
    i_min=index_min(numbers_list, N)

    !histogram
    call histograma(N,numbers_list,numbers_list(i_min)-0.5d0,&
     numbers_list(i_max)+0.5d0,Nb, xhis,yhis,dyhis,bs)

    open(1, file="gaussian_results.dat")
    do i=1,Nb
        write(1,"(3(e16.8,3X))")xhis(i),yhis(i),dyhis(i)
    enddo
    close(1)

    !3. Trajectory of one particle for tau=100dt
    Np=100 
    allocate(val_t(Np+1),val_x(2*(Np+1)))
    dt=0.1d0
    tau=Np*dt

    !simulation of the trajectory without PBC
    call recurrencia_v2(Np,2,1,[0.d0,0.d0],0.d0,tau,pas_euler,val_t,val_x)

    open(2, file="euler.dat")
    do i=1,2*(Np+1)-1,2
        write(2,"(3(e16.8,3X))")val_t((i-1)/2+1),val_x(i),val_x(i+1)
    enddo
    close(2)

    !4. PBC implementation
    L=50.d0 !(box from -50 to 50)
    N_part=1000
    allocate(val_x_PBC(2*(Np+1)))
    open(3, file="apartat4.dat") !all trajectories
    open(4, file="distribution.dat") !initial and final distributions
    open(10, file="gif_coord.dat") !extra: gif creation (for personal visualization)
    write(4,"(5(A16,3X))")"#N,""x_0","y_0","x_end","y_end"
    do i=1,N_part
        x_0=[genrand_real3()*2*L,genrand_real3()*2*L]-[50.d0,50.d0]
        !trajectory simulation with and without PBC
        call recurrencia_v2_PBC(Np,2,1,x_0,0.d0,tau,pas_euler_PBC,val_t,val_x,&
            val_x_PBC)
        write(3,*)
        write(3,*)
        write(3,*)"#Particula numero:", i
        write(3,"(5(A16,3X))")"#t","x","y","x_PBC","y_PBC"
        write(4,"(5(e16.8,3X))")dble(i),val_x(1),val_x(2),val_x_PBC(2*Np+1),&
        val_x_PBC(2*Np+2)
        do j=1,2*(Np+1)-1,2
            write(3,"(5(e16.8,3X))")val_t((j-1)/2+1),val_x(j),val_x(j+1),&
            val_x_PBC(j),val_x_PBC(j+1)
            write(10,"(3(f7.3,3X))")val_x_PBC(j),val_x_PBC(j+1),val_t((j-1)/2+1)
        enddo
    enddo
    close(3)
    close(4)
    close(10)

    !5. MSD
    gamma1_vec=[0.1d0,1.d0,3.d0,10.d0,30.d0]
    allocate(delta(Np+1),MSD(Np+1))
    open(5,file="delta_t.dat")
    open(8,file="D_gamma.dat")
    write(8,"(2(A16,3X))")"#D","Gamma"
    m=100
    do j=1,5
        do k=1,Np+1
            delta(k)=0.d0
        enddo
        gamma1=gamma1_vec(j)
        !average over particles
        do i=1,N_part
            do k=1,Np+1
                MSD(k)=0.d0
            enddo
            !expected value of MSD for each particle
            do i_2=1,m
                x_0=[genrand_real3()*L,genrand_real3()*L]
                call recurrencia_v2_PBC(Np,2,1,x_0,0.d0,tau,pas_euler_PBC,&
                    val_t,val_x,val_x_PBC)
                do k=1,2*(Np+1)-1,2
                    r2=(val_x(k)-x_0(1))**2.d0+(val_x(k+1)-x_0(2))**2.d0
                    MSD((k-1)/2+1)=MSD((k-1)/2+1)+r2/(m*1.d0)
                enddo
            enddo
            do i_2=1,Np+1
                delta(i_2)=delta(i_2)+MSD(i_2)/dble(N_part)
            enddo
        enddo
        write(5,*)
        write(5,*)
        write(5,*)"#Gamma=",gamma1
        write(5,"(2(A16,3X))")"#t","Delta^2"
        do k=1,Np+1
            write(5,"(2(e16.8,3X))")val_t(k),delta(k)
        enddo
        !diffusivity
        D(j)=delta(Np+1)/(4.d0*val_t(Np+1))

        write(8,"(2(e16.8,3X))")D(j),gamma1
    enddo

    write(*,*)"Gamma loop done."

    close(5)
    close(8)

    !D vs gamma relationship: D=C*Gamma^A, logD=Alog(Gamma)+B, B=logC
    do i=1,5
        log_D(i)=log(D(i))
        log_gamma_1(i)=log(gamma1_vec(i))
    enddo

    call calc_reg(log_gamma_1,log_D,5,a,da,b,db,r_reg)

    write(*,*)"A=",a,"+-",da
    write(*,*)"B=",b,"+-",db
    write(*,*)"r=",r_reg

    write(*,*)"S'ha arribat al final del programa."
end program

subroutine calc_reg(val_x,val_y,Ndat,m,dm,b,db,r)
    implicit none
    integer*4 :: Ndat,i
    real*8 :: m,dm,b,db,r,val_x(1:Ndat),val_y(1:Ndat)
    real*8 :: sumx,sumx2,sumy,sumy2,sumxy
    real*8 :: varx,vary,varxy,dyreg

    !calcul sumatoris
    sumx=0.d0
    sumx2=0.d0
    sumy=0.d0
    sumy2=0.d0
    sumxy=0.d0
    do i=1,Ndat 
        sumx=sumx+val_x(i)
        sumx2=sumx2+val_x(i)**2.d0
        sumy=sumy+val_y(i)
        sumy2=sumy2+val_y(i)**2.d0
        sumxy=sumxy+val_x(i)*val_y(i)
    enddo

    !mitjanes
    sumx=sumx/Ndat
    sumx2=sumx2/Ndat
    sumy=sumy/Ndat
    sumy2=sumy2/Ndat
    sumxy=sumxy/Ndat

    !variances
    varx=sumx2-sumx**2.d0
    vary=sumy2-sumy**2.d0
    varxy=sumxy-sumx*sumy

    !regressio
    m=varxy/varx 
    b=sumy-m*sumx 
    r=varxy/(varx*vary)**0.5d0 

    !incertesses
    dyreg=((vary*(1.d0-r**2.d0)*Ndat)/(Ndat-2.d0))**0.5d0
    dm=dyreg/(Ndat*varx)**0.5d0
    db=dyreg*(sumx2/(Ndat*varx))**0.5d0

    return
end subroutine

subroutine recurrencia_v2(np,nv,grau,y0,xmin,xmax,metode,valx,valy)
    implicit none
    integer :: np,nv,grau,i,j,k
    double precision :: xmax,xmin,y0(nv*grau),y_pas_in(nv*grau),&
    y_pas_out(nv),valx(np+1),valy(nv*(np+1)),h,xin(grau)
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

    do j=1,grau
        do k=1,nv 
            y_pas_in(k+(j-1)*nv)=y0(nv*(grau-j)+k)   
        enddo
    enddo
    do i=1,nv*grau
        valy(i)=y0(i)
    enddo

    do i=2,(np+1)
        valx(i)=xmin+h*(i-1)
    enddo
    do i=1,grau
        xin(i)=valx(grau+1-i)
    enddo

    do i=grau+1,(np+1) 

        call metode(xin,y_pas_in,nv,h,y_pas_out)
        do j=1,nv
            valy((i-1)*nv+j)=y_pas_out(j)
        enddo
        do j=1,grau
            xin(j)=valx(i+1-j)
            do k=1,nv
                y_pas_in(k+(j-1)*nv)=valy((i-j)*nv+k)
            enddo
        enddo
    enddo
    
    return
end subroutine

subroutine pas_euler(x,y_in,nvar,h,y_out)
    implicit none
    integer :: nvar
    double precision :: y_in(nvar),y_out(nvar),h,x,f_calc(nvar),gx,gy,&
    gamma1,gamma2,g(2),L
    common/parameters/gamma1,gamma2,L
    call boxmuller_2(1.d0,0.d0,gx,gy)
    g=[gx,gy]
    y_out=y_in+g*sqrt(2.d0*gamma1*h)/gamma2
    return
end subroutine


subroutine recurrencia_v2_PBC(np,nv,grau,y0,xmin,xmax,metode,valx,valy,valy_PBC)
    implicit none
    integer :: np,nv,grau,i,j,k
    double precision :: xmax,xmin,y0(nv*grau),y_pas_in(nv*grau),&
    y_pas_out(nv),valx(np+1),valy(nv*(np+1)),h,xin(grau),y_pas_out_PBC(nv),&
    y_pas_in_PBC(nv*grau),valy_PBC(nv*(np+1))
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
        call metode(xin,y_pas_in,y_pas_in_PBC,nv,h,y_pas_out,y_pas_out_PBC)
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

subroutine pas_euler_PBC(x,y_in,y_in_PBC,nvar,h,y_out,y_out_PBC)
    implicit none
    integer :: nvar,i
    double precision :: y_in(nvar),y_out(nvar),h,x,f_calc(nvar),gx,gy,&
    gamma1,gamma2,g(2),y_in_PBC(nvar),y_out_PBC(nvar),L,L_vec(1:nvar)
    common/parameters/gamma1,gamma2,L
    call boxmuller_2(1.d0,0.d0,gx,gy)
    g=[gx,gy]
    y_out=y_in+g*sqrt(2.d0*gamma1*h)/gamma2
    do i=1,nvar
        L_vec(i)=L
    enddo
    y_out_PBC=y_in_PBC+g*sqrt(2.d0*gamma1*h)/gamma2
    !PBC conditions
    do i=1,nvar
        if (y_out_PBC(i).lt.(-L_vec(i))) then
            y_out_PBC(i)=L_vec(i)+mod(y_out_PBC(i),L_vec(i))
        else if (y_out_PBC(i).gt.L_vec(i)) then
            y_out_PBC(i)=-L_vec(i)+mod(y_out_PBC(i),L_vec(i))
        endif
    enddo
    return
end subroutine

integer function index_max(xvalues,Ndat)
    implicit none 
    integer :: i,Ndat 
    real*8 :: xvalues(1:Ndat),xmax
    xmax=-10**(6) 
    do i=1,Ndat 
        if (xvalues(i).gt.xmax) then
            index_max=i
            xmax=xvalues(i)
        endif
    enddo
    return
end function


integer function index_min(xvalues,Ndat)
    implicit none 
    integer :: i,Ndat 
    real*8 :: xvalues(1:Ndat),xmin
    xmin=10**(6) 
    do i=1,Ndat 
        if (xvalues(i).lt.xmin) then
            index_min=i
            xmin=xvalues(i)
        endif
    enddo
    return
end function


subroutine histograma(ndat,xdata,xa,xb,nbox,xhis&
    ,vhis,errhis,boxsize)
    implicit none
    integer :: ndat,nbox,i,nb
    double precision :: xdata(ndat), xhis(nbox),&
    errhis(nbox),vhis(nbox),boxsize,xa,xb&
    ,valorsx(nbox+1),h,nk(nbox)
    !amplada dels intervals (boxsize=h), equiespaiats
    h=(xb-xa)/nbox
    !creació llista limits de les caixes i valor centre
    valorsx(nbox+1)=xb
    do i=0,nbox-1
        valorsx(i+1)=xa+i*h
        xhis(i+1)=xa+i*h+h/2.d0
    enddo   
    !write(*,*) "On estan els punts?"
    !on esta cada punt? calcul de nk
    do i=1,nbox
        nk(i)=0
    enddo
    !write(*,*) "1r bucle"
    do i=1,ndat
        if (xdata(i).ge.xa.and.xdata(i).le.xb) then
            nb=int(((xdata(i)-xa)/h)+1)
            if (nb.eq.nbox+1) then
                nb=nb-1
            endif
            nk(nb)=nk(nb)+1
        endif

    enddo
    !write(*,*) "2n bucle"
    !calcul de vhis i de l'error
    !write(*,*) "calcul de vhis i de l'error"
    do i=1,nbox
        vhis(i)=nk(i)/(ndat*h)
        errhis(i)=(((nk(i)/ndat)*(1-nk(i)/ndat))**0.5d0)&
        /(h*(ndat)**0.5d0)
    enddo

    boxsize=h
    !write(*,*) "valor boxsize"
    return
end subroutine

subroutine par_stats(ndades,numeros,xmed,var,desv)
    implicit none
    double precision :: numeros(ndades),xmed,x2med,&
    var,desv,s,s2,sumx,sumx2
    integer :: i,ndades
    sumx=0.d0
    sumx2=0.d0
    do i=1,ndades
        sumx=sumx+numeros(i)
        sumx2=sumx2+numeros(i)**2.d0
    enddo
    xmed=sumx/ndades
    x2med=sumx2/ndades
    s2=(x2med-xmed**2.d0)*ndades/(ndades-1)
    s=s2**0.5d0
    var=s2 
    desv=s
    !respecte el valor mitja xmed
    !var=s2/ndades
    !desv=(var)**0.5
    return
end subroutine

subroutine boxmuller(ndades,numeros,sigma,mu)
!retorna numeros aleatoris segons una dist. gaussiana
! amb valor esperat mu i desvacio estandard sigma
    implicit none
    integer :: ndades,comptar,i
    double precision :: numeros(ndades),sigma,mu
    double precision :: x,e1,e2,pi,y, genrand_real3
    pi=4.d0*atan(1.d0)
    do i=1,ndades
        e1=genrand_real3()
        e2=genrand_real3()
        x=(-2.d0*log(1.d0-e1))**0.5d0*cos(2.d0*pi*e2)
        y=(-2.d0*log(1.d0-e1))**0.5d0*sin(2.d0*pi*e2)
        numeros(i)=mu+sigma*x
    enddo
    return
end subroutine

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
