      program virus

      implicit real*8 (a-h,o-z)
      parameter (ndim=300)  ! dimension del reticulo 300 x 300 = 90000 people, 1000 x 1000 = 1000000
      common/physics/temp,rhosigma,r0,rmax
      common/recuperacion/a,b,c

      dimension ns(ndim,ndim),ni(ndim,ndim) ! casillas de susceptibles
      dimension nr(ndim,ndim),nd(ndim,ndim)  ! casillas de recovered
      dimension nti(ndim,ndim)   ! tiempo de infeccion, el dia que se infecto
      dimension nlock(ndim,ndim)   ! cuarentena total
      
      ndays=250   ! number of days que quieres calcular 
      rmax=100     ! no hay que cambiarlo, cada persona infecta a los que hay alrededor  
      ntot=ndim**2  ! number of people 300 x 300 

c     parameters
      pinfeccion=0.15  ! probabilidad de infeccion de una persona, ajustado al caso de Espanya. Vas cambiadolo hasta que se ajuste a la curva. Todos los parametros tienen que cambiarse con logica.
      pdeath=0.5  !  probabilidad de muerte, quiere decir que 30% muere, sino hubiera gente en cuarentena total, sino moririan casi todos. Algunos se salvan, si se quedan aislados, por causalidad.      

c parameters recuperacion      
      tr=20   ! tiempo de recuperacion para recuperarse o morir, relacionado con el tau
      dtr=5.    ! error del tiempo de recuperacion, funcion D meseta, mirar articulo, 20+-4 dias, si el tiempo es menor de 20-4 la posibilidad de muerte es 0 y es maxima cuando es 24+4

      temp=20   !  temperatura en el modelo de Planck, n fotones van con la distribucion de Planck en unidades arbitrarias, constante de Boltzmann e-kT/t, para cada persona infectada la temperatura te dice cuantas veces va y viene, como la energia cinetica media, 20 vueltas por el vecindario, parecidos a Espanya
      rhosigma=0.09  ! la densidad por la seccion eficaz, probabilidad de interaccion, en Espana era 0.2. A lo mejor menos densidad o menos seccion eficaz en SA. Al haber menos interaccion hay menos infecciones. El producto de temp x rhosigma te dice cuantas interacciones va a tener una persona en promedio. Ahora es 20x0.12=2 interacciones. 
      r0=2.9   ! la distancia o radio de accion del movimiento de una persona, 5 casillas en promedio, valor medio una funcion exponencial que disminuye con esa constante. 
            
c     lockdown day
      lock=10000        ! no lockdown; 10000th day is the day of lockdown (starting day)
      rlock=1.5     !  desde ese dia 1.5 casillas a la izda y derecha, comparado con Espanya se ha subido, no sirve para el caso de Sudafrica, esta efectivo durante el primer dia. En Espanya es 5 casillas hacia izda y derecha.

c    total lockdown cells
      plocktot=0.   ! 77% estan en total lockdown en forma aleatoria, no se infecta. La gente en pueblos perdidos o encerrados en su casa. El 23% restante se pueden infectar. El 0.77 afecta desde el dia que se empieza el lockdown (en este caso el dia decimo). Hay que aumentarlo si incrementamos la poblacion.

c     8M      
      m8=0 ! dia que se hace una manifestacion. Ahora mismo no hay efecto 8M. 
      rhosigma8=0.2   ! probabilidad de infeccion del 8M. Ahora mismo no funciona porque m8=0.
      
c write parameters      
      open(10,file='virus.dat',status='unknown')
      write(10,100)ndim,pinfeccion,pdeath,tr,dtr
 100  format('#ndim,pinfection,pdeath,tr,dtr=',I5,12(1x, f12.5))
      write(10,110)temp,rhosigma,r0
 110  format('#temp,rhosigma,r0=',3(1x, f12.5))
      write(10,115)lock,rlock,plocktot
 115  format('#lock,rlock,plocktot=',i3,1x,f12.5,1x,f12.5)
      write(10,120)m8,rhosigma8
 120  format('#m8,rhosigma8=',1x, i3, 1x,f12.5)

c initial infections      
      ini=1
      nstot=ntot-ini
      nitot=ini
      nrtot=0
      ndtot=0
      
c     matrix init
      do i=1,ndim
         do j=1,ndim
            ns(i,j)=1
            ni(i,j)=0
            nr(i,j)=0
            nd(i,j)=0
            nti(i,j)=0
            nlock(i,j)=0            
         enddo
      enddo

c     initial condition infected in the midle point

      i1=ndim/2
      j1=ndim/2
      ns(i1,j1)=0
      ni(i1,j1)=1
      
      rs=rhosigma
c loop over pandemic days      

      do 10  nt = 1, ndays

         if(nt.eq.m8) then
            rhosigma=rhosigma8
c            rhosigma=rs
         else
            rhosigma=rs
         endif   
         if(nt.gt.lock) r0=rlock ! if lockdown modify range

c        total lock cells
         if(nt.eq.lock) then
         do 40 i=1,ndim
            do 41 j=1,ndim
               x=ran()
               if(x.lt.plocktot) nlock(i,j)=1
 41         continue
 40         continue
         endif  ! total lock cells
               
         
         nsdia=0
         nidia=0
         nrdia=0
         nddia=0

        
c loop over cells         
         do 20 i=1,ndim
            do 21 j=1,ndim
            
               if(ni(i,j).eq.1) then    ! INFECTED
                  
c     random number of interactions 
                  nint = ninteract()
c                  write(*,*) nint
c     random range 
                  r = radius()        ! migration / distancia que se puede mover cada persona                         

c     if lock is total no interactions
                  if(nlock(i,j).eq.1) nint=0
                                    
c     interaction with nint indivials in the square of width r

                  do 30 k=1,nint

      if(nlock(i,j).eq.1)write(*,*)'lock is interacting'
                     
c choose a random cell to interact with
                  x=ran()
                  y=ran()
                  i2=i-R+2.0*R*x
                  j2=j-R+2.0*R*y

c if cell is out of lattice move it to the border                   
                  if(i2.lt.1) i2=1
                  if(j2.lt.1) j2=1
                  if(i2.gt.ndim) i2=ndim
                  if(j2.gt.ndim) j2=ndim

c if cell is susceptible                  
                  if(ns(i2,j2).eq.1) then

c cell is locked
                     if(nlock(i2,j2).eq.1) then
                           pk=1.0
                        else
                           pk=ran()
                      endif     

                     if(pk.lt.pinfeccion) then
c INFECTED
                        ns(i2,j2)=0
                        ni(i2,j2)=1
                        nti(i2,j2)=nt
                        nsdia=nsdia-1
                        nstot=nstot-1
                        nidia=nidia+1
                        nitot=nitot+1                        
                     endif
                     
                  endif
 30               continue
                  
c   CHECK IF REMOVED AT NIGHT
                  t=nt-nti(i,j)
c   probabilidad de recuperacion
                  pr=1.0/(1.0+exp((tr-t)/dtr))
c                  write(*,*)pr1,pr

                  x = ran()
                  if(x.lt.1.e-4)x=ran()  ! repeat if too small

                  if( x.lt. pr) then
c                     write(*,*)'RECOVERED'
                     ni(i,j)=0
                     nr(i,j)=1
                     nrdia=nrdia+1
                     nrtot=nrtot+1
                     nidia=nidia-1
                     nitot=nitot-1
                     y=ran()
                     if(y.lt.pdeath) then
                        nd(i,j)=1
                        nr(i,j)=0
c                        write(*,*)'DEATH',I,J,NDTOT
                        nddia=nddia+1
                        ndtot=ndtot+1
                        nrdia=nrdia-1
                        nrtot=nrtot-1
                     endif
                  endif   

                  endif  !  end if infected

 21         continue
 20      continue

c  WRITE DAYLY AND TOTAL S I R D          ! aqui se escibren los susceptibles,... totales     
         write(10,200)nt,nsdia,nidia,nrdia,nddia,   ! number of deaths del dia
     *                   nstot,nitot,nrtot,ndtot
         write(*,*)nt, nddia,ndtot
 200     format(i4,9i7)

 10   continue

      stop
      end


c----------------------------      
      function precupera(t)
c---------------------------------      
            implicit real*8 (a-h,o-z)

            common/recuperacion/a,b,c
            precupera= a/(c+exp(-t/b))
            
            return
            end

c--------------------------------------      
            function ninteract()
c--------------------------------------      
            implicit real*8 (a-h,o-z)
            common/physics/temp,rhosigma,r0,rmax

c  genera un numero del 1 al 100            
            do 10  i=1,100000
            n=100*ran()
            pn= exp(-n/temp)/(exp(1.d0/temp)-1)

            y=ran()
            if(y.lt.pn) goto 20

 10      continue
         
 20      ninteract=n*rhosigma
            return
            end


      

c----------------------------
            function radius()
c----------------------------      
            implicit real*8 (a-h,o-z)
            common/physics/temp,rhosigma,r0,rmax

c  genera un numero del 1 a rmax            
            do 10  i=1,100000
            r=rmax*ran()
            pr= exp(-r/r0)/r0
            y=ran()
            if(y.lt.pr) goto 20

 10      continue

 20      radius=r
            return
            end
