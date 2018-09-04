c calculates nonadditive 3-body energies from the He3 FCI fit (fitted to 253 points)
c reads in distances in bohr
      implicit real*8 (a-h,o-z)
      data E2K /315774.65d0/

c      unit=1.d6  ! microhartrees
      unit1=1.d0  ! hartrees
      unit2=E2K   ! kelvins

      write (*,'(a)')'     R1        R2        R3        unit1'
     >// '           unit2            '
      read (*,*) n
      do 10 i=1,n
      read (*,*) R1,R2,R3
      call He3 (R1,R2,R3, E3)
      write (*,'(3f10.6,3g14.6)') R1,R2,R3, E3*unit1,E3*unit2
 10   continue
      end
c---------------------------------------------------------------------
      subroutine He3 (R1,R2,R3, E3)
c R1,R2,R3: in bohr; E3: in hartree
      implicit real*8 (a-h,o-z)
      common /thetas/ theta1,theta2,theta3,costheta1,costheta2,costheta3

      costheta1 = 0.5d0*(r1*r1+r2*r2-r3*r3)/(r1*r2)
      costheta2 = 0.5d0*(r1*r1+r3*r3-r2*r2)/(r1*r3)
      costheta3 = 0.5d0*(r2*r2+r3*r3-r1*r1)/(r2*r3)
      If(costheta1.gt.1.0d0) costheta1 = 1.0d0
      If(costheta1.lt.-1.0d0) costheta1 = -1.0d0
      If(costheta2.gt.1.0d0) costheta2 = 1.0d0
      If(costheta2.lt.-1.0d0) costheta2 = -1.0d0
      If(costheta3.gt.1.0d0) costheta3 = 1.0d0
      If(costheta3.lt.-1.0d0) costheta3 = -1.0d0
      theta1 = dacos(costheta1)
      theta2 = dacos(costheta2)
      theta3 = dacos(costheta3)

      call doFCI (R1,R2,R3, E3)

      end
c---------------------------------------------------------------------
      subroutine doFCI (R1,R2,R3, corr)
      implicit real*8 (a-h,o-z)
      common /thetas/ theta1,theta2,theta3,costheta1,costheta2,costheta3
      dimension A(100),beta(100)
      save init
      data init /1/
      data A /100*0.d0/
      data beta /100*0.d0/

c exact values (Thakkar JCP 75, 4496) 
       data z111 /0.49311d0/
       data z112 /0.92372d0/
       data z113 /4.1241d0/
       data z122 /1.7377d0/
       data z222 /3.2839d0/

      if (init.eq.1) then
         call initab ('E3.dat',beta,46,A,41)
         init=0
      end if

      Kmax=4   ! nb=35+6
                      
      r13 = 1.d0/R1**3.0d0
      r23 = 1.d0/R2**3.0d0
      r33 = 1.d0/R3**3.0d0
      w111 = 3.0d0*r13*r23*r33*
     &       (1.0d0+3.0d0*
     &        costheta1*costheta2*costheta3)

       r24 = 1.d0/R2**4.0d0
       r34 = 1.d0/R3**4.0d0
       w112 = 0.1875d0*r13*r24*r34*(
     &        (9.0d0*costheta3-25.0d0*dcos(3.0d0*theta3))+
     &        6.0d0*dcos(theta1-theta2)*(3.0d0+5.0d0*
     &        dcos(2.0d0*theta3)))

       r14 = 1.d0/R1**4.d0
       w121 = 0.1875d0*r23*r14*r34*(
     &        (9.0d0*costheta2-25.0d0*dcos(3.0d0*theta2))+
     &        6.0d0*dcos(theta3-theta1)*(3.0d0+5.0d0*
     &        dcos(2.0d0*theta2)))

       w211 = 0.1875d0*r33*r24*r14*(
     &        (9.0d0*costheta1-25.0d0*dcos(3.0d0*theta1))+
     &        6.0d0*dcos(theta2-theta3)*(3.0d0+5.0d0*
     &        dcos(2.0d0*theta1)))

       r25 = 1.d0/R2**5.0d0
       r35 = 1.d0/R3**5.0d0
       w113 = 0.15625d0*r13*r25*r35
     $      * (9.0d0 + 8.0d0*dcos(2.0d0*theta3) 
     $         - 49.0d0*dcos(4.0d0*theta3)
     $         + 6.0d0*dcos(theta1-theta2)*
     $           (9.0d0*costheta3 + 7.0d0*dcos(3.0d0*theta3)))

       r15 = 1.d0/R1**5.0d0
       w131 = 0.15625d0*r33*r25*r15
     $      * (9.0d0 + 8.0d0*dcos(2.0d0*theta1) 
     $         - 49.0d0*dcos(4.0d0*theta1)
     $         + 6.0d0*dcos(theta3-theta2)*
     $           (9.0d0*costheta1 + 7.0d0*dcos(3.0d0*theta1)))

       w311 = 0.15625d0*r23*r15*r35
     $      * (9.0d0 + 8.0d0*dcos(2.0d0*theta2) 
     $         - 49.0d0*dcos(4.0d0*theta2)
     $         + 6.0d0*dcos(theta1-theta3)*
     $           (9.0d0*costheta2 + 7.0d0*dcos(3.0d0*theta2)))

       w122 = 0.234375d0*r14*r25*r34*(
     &        3.0d0*(costheta2+5.0d0*dcos(3.0d0*theta2)) +
     &    20.0d0*dcos(theta1-theta3)*(1.0d0-3.0d0
     &        *dcos(2.0d0*theta2))+
     &        70.0d0*dcos(2.0d0*(theta1-theta3))*costheta2)

       w212 = 0.234375d0*r24*r15*r34*(
     &        3.0d0*(costheta3+5.0d0*dcos(3.0d0*theta3)) +
     &    20.0d0*dcos(theta1-theta2)*(1.0d0-3.0d0
     &        *dcos(2.0d0*theta3))+
     &        70.0d0*dcos(2.0d0*(theta1-theta2))*costheta3)

       w221 = 0.234375d0*r24*r35*r14*(
     &        3.0d0*(costheta1+5.0d0*dcos(3.0d0*theta1)) +
     &    20.0d0*dcos(theta3-theta2)*(1.0d0-3.0d0
     &        *dcos(2.0d0*theta1))+
     &        70.0d0*dcos(2.0d0*(theta3-theta2))*costheta1)

       w222 = 0.1171875d0*r15*r25*r35*(
     &  -27.0d0+220.0d0*costheta1*costheta2*costheta3+
     & 490.0d0*dcos(2.0d0*theta1)*dcos(2.0d0*theta2)
     &       *dcos(2.0d0*theta3)
     &       +175.0d0*(dcos(2.0d0*(theta1-theta2))+
     &                 dcos(2.0d0*(theta2-theta3))+
     &                 dcos(2.0d0*(theta3-theta1))))
      
      a333 =  w111*z111  
      a344 =  w112*z112 
      a434 =  w121*z112 
      a443 =  w211*z112 
      a355 =  w113*z113 
      a553 =  w131*z113 
      a535 =  w311*z113 
      a454 =  w122*z122 
      a544 =  w212*z122 
      a445 =  w221*z122 
      a555 =  w222*z222 

      corr=0.d0
c--- exponential part

      k = 0
      j = 0
         do n1 = 0,Kmax
            do n2 =n1,Kmax
               do n3=n2,Kmax
               k = k + 1
               j = j + 1
               pn11 = p(n1,costheta1)
               pn12 = p(n1,costheta2)
               pn13 = p(n1,costheta3)
               pn21 = p(n2,costheta1)
               pn22 = p(n2,costheta2)
               pn23 = p(n2,costheta3)
               pn31 = p(n3,costheta1)
               pn32 = p(n3,costheta2)
               pn33 = p(n3,costheta3)
c
        term = dexp(-beta(k)*(R1+R2+R3))*(
     $          pn11*pn22*pn33
     $         +pn11*pn23*pn32
     $         +pn12*pn21*pn33
     $         +pn12*pn23*pn31
     $         +pn13*pn21*pn32
     $         +pn13*pn22*pn31)
        corr=corr +A(j)*term
               end do
            end do
         end do

c--- exponential part: end


c--- 3.order part

c damping functions .....
         k=k+1
        d333 = d3(3,3,3,beta(k),R1,R2,R3)
        k=k+1
        d344 = d3(3,4,4,beta(k),R1,R2,R3)
        d434 = d3(4,3,4,beta(k),R1,R2,R3)
        d443 = d3(4,4,3,beta(k),R1,R2,R3)
        k=k+1
        d355 = d3(3,5,5,beta(k),R1,R2,R3)
        d553 = d3(5,5,3,beta(k),R1,R2,R3)
        d535 = d3(5,3,5,beta(k),R1,R2,R3)
        k=k+1
        d454 = d3(4,5,4,beta(k),R1,R2,R3)
        d544 = d3(5,4,4,beta(k),R1,R2,R3)
        d445 = d3(4,4,5,beta(k),R1,R2,R3)
        k=k+1
        d555 = d3(5,5,5,beta(k),R1,R2,R3)

c damped 3.order
      term = d333*a333
      term = term + d344*a344 + d434*a434 + d443*a443
      term = term + d355*a355 + d553*a553 + d535*a535
      term = term + d454*a454 + d544*a544 + d445*a445
      term = term + d555*a555

      corr=corr +term

c--- 3.order part: end



c--------- 4th order part

      r16 = 1.d0/R1**6.d0
      r26 = 1.d0/R2**6.d0
      r36 = 1.d0/R3**6.d0

      w13 = 9.0d0*r16*r36*
     $ (1.0d0+costheta2**2)

      w23 = 9.0d0*r36*r26*
     $ (1.0d0+costheta3**2)

      w12 = 9.0d0*r26*r16*
     $ (1.0d0+costheta1**2)

      r23 = 1.d0/R2**3
      r33 = 1.d0/R3**3
      r13 = 1.d0/R1**3

      r17 = r16/R1
      r34 = 1.d0/R3**4
      w734 = 0.03125d0*r17*r23*r34
     $   *(-144.0d0*costheta2
     $     +36.0d0*dcos(theta1+theta3)
     $     +216.0d0*dcos(theta1-theta3)
     $     -120.0d0*dcos(3.0d0*theta2)
     $     -720.0d0*dcos(theta2-2.0d0*theta3) 
     $     -72.0d0*dcos(theta2-2.0d0*theta1))  

      r24 = 1.d0/R2**4
      w743 = 0.03125d0*r17*r24*r33
     $   *(-144.0d0*costheta1
     $     +36.0d0*dcos(theta2+theta3)
     $     +216.0d0*dcos(theta2-theta3)
     $     -120.0d0*dcos(3.0d0*theta1)
     $     -720.0d0*dcos(theta1-2.0d0*theta3) 
     $     -72.0d0*dcos(theta1-2.0d0*theta2))  

      r27 = 1.d0/R2**7
      w374 = 0.03125d0*r13*r27*r34
     $   *(-144.0d0*costheta3
     $     +36.0d0*dcos(theta1+theta2)
     $     +216.0d0*dcos(theta1-theta2)
     $     -120.0d0*dcos(3.0d0*theta3)
     $     -720.0d0*dcos(theta3-2.0d0*theta2) 
     $     -72.0d0*dcos(theta3-2.0d0*theta1))  

      r14 = 1.d0/R1**4
      w473 = 0.03125d0*r14*r27*r33
     $   *(-144.0d0*costheta1
     $     +36.0d0*dcos(theta2+theta3)
     $     +216.0d0*dcos(theta2-theta3)
     $     -120.0d0*dcos(3.0d0*theta1)
     $     -720.0d0*dcos(theta1-2.0d0*theta2) 
     $     -72.0d0*dcos(theta1-2.0d0*theta3))  

      r37 = 1.d0/R3**7
      w347 = 0.03125d0*r13*r24*r37
     $   *(-144.0d0*costheta3
     $     +36.0d0*dcos(theta1+theta2)
     $     +216.0d0*dcos(theta1-theta2)
     $     -120.0d0*dcos(3.0d0*theta3)
     $     -720.0d0*dcos(theta3-2.0d0*theta1) 
     $     -72.0d0*dcos(theta3-2.0d0*theta2))  

      w437 = 0.03125d0*r14*r23*r37
     $   *(-144.0d0*costheta2
     $     +36.0d0*dcos(theta1+theta3)
     $     +216.0d0*dcos(theta1-theta3)
     $     -120.0d0*dcos(3.0d0*theta2)
     $     -720.0d0*dcos(theta2-2.0d0*theta1) 
     $     -72.0d0*dcos(theta2-2.0d0*theta3))  

      w644 = 0.03125d0*r16*r24*r34
     $     * ( -111.0d0*costheta3
     $         -750.0d0*dcos(3.0d0*theta3)
     $         +180.0d0*dcos(theta1+theta2)
     $         +108.0d0*dcos(theta1-theta2)
     $         -90.0d0*dcos(theta3-2.0d0*theta1)
     $         -90.0d0*dcos(theta3-2.0d0*theta2))

      w464 = 0.03125d0*r14*r26*r34
     $     * ( -111.0d0*costheta2
     $         -750.0d0*dcos(3.0d0*theta2)
     $         +180.0d0*dcos(theta1+theta3)
     $         +108.0d0*dcos(theta1-theta3)
     $         -90.0d0*dcos(theta2-2.0d0*theta1)
     $         -90.0d0*dcos(theta2-2.0d0*theta3))

      w446 = 0.03125d0*r14*r24*r36
     $     * ( -111.0d0*costheta1
     $         -750.0d0*dcos(3.0d0*theta1)
     $         +180.0d0*dcos(theta2+theta3)
     $         +108.0d0*dcos(theta2-theta3)
     $         -90.0d0*dcos(theta1-2.0d0*theta2)
     $         -90.0d0*dcos(theta1-2.0d0*theta3))

      r18 = r17/R1
      w833 = -0.5d0*r18*r23*r33
     $     * ( 9.0d0*(dcos(2.0d0*theta1)+dcos(2.0d0*theta2))
     $       + 54.0d0*dcos(2.0d0*theta3))

      r28 = r27/R2
      w383 = -0.5d0*r13*r28*r33
     $     * ( 9.0d0*(dcos(2.0d0*theta1)+dcos(2.0d0*theta3))
     $       + 54.0d0*dcos(2.0d0*theta2))

      r38 = r37/R3 
      w338 = -0.5d0*r13*r23*r38
     $     * ( 9.0d0*(dcos(2.0d0*theta2)+dcos(2.0d0*theta3))
     $       + 54.0d0*dcos(2.0d0*theta1))

      w707 = 0.03125d0*r17*r37*(1485.0d0*costheta2
     $                       +384.0d0*dcos(3.0d0*theta2))

      w770 = 0.03125d0*r17*r27*(1485.0d0*costheta1
     $                       +384.0d0*dcos(3.0d0*theta1))

      w077 = 0.03125d0*r27*r37*(1485.0d0*costheta3
     $                       +384.0d0*dcos(3.0d0*theta3))

      w707= -0.5d0*w707
      w770= -0.5d0*w770
      w077= -0.5d0*w077



      w860 = 0.25d0*r18*r26*(369.0d0
     $                       +288.0d0*costheta1**2)

      w680 = 0.25d0*r16*r28*(369.0d0
     $                       +288.0d0*costheta1**2)

      w806 = 0.25d0*r18*r36*(369.0d0
     $                       +288.0d0*costheta2**2)

      w608 = 0.25d0*r16*r38*(369.0d0
     $                       +288.0d0*costheta2**2)

      w086 = 0.25d0*r28*r36*(369.0d0
     $                       +288.0d0*costheta3**2)

      w068 = 0.25d0*r26*r38*(369.0d0
     $                       +288.0d0*costheta3**2)


c damping functions....
      k = k + 1
      d660=d3(6,6,0,beta(k),R1,R2,R3)
      d606=d3(6,0,6,beta(k),R1,R2,R3)
      d066=d3(0,6,6,beta(k),R1,R2,R3)
      k = k + 1
      d743=d3(7,4,3,beta(k),R1,R2,R3)
      d734=d3(7,3,4,beta(k),R1,R2,R3)
      d374=d3(3,7,4,beta(k),R1,R2,R3)
      d473=d3(4,7,3,beta(k),R1,R2,R3)
      d347=d3(3,4,7,beta(k),R1,R2,R3)
      d437=d3(4,3,7,beta(k),R1,R2,R3)
      k = k + 1
      d644=d3(6,4,4,beta(k),R1,R2,R3)
      d464=d3(4,6,4,beta(k),R1,R2,R3)
      d446=d3(4,4,6,beta(k),R1,R2,R3)
      k = k + 1
      d833=d3(8,3,3,beta(k),R1,R2,R3)
      d383=d3(3,8,3,beta(k),R1,R2,R3)
      d338=d3(3,3,8,beta(k),R1,R2,R3)
      k = k + 1
      d770=d3(7,7,0,beta(k),R1,R2,R3)
      d707=d3(7,0,7,beta(k),R1,R2,R3)
      d077=d3(0,7,7,beta(k),R1,R2,R3)
      k = k + 1
      d860=d3(8,6,0,beta(k),R1,R2,R3)
      d680=d3(6,8,0,beta(k),R1,R2,R3)
      d806=d3(8,0,6,beta(k),R1,R2,R3)
      d608=d3(6,0,8,beta(k),R1,R2,R3)
      d086=d3(0,8,6,beta(k),R1,R2,R3)
      d068=d3(0,6,8,beta(k),R1,R2,R3)

c damped:
      j= j + 1
      corr=corr +A(j)* (d660*w12+d606*w13+d066*w23)
      j= j + 1
      corr=corr +A(j)* (d743*w743+d734*w734+d374*w374
     >          +d473*w473+d347*w347+d437*w437)
      j= j + 1
      corr=corr +A(j)*  (d644*w644+d464*w464+d446*w446)
      j= j + 1
      corr=corr +A(j)* (d833*w833+d383*w383+d338*w338)
      j= j + 1
      corr=corr +A(j)* (d770*w770+d707*w707+d077*w077)
      j= j + 1
      corr=corr +A(j)* (d860*w860+d680*w680+d806*w806
     >          +d608*w608+d086*w086+d068*w068)

c--------- 4th order part: end

      end
c---------------------------------------------------------------------
      function p(l,theta)
      implicit real*8 (a-h,o-z)
      if (l.eq.0) then
        p=1.0d0
      else if (l.eq.1) then
        p = theta
      else if (l.eq.2) then
        p=0.5d0*(3.0d0*theta*theta - 1.0d0)
      else if (l.eq.3) then
         p=0.5d0*(5.d0*theta*theta*theta - 3.d0*theta)
      else if (l.eq.4) then
        p=0.125d0*(35.d0*theta**4 - 30.d0*theta*theta + 3.d0)
      else
         write (*,*) 'l too large in P'
         stop
      endif
      end
c---------------------------------------------------------------------
      function d3(n1,n2,n3,beta,r1,r2,r3)
c
c     calculate the damping factor
c
      implicit real*8 (a-h,o-z)
      if (n1.eq.0) go to 10 
      br=beta*r1  
      sum=1.0d0
      term=1.0d0
      ncn=n1 
      do i=1,ncn
        term=term*br/i  
        sum=sum+term
      enddo
      dx=1.0d0 - dexp(-br)*sum
      go to 20 
10    dx = 1.0d0 
20    continue 
      if (n2.eq.0) go to 30 
      br=beta*r2  
      sum=1.0d0
      term=1.0d0
      ncn=n2 
      do i=1,ncn
        term=term*br/i  
        sum=sum+term
      enddo
      dy=1.0d0 - dexp(-br)*sum
      go to 40
30    dy = 1.0d0 
40    continue
      if (n3.eq.0) go to 50 
      br=beta*r3  
      sum=1.0d0
      term=1.0d0
      ncn=n3 
      do i=1,ncn
        term=term*br/i  
        sum=sum+term
      enddo
      dz=1.0d0 - dexp(-br)*sum
      go to 60 
50    dz = 1.0d0 
60    continue 
      d3 = dx*dy*dz
c     write(6,*) n,beta,r,d
      return
      end
c---------------------------------------------------------------------
      subroutine initab (fname,a1,l1,a2,l2)
      character*6 fname
      real*8 a1(*),a2(*)
      open (7,file=fname)
      read (7,*,err=999)
      read (7,*,err=999)
      read (7,*,err=999)
      do 10 i=1,l1
      read (7,*,err=999) ii,a1(i)
      if (i.ne.ii) goto 999
 10   continue
      read (7,*,err=999)
      do 20 i=1,l2
      read (7,*,err=999) ii,a2(i)
      if (i.ne.ii) goto 999
 20   continue
      close (7)
      return
 999  write (*,'(a,a)') 'Error reading file ',fname
      stop
      end
c---------------------------------------------------------------------
