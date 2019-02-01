module he3fci
implicit none
integer,parameter :: prec = selected_real_kind(P=15)
real(prec),parameter :: ZERO = 0._prec

private
public HE3 

contains

! calculates nonadditive 3-body energies from the He3 FCI fit (fitted to 253 points)
! reads in distances in bohr
      ! IMPLICIT NONE

!      REAL*8 e2k , e3 , r1 , r2 , r3 , unit1 , unit2
!      INTEGER i , n

!      DATA e2k/315774.65D0/
 
!      unit=1.d6  ! microhartrees
      ! unit1 = 1.D0
                  ! ! hartrees
      ! unit2 = e2k ! kelvins
 
      ! WRITE (*,'(a)') '     R1        R2        R3        unit1'//      &
                     ! &'           unit2            '
      ! READ (*,*) n
      ! DO i = 1 , n
         ! READ (*,*) r1 , r2 , r3
         ! CALL HE3(r1,r2,r3,e3)
         ! WRITE (*,'(3f10.6,3g14.6)') r1 , r2 , r3 , e3*unit1 , e3*unit2
      ! ENDDO
      ! END

!---------------------------------------------------------------------
      SUBROUTINE HE3(R1,R2,R3,E3)
! R1,R2,R3: in Angstrom; E3: in Cassandra atomic units
      IMPLICIT NONE

      REAL*8 COStheta1 , COStheta2 , COStheta3 , E3 , R1 , R2 , R3 , R1_bohr, R2_bohr, R3_bohr,    &
           & THEta1 , THEta2 , THEta3

      COMMON /THETAS/ THEta1 , THEta2 , THEta3 , COStheta1 , COStheta2 ,&
                    & COStheta3
 
      COStheta1 = 0.5D0*(R1*R1+R2*R2-R3*R3)/(R1*R2)
      COStheta2 = 0.5D0*(R1*R1+R3*R3-R2*R2)/(R1*R3)
      COStheta3 = 0.5D0*(R2*R2+R3*R3-R1*R1)/(R2*R3)
      IF ( COStheta1.GT.1.0D0 ) COStheta1 = 1.0D0
      IF ( COStheta1.LT.-1.0D0 ) COStheta1 = -1.0D0
      IF ( COStheta2.GT.1.0D0 ) COStheta2 = 1.0D0
      IF ( COStheta2.LT.-1.0D0 ) COStheta2 = -1.0D0
      IF ( COStheta3.GT.1.0D0 ) COStheta3 = 1.0D0
      IF ( COStheta3.LT.-1.0D0 ) COStheta3 = -1.0D0
      THEta1 = DACOS(COStheta1)
      THEta2 = DACOS(COStheta2)
      THEta3 = DACOS(COStheta3)
	  
	  !Convert Cassandra units (Angstrom) to bohr before passing into DOFCI
	  R1_bohr = R1 / 0.529177249_prec
	  R2_bohr = R2 / 0.529177249_prec
	  R3_bohr = R3 / 0.529177249_prec
 
      CALL DOFCI(R1_bohr,R2_bohr,R3_bohr,E3)
	  
	  E3 = E3 * 315775.13_prec * 0.8314472_prec !Return E3 in Cassandra atomic units
      
      END SUBROUTINE HE3

!---------------------------------------------------------------------
      SUBROUTINE DOFCI(R1,R2,R3,Corr)
! R1,R2,R3: in bohr; E3: in hartree
      IMPLICIT NONE
      
      REAL*8 a , a333 , a344 , a355 , a434 , a443 , a445 , a454 , a535 ,&
           & a544 , a553 , a555 , beta , Corr , COStheta1 , COStheta2 , &
           & COStheta3 , d066 , d068 , d077
      REAL*8 d086 , d333 , d338 , d344 , d347 , d355 , d374 ,      &
           & d383 , d434 , d437 , d443 , d445 , d446 , d454 , d464 ,    &
           & d473 , d535 , d544 , d553
      REAL*8 d555 , d606 , d608 , d644 , d660 , d680 , d707 , d734 ,    &
           & d743 , d770 , d806 , d833 , d860 , pn11 , pn12 , pn13 ,&
           & pn21 , pn22 , pn23
      REAL*8 pn31 , pn32 , pn33 , R1 , r13 , r14 , r15 , r16 , r17 ,    &
           & r18 , R2 , r23 , r24 , r25 , r26 , r27 , r28 , R3 , r33 ,  &
           & r34
      REAL*8 r35 , r36 , r37 , r38 , term , THEta1 , THEta2 , THEta3 ,  &
           & w068 , w077 , w086 , w111 , w112 , w113 , w12 , w121 ,     &
           & w122 , w13 , w131 , w211
      REAL*8 w212 , w221 , w222 , w23 , w311 , w338 , w347 , w374 ,     &
           & w383 , w437 , w446 , w464 , w473 , w608 , w644 , w680 ,    &
           & w707 , w734 , w743 , w770
      REAL*8 w806 , w833 , w860 , z111 , z112 , z113 , z122 , z222
      INTEGER init , j , k , kmax , n1 , n2 , n3
!*** End of declarations inserted by SPAG
      COMMON /THETAS/ THEta1 , THEta2 , THEta3 , COStheta1 , COStheta2 ,&
                    & COStheta3
      DIMENSION a(100) , beta(100)
      SAVE init
      DATA init/1/
      DATA a/100*0.D0/
      DATA beta/100*0.D0/
 
! exact values (Thakkar JCP 75, 4496)
      DATA z111/0.49311D0/
      DATA z112/0.92372D0/
      DATA z113/4.1241D0/
      DATA z122/1.7377D0/
      DATA z222/3.2839D0/

      IF ( init.EQ.1 ) THEN
         CALL INITAB('E3.dat',beta,46,a,41)
         init = 0
      ENDIF
 
      kmax = 4 ! nb=35+6
 
      r13 = 1.D0/R1**3.0D0
      r23 = 1.D0/R2**3.0D0
      r33 = 1.D0/R3**3.0D0
      w111 = 3.0D0*r13*r23*r33*                                         &
           & (1.0D0+3.0D0*COStheta1*COStheta2*COStheta3)
 
      r24 = 1.D0/R2**4.0D0
      r34 = 1.D0/R3**4.0D0
      w112 = 0.1875D0*r13*r24*r34*                                      &
           & ((9.0D0*COStheta3-25.0D0*DCOS(3.0D0*THEta3))               &
           & +6.0D0*DCOS(THEta1-THEta2)*(3.0D0+5.0D0*DCOS(2.0D0*THEta3))&
           & )
 
      r14 = 1.D0/R1**4.D0
      w121 = 0.1875D0*r23*r14*r34*                                      &
           & ((9.0D0*COStheta2-25.0D0*DCOS(3.0D0*THEta2))               &
           & +6.0D0*DCOS(THEta3-THEta1)*(3.0D0+5.0D0*DCOS(2.0D0*THEta2))&
           & )
 
      w211 = 0.1875D0*r33*r24*r14*                                      &
           & ((9.0D0*COStheta1-25.0D0*DCOS(3.0D0*THEta1))               &
           & +6.0D0*DCOS(THEta2-THEta3)*(3.0D0+5.0D0*DCOS(2.0D0*THEta1))&
           & )
 
      r25 = 1.D0/R2**5.0D0
      r35 = 1.D0/R3**5.0D0
      w113 = 0.15625D0*r13*r25*r35*(9.0D0+8.0D0*DCOS(2.0D0*THEta3)      &
           & -49.0D0*DCOS(4.0D0*THEta3)+6.0D0*DCOS(THEta1-THEta2)       &
           & *(9.0D0*COStheta3+7.0D0*DCOS(3.0D0*THEta3)))
 
      r15 = 1.D0/R1**5.0D0
      w131 = 0.15625D0*r33*r25*r15*(9.0D0+8.0D0*DCOS(2.0D0*THEta1)      &
           & -49.0D0*DCOS(4.0D0*THEta1)+6.0D0*DCOS(THEta3-THEta2)       &
           & *(9.0D0*COStheta1+7.0D0*DCOS(3.0D0*THEta1)))
 
      w311 = 0.15625D0*r23*r15*r35*(9.0D0+8.0D0*DCOS(2.0D0*THEta2)      &
           & -49.0D0*DCOS(4.0D0*THEta2)+6.0D0*DCOS(THEta1-THEta3)       &
           & *(9.0D0*COStheta2+7.0D0*DCOS(3.0D0*THEta2)))
 
      w122 = 0.234375D0*r14*r25*r34*                                    &
           & (3.0D0*(COStheta2+5.0D0*DCOS(3.0D0*THEta2))                &
           & +20.0D0*DCOS(THEta1-THEta3)                                &
           & *(1.0D0-3.0D0*DCOS(2.0D0*THEta2))                          &
           & +70.0D0*DCOS(2.0D0*(THEta1-THEta3))*COStheta2)
 
      w212 = 0.234375D0*r24*r15*r34*                                    &
           & (3.0D0*(COStheta3+5.0D0*DCOS(3.0D0*THEta3))                &
           & +20.0D0*DCOS(THEta1-THEta2)                                &
           & *(1.0D0-3.0D0*DCOS(2.0D0*THEta3))                          &
           & +70.0D0*DCOS(2.0D0*(THEta1-THEta2))*COStheta3)
 
      w221 = 0.234375D0*r24*r35*r14*                                    &
           & (3.0D0*(COStheta1+5.0D0*DCOS(3.0D0*THEta1))                &
           & +20.0D0*DCOS(THEta3-THEta2)                                &
           & *(1.0D0-3.0D0*DCOS(2.0D0*THEta1))                          &
           & +70.0D0*DCOS(2.0D0*(THEta3-THEta2))*COStheta1)
 
      w222 = 0.1171875D0*r15*r25*r35*                                   &
           & (-27.0D0+220.0D0*COStheta1*COStheta2*COStheta3+            &
           & 490.0D0*DCOS(2.0D0*THEta1)*DCOS(2.0D0*THEta2)              &
           & *DCOS(2.0D0*THEta3)                                        &
           & +175.0D0*(DCOS(2.0D0*(THEta1-THEta2))+DCOS                 &
           & (2.0D0*(THEta2-THEta3))+DCOS(2.0D0*(THEta3-THEta1))))
 
      a333 = w111*z111
      a344 = w112*z112
      a434 = w121*z112
      a443 = w211*z112
      a355 = w113*z113
      a553 = w131*z113
      a535 = w311*z113
      a454 = w122*z122
      a544 = w212*z122
      a445 = w221*z122
      a555 = w222*z222
 
      Corr = 0.D0
!--- exponential part
 
      k = 0
      j = 0
      DO n1 = 0 , kmax
         DO n2 = n1 , kmax
            DO n3 = n2 , kmax
               k = k + 1
               j = j + 1
               pn11 = P(n1,COStheta1)
               pn12 = P(n1,COStheta2)
               pn13 = P(n1,COStheta3)
               pn21 = P(n2,COStheta1)
               pn22 = P(n2,COStheta2)
               pn23 = P(n2,COStheta3)
               pn31 = P(n3,COStheta1)
               pn32 = P(n3,COStheta2)
               pn33 = P(n3,COStheta3)
!
               term = DEXP(-beta(k)*(R1+R2+R3))                         &
                    & *(pn11*pn22*pn33+pn11*pn23*pn32+pn12*pn21*pn33+   &
                    & pn12*pn23*pn31+pn13*pn21*pn32+pn13*pn22*pn31)
               Corr = Corr + a(j)*term
            ENDDO
         ENDDO
      ENDDO
 
!--- exponential part: end
 
 
!--- 3.order part
 
! damping functions .....
      k = k + 1
      d333 = D3(3,3,3,beta(k),R1,R2,R3)
      k = k + 1
      d344 = D3(3,4,4,beta(k),R1,R2,R3)
      d434 = D3(4,3,4,beta(k),R1,R2,R3)
      d443 = D3(4,4,3,beta(k),R1,R2,R3)
      k = k + 1
      d355 = D3(3,5,5,beta(k),R1,R2,R3)
      d553 = D3(5,5,3,beta(k),R1,R2,R3)
      d535 = D3(5,3,5,beta(k),R1,R2,R3)
      k = k + 1
      d454 = D3(4,5,4,beta(k),R1,R2,R3)
      d544 = D3(5,4,4,beta(k),R1,R2,R3)
      d445 = D3(4,4,5,beta(k),R1,R2,R3)
      k = k + 1
      d555 = D3(5,5,5,beta(k),R1,R2,R3)
 
! damped 3.order
      term = d333*a333
      term = term + d344*a344 + d434*a434 + d443*a443
      term = term + d355*a355 + d553*a553 + d535*a535
      term = term + d454*a454 + d544*a544 + d445*a445
      term = term + d555*a555
 
      Corr = Corr + term
 
!--- 3.order part: end
 
 
 
!--------- 4th order part
 
      r16 = 1.D0/R1**6.D0
      r26 = 1.D0/R2**6.D0
      r36 = 1.D0/R3**6.D0
 
      w13 = 9.0D0*r16*r36*(1.0D0+COStheta2**2)
 
      w23 = 9.0D0*r36*r26*(1.0D0+COStheta3**2)
 
      w12 = 9.0D0*r26*r16*(1.0D0+COStheta1**2)
 
      r23 = 1.D0/R2**3
      r33 = 1.D0/R3**3
      r13 = 1.D0/R1**3
 
      r17 = r16/R1
      r34 = 1.D0/R3**4
      w734 = 0.03125D0*r17*r23*r34*(-144.0D0*COStheta2+                 &
           & 36.0D0*DCOS(THEta1+THEta3)+216.0D0*DCOS(THEta1-THEta3)     &
           & -120.0D0*DCOS(3.0D0*THEta2)                                &
           & -720.0D0*DCOS(THEta2-2.0D0*THEta3)                         &
           & -72.0D0*DCOS(THEta2-2.0D0*THEta1))
 
      r24 = 1.D0/R2**4
      w743 = 0.03125D0*r17*r24*r33*(-144.0D0*COStheta1+                 &
           & 36.0D0*DCOS(THEta2+THEta3)+216.0D0*DCOS(THEta2-THEta3)     &
           & -120.0D0*DCOS(3.0D0*THEta1)                                &
           & -720.0D0*DCOS(THEta1-2.0D0*THEta3)                         &
           & -72.0D0*DCOS(THEta1-2.0D0*THEta2))
 
      r27 = 1.D0/R2**7
      w374 = 0.03125D0*r13*r27*r34*(-144.0D0*COStheta3+                 &
           & 36.0D0*DCOS(THEta1+THEta2)+216.0D0*DCOS(THEta1-THEta2)     &
           & -120.0D0*DCOS(3.0D0*THEta3)                                &
           & -720.0D0*DCOS(THEta3-2.0D0*THEta2)                         &
           & -72.0D0*DCOS(THEta3-2.0D0*THEta1))
 
      r14 = 1.D0/R1**4
      w473 = 0.03125D0*r14*r27*r33*(-144.0D0*COStheta1+                 &
           & 36.0D0*DCOS(THEta2+THEta3)+216.0D0*DCOS(THEta2-THEta3)     &
           & -120.0D0*DCOS(3.0D0*THEta1)                                &
           & -720.0D0*DCOS(THEta1-2.0D0*THEta2)                         &
           & -72.0D0*DCOS(THEta1-2.0D0*THEta3))
 
      r37 = 1.D0/R3**7
      w347 = 0.03125D0*r13*r24*r37*(-144.0D0*COStheta3+                 &
           & 36.0D0*DCOS(THEta1+THEta2)+216.0D0*DCOS(THEta1-THEta2)     &
           & -120.0D0*DCOS(3.0D0*THEta3)                                &
           & -720.0D0*DCOS(THEta3-2.0D0*THEta1)                         &
           & -72.0D0*DCOS(THEta3-2.0D0*THEta2))
 
      w437 = 0.03125D0*r14*r23*r37*(-144.0D0*COStheta2+                 &
           & 36.0D0*DCOS(THEta1+THEta3)+216.0D0*DCOS(THEta1-THEta3)     &
           & -120.0D0*DCOS(3.0D0*THEta2)                                &
           & -720.0D0*DCOS(THEta2-2.0D0*THEta1)                         &
           & -72.0D0*DCOS(THEta2-2.0D0*THEta3))
 
      w644 = 0.03125D0*r16*r24*r34*(-111.0D0*COStheta3-                 &
           & 750.0D0*DCOS(3.0D0*THEta3)+180.0D0*DCOS(THEta1+THEta2)     &
           & +108.0D0*DCOS(THEta1-THEta2)                               &
           & -90.0D0*DCOS(THEta3-2.0D0*THEta1)                          &
           & -90.0D0*DCOS(THEta3-2.0D0*THEta2))
 
      w464 = 0.03125D0*r14*r26*r34*(-111.0D0*COStheta2-                 &
           & 750.0D0*DCOS(3.0D0*THEta2)+180.0D0*DCOS(THEta1+THEta3)     &
           & +108.0D0*DCOS(THEta1-THEta3)                               &
           & -90.0D0*DCOS(THEta2-2.0D0*THEta1)                          &
           & -90.0D0*DCOS(THEta2-2.0D0*THEta3))
 
      w446 = 0.03125D0*r14*r24*r36*(-111.0D0*COStheta1-                 &
           & 750.0D0*DCOS(3.0D0*THEta1)+180.0D0*DCOS(THEta2+THEta3)     &
           & +108.0D0*DCOS(THEta2-THEta3)                               &
           & -90.0D0*DCOS(THEta1-2.0D0*THEta2)                          &
           & -90.0D0*DCOS(THEta1-2.0D0*THEta3))
 
      r18 = r17/R1
      w833 = -0.5D0*r18*r23*r33*(9.0D0*(DCOS(2.0D0*THEta1)+DCOS(2.0D0*  &
           & THEta2))+54.0D0*DCOS(2.0D0*THEta3))
 
      r28 = r27/R2
      w383 = -0.5D0*r13*r28*r33*(9.0D0*(DCOS(2.0D0*THEta1)+DCOS(2.0D0*  &
           & THEta3))+54.0D0*DCOS(2.0D0*THEta2))
 
      r38 = r37/R3
      w338 = -0.5D0*r13*r23*r38*(9.0D0*(DCOS(2.0D0*THEta2)+DCOS(2.0D0*  &
           & THEta3))+54.0D0*DCOS(2.0D0*THEta1))
 
      w707 = 0.03125D0*r17*r37*(1485.0D0*COStheta2+                     &
           & 384.0D0*DCOS(3.0D0*THEta2))
 
      w770 = 0.03125D0*r17*r27*(1485.0D0*COStheta1+                     &
           & 384.0D0*DCOS(3.0D0*THEta1))
 
      w077 = 0.03125D0*r27*r37*(1485.0D0*COStheta3+                     &
           & 384.0D0*DCOS(3.0D0*THEta3))
 
      w707 = -0.5D0*w707
      w770 = -0.5D0*w770
      w077 = -0.5D0*w077
 
 
 
      w860 = 0.25D0*r18*r26*(369.0D0+288.0D0*COStheta1**2)
 
      w680 = 0.25D0*r16*r28*(369.0D0+288.0D0*COStheta1**2)
 
      w806 = 0.25D0*r18*r36*(369.0D0+288.0D0*COStheta2**2)
 
      w608 = 0.25D0*r16*r38*(369.0D0+288.0D0*COStheta2**2)
 
      w086 = 0.25D0*r28*r36*(369.0D0+288.0D0*COStheta3**2)
 
      w068 = 0.25D0*r26*r38*(369.0D0+288.0D0*COStheta3**2)
 
 
! damping functions....
      k = k + 1
      d660 = D3(6,6,0,beta(k),R1,R2,R3)
      d606 = D3(6,0,6,beta(k),R1,R2,R3)
      d066 = D3(0,6,6,beta(k),R1,R2,R3)
      k = k + 1
      d743 = D3(7,4,3,beta(k),R1,R2,R3)
      d734 = D3(7,3,4,beta(k),R1,R2,R3)
      d374 = D3(3,7,4,beta(k),R1,R2,R3)
      d473 = D3(4,7,3,beta(k),R1,R2,R3)
      d347 = D3(3,4,7,beta(k),R1,R2,R3)
      d437 = D3(4,3,7,beta(k),R1,R2,R3)
      k = k + 1
      d644 = D3(6,4,4,beta(k),R1,R2,R3)
      d464 = D3(4,6,4,beta(k),R1,R2,R3)
      d446 = D3(4,4,6,beta(k),R1,R2,R3)
      k = k + 1
      d833 = D3(8,3,3,beta(k),R1,R2,R3)
      d383 = D3(3,8,3,beta(k),R1,R2,R3)
      d338 = D3(3,3,8,beta(k),R1,R2,R3)
      k = k + 1
      d770 = D3(7,7,0,beta(k),R1,R2,R3)
      d707 = D3(7,0,7,beta(k),R1,R2,R3)
      d077 = D3(0,7,7,beta(k),R1,R2,R3)
      k = k + 1
      d860 = D3(8,6,0,beta(k),R1,R2,R3)
      d680 = D3(6,8,0,beta(k),R1,R2,R3)
      d806 = D3(8,0,6,beta(k),R1,R2,R3)
      d608 = D3(6,0,8,beta(k),R1,R2,R3)
      d086 = D3(0,8,6,beta(k),R1,R2,R3)
      d068 = D3(0,6,8,beta(k),R1,R2,R3)
 
! damped:
      j = j + 1
      Corr = Corr + a(j)*(d660*w12+d606*w13+d066*w23)
      j = j + 1
      Corr = Corr + a(j)                                                &
           & *(d743*w743+d734*w734+d374*w374+d473*w473+d347*w347+       &
           & d437*w437)
      j = j + 1
      Corr = Corr + a(j)*(d644*w644+d464*w464+d446*w446)
      j = j + 1
      Corr = Corr + a(j)*(d833*w833+d383*w383+d338*w338)
      j = j + 1
      Corr = Corr + a(j)*(d770*w770+d707*w707+d077*w077)
      j = j + 1
      Corr = Corr + a(j)                                                &
           & *(d860*w860+d680*w680+d806*w806+d608*w608+d086*w086+       &
           & d068*w068)
 
!--------- 4th order part: end
 
      END SUBROUTINE DOFCI

!---------------------------------------------------------------------
      FUNCTION P(L,Theta)
      IMPLICIT NONE

      INTEGER L
      REAL*8 P , Theta

      IF ( L.EQ.0 ) THEN
         P = 1.0D0
      ELSEIF ( L.EQ.1 ) THEN
         P = Theta
      ELSEIF ( L.EQ.2 ) THEN
         P = 0.5D0*(3.0D0*Theta*Theta-1.0D0)
      ELSEIF ( L.EQ.3 ) THEN
         P = 0.5D0*(5.D0*Theta*Theta*Theta-3.D0*Theta)
      ELSEIF ( L.EQ.4 ) THEN
         P = 0.125D0*(35.D0*Theta**4-30.D0*Theta*Theta+3.D0)
      ELSE
         WRITE (*,*) 'l too large in P'
         STOP
      ENDIF
      END FUNCTION P

!---------------------------------------------------------------------
      FUNCTION D3(N1,N2,N3,Beta,R1,R2,R3)
!
!     calculate the damping factor
!
      IMPLICIT NONE

      REAL*8 Beta , br , D3 , dx , dy , dz , R1 , R2 , R3 , sum , term
      INTEGER i , N1 , N2 , N3 , ncn

      IF ( N1.EQ.0 ) THEN
         dx = 1.0D0
      ELSE
         br = Beta*R1
         sum = 1.0D0
         term = 1.0D0
         ncn = N1
         DO i = 1 , ncn
            term = term*br/i
            sum = sum + term
         ENDDO
         dx = 1.0D0 - DEXP(-br)*sum
      ENDIF
      IF ( N2.EQ.0 ) THEN
         dy = 1.0D0
      ELSE
         br = Beta*R2
         sum = 1.0D0
         term = 1.0D0
         ncn = N2
         DO i = 1 , ncn
            term = term*br/i
            sum = sum + term
         ENDDO
         dy = 1.0D0 - DEXP(-br)*sum
      ENDIF
      IF ( N3.EQ.0 ) THEN
         dz = 1.0D0
      ELSE
         br = Beta*R3
         sum = 1.0D0
         term = 1.0D0
         ncn = N3
         DO i = 1 , ncn
            term = term*br/i
            sum = sum + term
         ENDDO
         dz = 1.0D0 - DEXP(-br)*sum
      ENDIF
      D3 = dx*dy*dz
!     write(6,*) n,beta,r,d
      END FUNCTION D3
!---------------------------------------------------------------------
      SUBROUTINE INITAB(Fname,A1,L1,A2,L2)
      IMPLICIT NONE

      INTEGER i , ii , L1 , L2
      CHARACTER*6 Fname
      REAL*8 A1(*) , A2(*)
      OPEN (7,FILE=Fname)
      READ (7,*,ERR=100)
      READ (7,*,ERR=100)
      READ (7,*,ERR=100)
      DO i = 1 , L1
         READ (7,*,ERR=100) ii , A1(i)
         IF ( i.NE.ii ) GOTO 100
      ENDDO
      READ (7,*,ERR=100)
      DO i = 1 , L2
         READ (7,*,ERR=100) ii , A2(i)
         IF ( i.NE.ii ) GOTO 100
      ENDDO
      CLOSE (7)
      RETURN
 100  WRITE (*,'(a,a)') 'Error reading file ' , Fname
      STOP
      END SUBROUTINE INITAB
!---------------------------------------------------------------------


end module he3fci