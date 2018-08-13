module potentials
implicit none
integer,parameter :: prec = selected_real_kind(P=15)
real(prec),parameter :: ZERO = 0._prec

private
public V_BO,V_AD,V_REL,V_QED,V,V_diff
public mu2V_NA,mu2Wr,mu2Wv,diff_mu2Wv
public sigma_BO,sigma_AD,sigma_REL,sigma_QED,sigma

!-------------------------------------------------------------------------------

integer,parameter :: M_BO = 3
integer,parameter :: I0_BO = -1, I1_BO = 2
integer,parameter :: N0_BO =  6, N1_BO = 16

real(prec),parameter :: a_BO(M_BO) = &
     (/&
     2.19076900924970_prec,&
     3.91558362097102_prec,&
     8.87314362466591_prec &
     /)
real(prec),parameter :: P_BO(I0_BO:I1_BO,M_BO) = reshape(source=&
     (/&
      2.69150966652494e1_prec,&
      1.74015968871153e1_prec,&
     -3.26566208189471e0_prec,&
      1.46700433938961e-1_prec,&
     -2.62457650950449e2_prec,&
      4.57300508401525e2_prec,&
      6.66624947516144e2_prec,&
     -1.01541164656195e2_prec,&
      2.39542554285199e2_prec,&
      6.73223355754347e2_prec,&
      8.62651807552327e2_prec,&
      7.86212667513655e2_prec &
     /),&
     shape=(/I1_BO-I0_BO+1,M_BO/))

real(prec),parameter :: zeta_BO = &
     4.45565354034820_prec
real(prec),parameter :: C_BO(N0_BO:N1_BO) = &
     (/&
       1.460977837725_prec,&
       ZERO,&
      14.11785737_prec,&
       ZERO,&
     183.691075_prec,&
     -76.72571_prec,&
       3.372e3_prec,&
   -3808.3254_prec,&
       8.534e4_prec,&
      -1.707e5_prec,&
       2.86e6_prec &
     /)

!-------------------------------------------------------------------------------

integer,parameter :: M_AD = 3
integer,parameter :: I0_AD = 0, I1_AD = 2
integer,parameter :: N0_AD = 6, N1_AD = 10

real(prec),parameter :: a_AD(M_AD) = &
     (/&
     2.75067513975782_prec,&
     5.36085266623527_prec,&
     7.13731208196370_prec &
     /)
real(prec),parameter :: P_AD(I0_AD:I1_AD,M_AD) = reshape(source=&
     (/&
      1.60034916931265e-2_prec,&
     -5.17915211148225e-3_prec,&
      3.20021447656502e-3_prec,&
     -5.02322894521138e-1_prec,&
     -5.28178298205425e-2_prec,&
      2.93999790052146e-1_prec,&
     -7.80222735840106e-2_prec,&
      1.83457620987126e0_prec,&
      9.46016993696249e-1_prec &
     /),&
     shape=(/I1_AD-I0_AD+1,M_AD/))

real(prec),parameter :: zeta_AD = &
     4.80971750762740_prec
real(prec),parameter :: C_AD(N0_AD:N1_AD) = &
     (/&
     1.1445e-3_prec,&
     ZERO,&
     6.519e-3_prec,&
     ZERO,&
     6.68e-2_prec &
     /)

!-------------------------------------------------------------------------------

integer,parameter :: M_REL = 3
integer,parameter :: I0_REL = 0, I1_REL = 2
integer,parameter :: N0_REL = 4, N1_REL = 8

real(prec),parameter :: a_REL(M_REL) = &
     (/&
     2.55442787544336_prec,&
     4.26822206589066_prec,&
     6.31801545306225_prec &
     /)
real(prec),parameter :: P_REL(I0_REL:I1_REL,M_REL) = reshape(source=&
     (/&
      1.16259940182396e-2_prec,&
     -7.00805799572921e-3_prec,&
      7.83673396684808e-4_prec,&
     -5.73711010584383e-1_prec,&
      3.97167157865319e-1_prec,&
     -1.10854901515699e-1_prec,&
      6.22980694403664e-1_prec,&
      4.87089036456869e-1_prec,&
      6.81223617553767e-1_prec &
     /),&
     shape=(/I1_REL-I0_REL+1,M_REL/))

real(prec),parameter :: zeta_REL = &
     5.00212547573053_prec
real(prec),parameter :: C_REL(N0_REL:N1_REL) = &
     (/&
     -3.5322e-5_prec,&
      ZERO,&
     -3.434e-4_prec,&
      ZERO,&
     -3.55798027534830e-3_prec &
     /)

!-------------------------------------------------------------------------------

integer,parameter :: M_QED = 2
integer,parameter :: I0_QED = 0, I1_QED = 2
integer,parameter :: N0_QED = 3, N1_QED = 6

real(prec),parameter :: a_QED(M_QED) = &
     (/&
     2.71385934118196_prec,&
     4.85539104988208_prec &
     /)
real(prec),parameter :: P_QED(I0_QED:I1_QED,M_QED) = reshape(source=&
     (/&
      1.34952154947189e-3_prec,&
     -1.02362538302740e-3_prec,&
      3.56862380020579e-4_prec,&
     -1.96866479382540e-4_prec,&
     -9.61841561868819e-4_prec,&
     -1.08014679908054e-3_prec &
     /),&
     shape=(/I1_QED-I0_QED+1,M_QED/))

real(prec),parameter :: zeta_QED = &
     5.43735514285525_prec
real(prec),parameter :: C_QED(N0_QED:N1_QED) = &
     (/&
     5.772353e-7_prec,&
     ZERO,&
     1.377841e-6_prec,&
     7.61187886972970e-5_prec &
     /)

!-------------------------------------------------------------------------------

integer,parameter :: M_NA = 3
integer,parameter :: I0_NA = 0, I1_NA = 2
integer,parameter :: N0_NA = 6, N1_NA = 8

real(prec),parameter :: a_NA(M_NA) = &
     (/&
     2.53831898631499_prec,&
     2.67481454650236_prec,&
     4.93815115302358_prec &
     /)
real(prec),parameter :: P_NA(I0_NA:I1_NA,M_NA) = reshape(source=&
     (/&
      1.47701736144261e2_prec,&
     -2.58061851134832e2_prec,&
     -1.11536126798668e2_prec,&
      3.46349558392539e2_prec,&
     -2.83069102118783e2_prec,&
      3.33770327404411e2_prec,&
      6.10470730307975e2_prec,&
     -8.69955973665686e1_prec,&
     -7.00699585484731e2_prec &
     /),&
     shape=(/I1_NA-I0_NA+1,M_NA/))

real(prec),parameter :: zeta_NA = &
     1.01792095068557_prec
real(prec),parameter :: C_NA(N0_NA:N1_NA) = &
     (/&
     3.68e0_prec,&
     ZERO,&
     1.11e2_prec &
     /)

!-------------------------------------------------------------------------------

integer,parameter :: M_Wr = 2
integer,parameter :: I0_Wr = 0, I1_Wr = 2
integer,parameter :: N0_Wr = 8, N1_Wr = 8
real(prec),parameter :: a_Wr(M_Wr) = &
     (/&
     2.54103162649609_prec,&
     4.88260628664777_prec &
     /)
real(prec),parameter :: P_Wr(I0_Wr:I1_Wr,M_Wr) = reshape(source=&
     (/&
      5.14574051998990e1_prec,&
     -3.77595327620630e1_prec,&
     -2.30240580489152e-6_prec,&
     -2.37370348095298e3_prec,&
      2.69532593707433e3_prec,&
     -1.37210463348773e3_prec &
     /),&
     shape=(/I1_Wr-I0_Wr+1,M_Wr/))

real(prec),parameter :: zeta_Wr = &
     3.37428444329544_prec
real(prec),parameter :: C_Wr(N0_Wr:N1_Wr) = &
     (/&
     -2.49e1_prec &
     /)

!-------------------------------------------------------------------------------

integer,parameter :: M_Wv = 2
integer,parameter :: I0_Wv = 0, I1_Wv = 3
integer,parameter :: N0_Wv = 8, N1_Wv = 8
real(prec),parameter :: a_Wv(M_Wv) = &
     (/&
     2.22718808355432_prec,&
     5.19023925799186_prec &
     /)
real(prec),parameter :: P_Wv(I0_Wv:I1_Wv,M_Wv) = reshape(source=&
     (/&
      8.01042582721457e1_prec,&
     -7.97402957505628e1_prec,&
      3.16348688660220e1_prec,&
     -2.80107566421106e0_prec,&
     -1.81464982364927e3_prec,&
      2.33555265113313e3_prec,&
     -1.91611678345999e3_prec,&
      5.32423263513008e2_prec &
     /),&
     shape=(/I1_Wv-I0_Wv+1,M_Wv/))

real(prec),parameter :: zeta_Wv = &
     3.18684428430590_prec
real(prec),parameter :: C_Wv(N0_Wv:N1_Wv) = &
     (/&
     1.49e2_prec &
     /)

!-------------------------------------------------------------------------------

! first  value - exponent
! second value - coefficient

integer,parameter :: nS_BO = 3
real(prec),parameter :: aS_BO(2,0:nS_BO) = reshape(source=&
     (/&
     0.80740_prec , 1.2195e-6_prec ,&
     1.2603_prec  , 2.6076e-5_prec ,&
     0.13365_prec ,-1.6989e-7_prec ,&
     0.057123_prec,-5.9399e-8_prec  &
     /),&
     shape=(/2,nS_BO+1/))

integer,parameter :: nS_AD = 3
real(prec),parameter :: aS_AD(2,0:nS_AD) = reshape(source=&
     (/&
     0.30940_prec , 7.5317e-12_prec,&
     1.8485_prec  , 5.8720e-6_prec ,&
     0.27382_prec , 4.8110e-8_prec ,&
     0.047007_prec, 1.0593e-10_prec &
     /),&
     shape=(/2,nS_AD+1/))

integer,parameter :: nS_REL = 4
real(prec),parameter :: aS_REL(2,0:nS_REL) = reshape(source=&
     (/&
     0.20535_prec , 3.0628e-12_prec,&
     1.6829_prec  , 2.3492e-7_prec ,&
     0.42205_prec ,-1.2725e-6_prec ,&
     0.39672_prec , 1.2107e-6_prec ,&
     0.064498_prec, 3.8763e-10_prec &
     /),&
     shape=(/2,nS_REL+1/))

integer,parameter :: nS_QED = 3
real(prec),parameter :: aS_QED(2,0:nS_QED) = reshape(source=&
     (/&
     0.80740_prec , 3.3214e-7_prec ,&
     0.50144_prec ,-1.0022e-7_prec ,&
     0.16251_prec ,-6.5718e-8_prec ,&
     0.062390_prec,-2.1856e-8_prec  &
     /),&
     shape=(/2,nS_QED+1/))

!-------------------------------------------------------------------------------

contains

function V_BO(R)
implicit none
real(prec) :: V_BO
real(prec),intent(in) :: R
V_BO = 0._prec
V_BO = V_BO + formShort(R,M_BO,I0_BO,I1_BO,a_BO,P_BO)
V_BO = V_BO + formLong(R,0,N0_BO,N1_BO,zeta_BO,C_BO)
end function V_BO

function V_AD(R)
implicit none
real(prec) :: V_AD
real(prec),intent(in) :: R
V_AD = 0._prec
V_AD = V_AD + formShort(R,M_AD,I0_AD,I1_AD,a_AD,P_AD)
V_AD = V_AD + formLong(R,0,N0_AD,N1_AD,zeta_AD,C_AD)
end function V_AD

function V_REL(R)
implicit none
real(prec) :: V_REL
real(prec),intent(in) :: R
V_REL = 0._prec
V_REL = V_REL + formShort(R,M_REL,I0_REL,I1_REL,a_REL,P_REL)
V_REL = V_REL + formLong(R,0,N0_REL,N1_REL,zeta_REL,C_REL)
end function V_REL

function V_QED(R)
implicit none
real(prec) :: V_QED
real(prec),intent(in) :: R
V_QED = 0._prec
V_QED = V_QED + formShort(R,M_QED,I0_QED,I1_QED,a_QED,P_QED)
V_QED = V_QED + formLong(R,0,N0_QED,N1_QED,zeta_QED,C_QED)
end function V_QED

function V(R_temp,ret)
implicit none
real(prec) :: V
real(prec),intent(in) :: R_temp
logical,intent(in),optional :: ret
logical :: do_ret
real(prec) :: elms(4)
real(prec) :: R
R = R_temp / 0.529177249_prec
do_ret = .false.
if(present(ret)) do_ret = ret
elms = 0._prec
elms(1) = elms(1) + formShort(R,M_BO,I0_BO,I1_BO,a_BO,P_BO)
elms(2) = elms(2) + formShort(R,M_AD,I0_AD,I1_AD,a_AD,P_AD)
elms(3) = elms(3) + formShort(R,M_REL,I0_REL,I1_REL,a_REL,P_REL)
elms(4) = elms(4) + formShort(R,M_QED,I0_QED,I1_QED,a_QED,P_QED)
if(do_ret) then
   elms(1) = elms(1) + formLong(R,1,N0_BO,N1_BO,zeta_BO,C_BO)
   elms(2) = elms(2) + formLong(R,0,N0_AD,N1_AD,zeta_AD,C_AD)
   elms(3) = elms(3) + formLong(R,2,N0_REL,N1_REL,zeta_REL,C_REL)
   elms(4) = elms(4) + formLong(R,2,N0_QED,N1_QED,zeta_QED,C_QED)
else
   elms(1) = elms(1) + formLong(R,0,N0_BO,N1_BO,zeta_BO,C_BO)
   elms(2) = elms(2) + formLong(R,0,N0_AD,N1_AD,zeta_AD,C_AD)
   elms(3) = elms(3) + formLong(R,0,N0_REL,N1_REL,zeta_REL,C_REL)
   elms(4) = elms(4) + formLong(R,0,N0_QED,N1_QED,zeta_QED,C_QED)
endif
V = sum(elms) * 315775.13_prec
end function V

! I added this function, meant to be the energy derivative.
! It cannot handle a value of ret which does not equal 0,
! but I allowed it to be passed for consistency and future
! changes.
function V_diff(R_temp, ret)
implicit none
real(prec) :: V_diff
real(prec), intent(in) :: R_temp
logical, intent(in), optional :: ret
real(prec) :: elms(4)
real(prec) :: R
elms= 0._prec
R = R_temp / 0.529177249_prec
! Assumes ret is false for now
elms(1) = elms(1) + diffShort(R,M_BO,I0_BO,I1_BO,a_BO,P_BO) 
elms(2) = elms(2) + diffShort(R,M_AD,I0_AD,I1_AD,a_AD,P_AD)
elms(3) = elms(3) + diffShort(R,M_REL,I0_REL,I1_REL,a_REL,P_REL)
elms(4) = elms(4) + diffShort(R,M_QED,I0_QED,I1_QED,a_QED,P_QED)
elms(1) = elms(1) + diffLong(R,0,N0_BO,N1_BO,zeta_BO,C_BO)
elms(2) = elms(2) + diffLong(R,0,N0_AD,N1_AD,zeta_AD,C_AD)
elms(3) = elms(3) + diffLong(R,0,N0_REL,N1_REL,zeta_REL,C_REL)
elms(4) = elms(4) + diffLong(R,0,N0_QED,N1_QED,zeta_QED,C_QED)
V_diff = sum(elms)
V_diff = V_diff * 315775.13_prec * 1.8897261254535_prec
end function V_diff

function mu2V_NA(R)
implicit none
real(prec) :: mu2V_NA
real(prec),intent(in) :: R
mu2V_NA = 0._prec
mu2V_NA = mu2V_NA + formShort(R,M_NA,I0_NA,I1_NA,a_NA,P_NA)
mu2V_NA = mu2V_NA + formLong(R,0,N0_NA,N1_NA,zeta_NA,C_NA)
end function mu2V_NA

function mu2Wr(R)
implicit none
real(prec) :: mu2Wr
real(prec),intent(in) :: R
mu2Wr = 0._prec
mu2Wr = mu2Wr + formShort(R,M_Wr,I0_Wr,I1_Wr,a_Wr,P_Wr)
mu2Wr = mu2Wr + formLong(R,0,N0_Wr,N1_Wr,zeta_Wr,C_Wr)
end function mu2Wr

function mu2Wv(R)
implicit none
real(prec) :: mu2Wv
real(prec),intent(in) :: R
mu2Wv = 0._prec
mu2Wv = mu2Wv + formShort(R,M_Wv,I0_Wv,I1_Wv,a_Wv,P_Wv)
mu2Wv = mu2Wv + formLong(R,0,N0_Wv,N1_Wv,zeta_Wv,C_Wv)
end function mu2Wv

function diff_mu2Wv(R)
implicit none
real(prec) :: diff_mu2Wv
real(prec),intent(in) :: R
diff_mu2Wv = 0._prec
diff_mu2Wv = diff_mu2Wv + diffShort(R,M_Wv,I0_Wv,I1_Wv,a_Wv,P_Wv)
diff_mu2Wv = diff_mu2Wv + diffLong(R,0,N0_Wv,N1_Wv,zeta_Wv,C_Wv)
end function diff_mu2Wv

function sigma_BO(R)
implicit none
real(prec) :: sigma_BO
real(prec),intent(in) :: R
sigma_BO = formSigma(R,nS_BO,aS_BO)
end function sigma_BO

function sigma_AD(R)
implicit none
real(prec) :: sigma_AD
real(prec),intent(in) :: R
sigma_AD = formSigma(R,nS_AD,aS_AD)
end function sigma_AD

function sigma_REL(R)
implicit none
real(prec) :: sigma_REL
real(prec),intent(in) :: R
sigma_REL = formSigma(R,nS_REL,aS_REL)
end function sigma_REL

function sigma_QED(R)
implicit none
real(prec) :: sigma_QED
real(prec),intent(in) :: R
sigma_QED = formSigma(R,nS_QED,aS_QED)
end function sigma_QED

function sigma(R)
implicit none
real(prec) :: sigma
real(prec),intent(in) :: R
real(prec) :: elms(4),maxelm
elms(1) = formSigma(R,nS_BO,aS_BO)
elms(2) = formSigma(R,nS_AD,aS_AD)
elms(3) = formSigma(R,nS_REL,aS_REL)
elms(4) = formSigma(R,nS_QED,aS_QED)
maxelm = maxval(elms)
sigma = sqrt(sum((elms*(1._prec/maxelm))**2))*maxelm
end function sigma

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

function formShort(R,M,I0,I1,a,P)
!short-range formula
!
! sum_{k=1}^{M} exp(-a(k) R) sum_{i=I0}^{I1} P(i,k) R^i
!
implicit none
real(prec) :: formShort
real(prec),intent(in) :: R
integer,intent(in) :: M,I0,I1
real(prec),intent(in) :: a(M)
real(prec),intent(in) :: P(I0:I1,M)
real(prec) :: term
integer :: i,k
formShort = 0._prec
do k=1,M
   term = P(I1,k)
   do i=I1-1,I0,-1
      term = term*R + P(i,k)
   enddo
   formShort = formShort + exp(-a(k)*R)*term
enddo
if(I0/=0) formShort = formShort*R**I0
end function formShort

function diffShort(R,M,I0,I1,a,P)
!derivative of the short-range formula
!
! sum_{k=1}^{M} exp(-a(k) R) sum_{i=I0}^{I1} P(i,k) (i/R - a(k)) R^i
!
implicit none
real(prec) :: diffShort
real(prec),intent(in) :: R
integer,intent(in) :: M,I0,I1
real(prec),intent(in) :: a(M)
real(prec),intent(in) :: P(I0:I1,M)
real(prec) :: iR,term
integer :: i,k
iR = 1._prec/R
diffShort = 0._prec
do k=1,M
   term = P(I1,k)*(I1*iR - a(k))
   do i=I1-1,I0,-1
      term = term*R + P(i,k)*(i*iR - a(k))
   enddo
   diffShort = diffShort + exp(-a(k)*R)*term
enddo
if(I0/=0) diffShort = diffShort*R**I0
end function diffShort

function formLong(R,ret_type,N0,N1,zeta,C)
!long-range formula
!
! - sum_{n=N0}^{N1} C(n) dampTT(n,zeta R)/R^n
!
! damping of the leading term
! ret_type = 0  ->   1     + shiftTT = dampTT (no retardation)
! ret_type = 1  ->   retar + shiftTT          (BO [C6] potential)
! ret_type = 2  ->   0     + shiftTT          (REL [C4] and QED [C3] potentials)
!
implicit none
real(prec) :: formLong
real(prec),intent(in) :: R
integer,intent(in) :: ret_type,N0,N1
real(prec),intent(in) :: zeta
real(prec),intent(in) :: C(N0:N1)
real(prec) :: iR,zetaR
integer :: n
iR = 1._prec/R
zetaR = zeta*R
if(N0<N1) then
   formLong = - C(N1)*dampTT(N1,zetaR)
   do n=N1-1,N0+1,-1
      if(C(n)/=ZERO) then
         formLong = formLong*iR - C(n)*dampTT(n,zetaR)
      else
         formLong = formLong*iR
      endif
   enddo
else
   formLong = 0._prec
endif
select case(ret_type)
case(0)
   formLong = formLong*iR - C(N0)*dampTT(N0,zetaR)
case(1)
   formLong = formLong*iR - C(N0)*(retar(R)+shftTT(N0,zetaR))
case(2)
   formLong = formLong*iR - C(N0)*shftTT(N0,zetaR)
end select
formLong = formLong*iR**N0
end function formLong

function diffLong(R,ret_type,N0,N1,zeta,C)
!derivative of the long-range formula
!
! + sum_{n=N0}^{N1} C(n) (n*dampTT(n,zeta R) - zeta R diffTT(n,zeta R))/R^(n+1)
!
implicit none
real(prec) :: diffLong
real(prec),intent(in) :: R
integer,intent(in) :: ret_type,N0,N1
real(prec),intent(in) :: zeta
real(prec),intent(in) :: C(N0:N1)
real(prec) :: iR,zetaR
integer :: n
if(ret_type/=0) then
   write(*,*) 'ERROR! diffLong called with non-zero ret_type'
   stop
endif
iR = 1._prec/R
zetaR = zeta*R
diffLong = C(N1)*(N1*dampTT(N1,zetaR) - zetaR*diffTT(N1,zetaR))
do n=N1-1,N0,-1
   if(C(n)/=ZERO) then
      diffLong = diffLong*iR + C(n)*(n*dampTT(n,zetaR) - zetaR*diffTT(n,zetaR))
   else
      diffLong = diffLong*iR
   endif
enddo
diffLong = diffLong*iR**(N0+1)
end function diffLong

function formSigma(R,n,aS)
!sigma formula
!
! exp(-a(0) R) s(0) + sum_{i=1}^{n} exp(-a(i) R^2) s(i)
!
implicit none
real(prec) :: formSigma
real(prec),intent(in) :: R
integer,intent(in) :: n
real(prec),intent(in) :: aS(2,0:n)
real(prec) :: R2
integer :: i
R2 = R*R
formSigma = exp(-aS(1,0)*R)*aS(2,0)
do i=1,n
   formSigma = formSigma + exp(-aS(1,i)*R2)*aS(2,i)
enddo
end function formSigma

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

function dampTT(n,x)
!Tang-Toennies damping function
!
! dampTT(n,x) = 1 - exp(-x) sum_{i=0}^{n} x^i/i!
!
!             = exp(-x) sum_{i=n+1}^{infty} x^i/i!
!
implicit none
real(prec) :: dampTT
integer,intent(in) :: n
real(prec),intent(in) :: x
real(prec),parameter :: dampTT_MINlimit = tiny(0._prec)
real(prec),parameter :: dampTT_MIDlimit(0:16) = (/&
     0.10536052_prec,&
     0.53181161_prec, 1.1020653_prec, 1.7447696_prec, 2.4325910_prec,&
     3.1518980_prec , 3.8947668_prec, 4.6561182_prec, 5.4324681_prec,&
     6.2213046_prec , 7.0207466_prec, 7.8293420_prec, 8.6459425_prec,&
     9.4696212_prec ,10.299617_prec ,11.135297_prec ,11.976127_prec  /)
real(prec),parameter :: dampTT_MAXlimit(0:16) = (/&
     35.350506_prec,&
     39.040395_prec, 42.189070_prec, 45.049320_prec, 47.719353_prec,&
     50.250606_prec, 52.674332_prec, 55.011322_prec, 57.276297_prec,&
     59.480153_prec, 61.631241_prec, 63.736139_prec, 65.800138_prec,&
     67.827579_prec, 69.822072_prec, 71.786664_prec, 73.723952_prec /)
real(prec) :: term,suma
integer :: i
if(x<dampTT_MINlimit) then
   dampTT = 0._prec
elseif(x<dampTT_MIDlimit(n)) then
   term = 1._prec
   do i=1,n
      term = term*x/i
   enddo
   suma = 0._prec
   do i=n+1,2*n+17
      term = term*x/i
      suma = suma + term
   enddo
   dampTT = exp(-x)*suma
elseif(x<dampTT_MAXlimit(n)) then
   term = 1._prec
   suma = 1._prec
   do i=1,n
      term = term*x/i
      suma = suma + term
   enddo
   dampTT = 1._prec - exp(-x)*suma
else
   dampTT = 1._prec
endif
end function dampTT

function shftTT(n,x)
!shifted Tang-Toennies damping function
!
! shftTT(n,x) = dampTT(n,x) - 1
!
!             = - exp(-x) sum_{i=0}^{n} x^i/i!
!
implicit none
real(prec) :: shftTT
integer,intent(in) :: n
real(prec),intent(in) :: x
real(prec),parameter :: shftTT_MAXlimit = 709.782712893384_prec
real(prec) :: term,suma
integer :: i
if(x<shftTT_MAXlimit) then
   term = 1._prec
   suma = 1._prec
   do i=1,n
      term = term*x/i
      suma = suma + term
   enddo
   shftTT = -exp(-x)*suma
else
   shftTT = 0._prec
endif
end function shftTT

function diffTT(n,x)
!derivative of the Tang-Toennies damping function
!
! diffTT(n,x) = exp(-x) x^n/n!
!
implicit none
real(prec) :: diffTT
integer,intent(in) :: n
real(prec),intent(in) :: x
real(prec),parameter :: diffTT_MAXlimit = 709.782712893384_prec
real(prec) :: term
integer :: i
if(x<diffTT_MAXlimit) then
   term = 1._prec
   do i=1,n
      term = term*x/i
   enddo
   diffTT = exp(-x)*term
else
   diffTT = 0._prec
endif
end function diffTT

function retar(R)
implicit none
real(prec) :: retar
real(prec),intent(in) :: R
real(prec),parameter :: A(5) = (/&
     6.169870091028090e-2_prec ,&
     8.281954197590046e-4_prec ,&
     2.936387935854446e-6_prec ,&
     4.020288241131003e-9_prec ,&
     2.948899823358023e-12_prec/)
real(prec),parameter :: B(6) = (/&
     6.169870091028090e-2_prec ,&
     8.523723778811090e-4_prec ,&
     4.032972820432622e-6_prec ,&
     9.969788378653687e-9_prec ,&
     1.224004819647673e-11_prec,&
     8.978131367598067e-15_prec/)
real(prec) :: num,den
num =         A(5)
num = num*R + A(4)
num = num*R + A(3)
num = num*R + A(2)
num = num*R + A(1)
num = num*R + 1._prec
den =         B(6)
den = den*R + B(5)
den = den*R + B(4)
den = den*R + B(3)
den = den*R + B(2)
den = den*R + B(1)
den = den*R + 1._prec
retar = num/den
end function retar

end module potentials
