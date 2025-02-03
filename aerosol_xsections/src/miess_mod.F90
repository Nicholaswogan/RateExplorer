module miess_mod

  use system_mod
  use adgaquad_types_mod,  only: nf

  implicit none

  private
  
  public :: miess

  contains

  !! This subroutine computes mie scattering by a stratified sphere,
  !! i.e. a particle consisting of a spherical core surrounded by a
  !! Spherical shell.  The basic code used was that described in the
  !! report: " Subroutines for computing the parameters of the
  !! electromagnetic radiation scattered by a sphere " J.V. Dave,
  !! IBM Scientific Center, Palo Alto , California.
  !! Report No. 320 - 3236 .. May 1968 .
  !!
  !! The modifications for stratified spheres are described in
  !!   Toon and Ackerman, Appl. Optics, in press, 1981
  !!
  !! The definitions for the output parameters can be found in "Light
  !! scattering by small particles, H.C.Van de Hulst, John Wiley '
  !! Sons, Inc., New York, 1957".
  !!
  !! The latest version of this program is available by anonymous ftp 
  !! from climate.gsfc.nasa.gov in directory pub/wiscombe.
  !!  
  !! This code has been modified by E.T. Wolf (2013) to be joined 
  !! with the fractal meanfield scattering code of Botet (1997).
  !! Here primary code outputs are the Mie wave coefficients an,bn
  subroutine miess( RCORE, RSHELL, WVNO, RINDSH, RINDCO, MU, &
     &                   NUMANG, MAXANG, ACOE_OUT, BCOE_OUT, NSTOP) 

    ! Rework arguments, an, bn arrays output, nstop out 

    ! Input Arguments
    real(kind=f), intent(in)     ::  RCORE
    real(kind=f), intent(in)     ::  RSHELL
    real(kind=f), intent(in)     ::  WVNO
    complex(kind=f), intent(in)  ::  RINDSH
    complex(kind=f), intent(in)  ::  RINDCO
    integer, intent(in)          ::  NUMANG
    integer, intent(in)          ::  MAXANG
    real(kind=f), intent(in)     ::  MU(NUMANG)
    complex(kind=f),intent(out)   :: ACOE_OUT(nf)
    complex(kind=f),intent(out)   :: BCOE_OUT(nf)
    integer, intent(out)          :: NSTOP

    ! Local declarations
    integer, parameter      :: MXANG = 100
    integer, parameter      :: LL = 1000
    real(kind=f), parameter :: ZERO = 0.0_f
    real(kind=f), parameter :: ONE = 1.0_f
    real(kind=f), parameter :: TWO = 2.0_f
!    logical                :: INPERR, PASS1
    integer                 :: J, K, M, N, NMX1, NMX2, NN
    real(kind=f)            :: AA, AIM, AM1IM, AM1RE, ARE, BB, BIM, BM1IM
    real(kind=f)            :: BM1RE, BRE, CC, COSX1, COSX4, DD, DENOM
    real(kind=f)            :: GQSC, QBS, QEXT, QSCA
    real(kind=f)            :: D21(MAXANG,2)
    real(kind=f)            :: S21(MAXANG,2)
    real(kind=f)            :: M1(MAXANG,2)
    real(kind=f)            :: M2(MAXANG,2)
    real(kind=f)            :: DGQSC, DQEXT, DQSCA, E2Y1
    real(kind=f)            :: EY1, EY1MY4, EY1PY4, EY4, FOURPI, PINUM
    real(kind=f)            :: RMM, RX, SINX1, SINX4, TOLER, X1, X4
    real(kind=f)            :: XCORE, XSHELL, Y1, Y4
    complex(kind=f)         :: AC, ACOE, ACOEM1, BC, BCOE, BCOEM1, CI, CZERO
    complex(kind=f)         :: DH1, DH2, DH4, DUMMY, DUMSQ, K1, K2, K3
    complex(kind=f)         :: P24H21, P24H24, RRFX, SBACK, WM1
    real(kind=f)            :: PI(MXANG,3)
    real(kind=f)            :: SI2THT(MXANG )
    real(kind=f)            :: T(5)
    real(kind=f)            :: TA(4)
    real(kind=f)            :: TAU(MXANG,3)
    complex(kind=f)         :: ACAP(LL)
    complex(kind=f)         :: S1(MXANG,2)
    complex(kind=f)         :: S2(MXANG,2)      
    complex(kind=f)         :: U(8) 
    complex(kind=f)         :: W(3,LL)
    complex(kind=f)         :: WFN(2)
    complex(kind=f)         :: Z(4)

!    SAVE  PINUM, PASS1

 !   data PASS1 / .True. /
    data TOLER / 1.e-6_f /
    data CZERO / (0.0_f, 0.0_f) /
    data CI    / (0.0_f, 1.0_f) /

    !
    ! Start Code
    !
    !IF (PASS1) THEN
      PINUM  = TWO*ASIN(ONE)
   !   PASS1  = .False.
   ! ENDIF

    XSHELL = RSHELL*WVNO
    XCORE  = RCORE*WVNO
    T(1) = XSHELL*ABS(RINDSH)
    NMX1   = 1.1_f*T(1)
    NMX2   = T(1)

    IF (NMX1 .LE. 150) THEN
      NMX1   = 150
      NMX2   = 135
    END IF

    !
    ! Check input arguments for gross errors
    !
    IF (WVNO .LE. 0.0) write(*,*) "ERROR miess_mod: wavenumber less than zero"
    IF (RSHELL .LE. 0.0) write(*,*) "ERROR miess_mod: RSHELL less than zero"
    IF (RCORE .LE. 0.0)  write(*,*) "ERROR miess_mod: RCORE less than zero"
    IF (RCORE .GT. RSHELL) write(*,*) "ERROR miess_mod: RCORE greater than RSHELL"
    IF (REAL(RINDSH) .LE. 0.0) write(*,*) "ERROR miess_mod: RINDSH REAL less than zero"
    IF (AIMAG(RINDSH) .GT. 0.0) write(*,*) "ERROR miess_mod: RINDSH IMAGINARY greater than zero"
    IF (REAL(RINDCO) .LE. 0.0) write(*,*) "ERROR miess_mod: RINDCO REAL less than zero"
    IF (AIMAG(RINDCO) .GT. 0.0) write(*,*) "ERROR miess_mod: RINDCO IMAGINARY greater than zero"
    IF (NUMANG .LT. 0 ) write(*,*) "ERROR miess_mod: NUMANG less than zero"
    IF (NUMANG .GT. MAXANG ) write(*,*) "ERROR miess_mod: NUMANG greater than MAXANG"
    IF (NMX1+1 .GT. LL ) write(*,*) "ERROR miess_mod: NMX1+1 greater than LL"
    
    do J=1, NUMANG
      IF (MU(J) .LT. -TOLER .OR. MU(J) .GT. 1.0_f+TOLER) write(*,*) "ERROR miess_mod: MU invalid"
    enddo

    K1     = RINDCO*WVNO
    K2     = RINDSH*WVNO
    K3     = DCMPLX(WVNO)
    Z(1)   = RINDSH*XSHELL
    Z(2)   = XSHELL
    Z(3)   = RINDCO*XCORE
    Z(4)   = RINDSH*XCORE
    X1     = DBLE(Z(1))
    Y1     = DIMAG(Z(1))
    X4     = DBLE(Z(4))
    Y4     = DIMAG(Z(4))
    RX     = ONE/XSHELL

    ! Down-recurrence for A function
    ACAP(NMX1+1) = CZERO
    do M=1,3
      W(M,NMX1+ 1) = CZERO
    enddo

    RRFX = ONE/(RINDSH*XSHELL)
    do NN=NMX1,1,- 1
      ACAP(NN) = ((NN+1)*RRFX)-ONE/(((NN+1)*RRFX)+ACAP(NN+1))
      do M = 1, 3
        W(M,NN) = ((NN+1)/Z(M+1))-ONE/(((NN+1)/Z(M+1))+W(M,NN+1))
      enddo
    enddo


    do J=1,NUMANG
      SI2THT(J) = ONE-MU(J)**2
      PI(J,1)   = ZERO
      PI(J,2)   = ONE
      TAU(J,1)  = ZERO
      TAU(J,2)  = MU(J)
    enddo

    ! Initialization of homogeneous sphere
    T(1)   = COS(XSHELL)
    T(2)   = SIN(XSHELL)
    WM1    = DCMPLX(T(1),-T(2))
    WFN(1) = DCMPLX(T(2),T(1))
    TA(1)  = T(2)
    TA(2)  = T(1)
    WFN(2) = RX*WFN(1)-WM1
    TA(3)  = DBLE(WFN(2))
    TA(4)  = DIMAG(WFN(2))

    ! Initialization procedure for stratified sphere
    N      = 1
    SINX1  = SIN(X1)
    SINX4  = SIN(X4)
    COSX1  = COS(X1)
    COSX4  = COS(X4)
    EY1    = EXP(Y1)
    E2Y1   = EY1**2
    EY4    = EXP(Y4)
    EY1MY4 = EXP(Y1-Y4)
    EY1PY4 = EY1*EY4
    AA     = SINX4*(EY1PY4+EY1MY4)
    BB     = COSX4*(EY1PY4-EY1MY4)
    CC     = SINX1*(E2Y1+ONE)
    DD     = COSX1*(E2Y1-ONE)
    DENOM  = ONE+E2Y1*(4.0_f*SINX1**2-TWO+E2Y1)
    DUMMY  = DCMPLX((AA*CC+BB*DD)/DENOM, (BB*CC-AA*DD)/DENOM)
    DUMMY  = DUMMY*(ACAP(N)+N/Z(1))/(W(3,N)+N/Z(4))
    DUMSQ  = DUMMY**2

    P24H24 = 0.5_f+DCMPLX(SINX4**2-0.5_f, COSX4*SINX4)*EY4**2
    P24H21 = 0.5_f*DCMPLX(SINX1*SINX4-COSX1*COSX4, SINX1*COSX4+COSX1*SINX4)*EY1PY4   & 
           + 0.5_f*DCMPLX(SINX1*SINX4+COSX1*COSX4, -SINX1*COSX4+COSX1*SINX4)*EY1MY4
    DH1    = Z(1)/(ONE+CI*Z(1))-ONE/Z(1)
    DH2    = Z(2)/(ONE+CI*Z(2))-ONE/Z(2)
    DH4    = Z(4)/(ONE+CI*Z(4))-ONE/Z(4)
    P24H24 = P24H24/((DH4+N/Z(4))*(W(3,N)+N/Z(4)))
    P24H21 = P24H21/((DH1+N/Z(1))*(W(3,N)+N/Z(4)))

    U(1) = K3*ACAP(N)-K2*W(1,N)
    U(2) = K3*ACAP(N)-K2*DH2
    U(3) = K2*ACAP(N)-K3*W(1,N)
    U(4) = K2*ACAP(N)-K3*DH2
    U(5) = K1*W(3,N)-K2*W(2,N)
    U(6) = K2*W(3,N)-K1*W(2,N)
    U(7) = -CI*(DUMMY*P24H21-P24H24)
    U(8) = TA(3)/WFN(2)

    ACOE  = U(8)*(U(1)*U(5)*U(7)+K1*U(1)-DUMSQ*K3*U(5))  &
                /(U(2)*U(5)*U(7)+K1*U(2)-DUMSQ*K3*U(5))

    BCOE  = U(8)*(U(3)*U(6)*U(7)+K2*U(3)-DUMSQ*K2*U(6))  &
                /(U(4)*U(6)*U(7)+K2*U(4)-DUMSQ*K2*U(6))

    ACOE_OUT(N) = ACOE
    BCOE_OUT(N) = BCOE

    ACOEM1 = ACOE
    BCOEM1 = BCOE
    ARE    =  DBLE(ACOE)
    AIM    = DIMAG(ACOE)
    BRE    =  DBLE(BCOE)
    BIM    = DIMAG(BCOE)

    DQEXT  = 3.0_f*(ARE+BRE)
    DQSCA  = 3.0_f*(ARE**2+AIM**2+BRE**2+BIM**2)
    DGQSC  = ZERO
    SBACK  = 3.0_f*(ACOE-BCOE)
    RMM    = ONE

    AC  = 1.5_f*ACOE
    BC  = 1.5_f*BCOE
   
    do J=1,NUMANG
      S1(J,1) = AC*PI(J,2) + BC*TAU(J,2)
      S1(J,2) = AC*PI(J,2) - BC*TAU(J,2)
      S2(J,1) = BC*PI(J,2) + AC*TAU(J,2)
      S2(J,2) = BC*PI(J,2) - AC*TAU(J,2)
    enddo

    ! ***************** Start of Mie summing loop ******************
    N  = 2
 70 CONTINUE
 
    ! Recurrences for functions little-pi,little-tau of Mie theory
    T(1) = 2*N-1
    T(2) = N-1
    do J=1,NUMANG
      PI(J,3) = (T(1)*PI(J,2)*MU(J)-N*PI(J,1))/T(2)
      TAU(J,3) = MU(J)*(PI(J,3)-PI(J,1))-T(1)*SI2THT(J)*PI(J,2)+TAU(J,1)
    enddo

    ! Here set up homogeneous sphere
    WM1    = WFN(1)
    WFN(1) = WFN(2)
    WFN(2) = T(1)*RX*WFN(1)-WM1
    TA(1) =  DBLE(WFN(1))
    TA(2) = DIMAG(WFN(1))
    TA(3) =  DBLE(WFN(2))
    TA(4) = DIMAG(WFN(2))

    ! Here set up stratified sphere
    DH1    = -N/Z(1)+ONE/(N/Z(1)-DH1)
    DH2    = -N/Z(2)+ONE/(N/Z(2)-DH2)
    DH4    = -N/Z(4)+ONE/(N/Z(4)-DH4)
    P24H24 = P24H24/((DH4+N/Z(4))*(W(3,N)+N/Z(4)))
    P24H21 = P24H21/((DH1+N/Z(1))*(W(3,N)+N/Z(4)))
    DUMMY  = DUMMY*(ACAP(N)+N/Z(1))/(W(3,N)+N/Z(4))
    DUMSQ  = DUMMY**2

    U(1) = K3*ACAP(N)-K2*W(1,N)
    U(2) = K3*ACAP(N)-K2*DH2
    U(3) = K2*ACAP(N)-K3*W(1,N)
    U(4) = K2*ACAP(N)-K3*DH2
    U(5) = K1*W(3,N)-K2*W(2,N)
    U(6) = K2*W(3,N)-K1*W(2,N)
    U(7) = -CI*(DUMMY*P24H21-P24H24)
    U(8) = TA(3)/WFN(2)

    ACOE  = U(8)*(U(1)*U(5)*U(7)+K1*U(1)-DUMSQ*K3*U(5))   & 
                /(U(2)*U(5)*U(7)+K1*U(2)-DUMSQ*K3*U(5))

    BCOE  = U(8)*(U(3)*U(6)*U(7)+K2*U(3)-DUMSQ*K2*U(6))    &
                /(U(4)*U(6)*U(7)+K2*U(4)-DUMSQ*K2*U(6))

    ACOE_OUT(N) = ACOE
    BCOE_OUT(N) = BCOE

    ARE  = DBLE(ACOE)
    AIM  = DIMAG(ACOE)
    BRE  = DBLE(BCOE)
    BIM  = DIMAG(BCOE)

    ! Increment sums for efficiency factors
    AM1RE  = DBLE(ACOEM1)
    AM1IM  = DIMAG(ACOEM1)
    BM1RE  = DBLE(BCOEM1)
    BM1IM  = DIMAG(BCOEM1)
    T(4)   = (2*N-ONE)/(N*(N-ONE))
    T(2)   = (N-ONE)*(N+ONE)/N
    DGQSC  = DGQSC+T(2)*(AM1RE*ARE+AM1IM*AIM   &
           + BM1RE*BRE+BM1IM*BIM)   &
           + T(4)*(AM1RE*BM1RE+AM1IM*BM1IM)

    T(3)    = 2*N+1
    DQEXT   = DQEXT+T(3)*(ARE+BRE)
    T(4)  = ARE**2+AIM**2+BRE**2+BIM**2
    DQSCA   = DQSCA+T(3)*T(4)
    RMM     = -RMM
    SBACK   = SBACK+T(3)*RMM*(ACOE-BCOE)

    T(2) = N*(N+1)
    T(1) = T(3)/T(2)

    AC  = T(1)*ACOE
    BC  = T(1)*BCOE
    do J=1,NUMANG
      S1(J,1) = S1(J,1)+AC*PI(J,3)+BC*TAU(J,3)
      S2(J,1) = S2(J,1)+BC*PI(J,3)+AC*TAU(J,3)
    enddo

    ! Scattering matrix elements for supplements of 0-90 degree scattering angles submitted by user
    IF (MOD(N,2) .EQ. 0) THEN
      do J=1,NUMANG
        S1(J,2) = S1(J,2)-AC*PI(J,3)+BC*TAU(J,3)
        S2(J,2) = S2(J,2)-BC*PI(J,3)+AC*TAU(J,3)
      enddo
    ELSE
      do J=1,NUMANG
        S1(J,2) = S1(J,2)+AC*PI(J,3)-BC*TAU(J,3)
        S2(J,2) = S2(J,2)+BC*PI(J,3)-AC*TAU(J,3)
      enddo
    END IF

    ! Test for convergence of sums
    IF (T(4) .GE. 1.0e-14_f) THEN
      N=N+1
      !IF (N .GT. NMX2) CALL ERRMSG(
     !&       'MIELAY--Dimensions for W,ACAP not enough. Suggest'//
     !&       ' get detailed output, modify routine', .True. )
       IF (N .GT. NMX2) write(*,*) "N>NMX2"
       do J=1,NUMANG
         PI(J,1) = PI(J,2)
         PI(J,2) = PI(J,3)
         TAU(J,1) = TAU(J,2)
         TAU(J,2) = TAU(J,3)
       enddo

       ACOEM1 = ACOE
       BCOEM1 = BCOE

       GO TO 70

    END IF
    NSTOP = N
    ! ***************** End of summing loop ******************

    !Transform complex scattering amplitudes into elements of real scattering matrix

    do J=1,NUMANG
      do K=1,2
        M1(J,K)  = DBLE(S1(J,K))**2+DIMAG(S1(J,K))**2
        M2(J,K)  = DBLE(S2(J,K))**2+DIMAG(S2(J,K))**2
        S21(J,K) = DBLE(S1(J,K))*DBLE(S2(J,K))+DIMAG(S1(J,K))*DIMAG(S2(J,K))
        D21(J,K) = DIMAG(S1(J,K))*DBLE(S2(J,K))-DIMAG(S2(J,K))*DBLE(S1(J,K))
      enddo
    enddo

    T(1)   = TWO*RX**2
    QEXT   = T(1)*DQEXT
    QSCA   = T(1)*DQSCA
    GQSC   = TWO*T(1)*DGQSC
    SBACK  = 0.5_f*SBACK
    QBS    = (DBLE(SBACK)**2+DIMAG(SBACK)**2)/(PINUM*XSHELL**2)
  end subroutine miess

end module