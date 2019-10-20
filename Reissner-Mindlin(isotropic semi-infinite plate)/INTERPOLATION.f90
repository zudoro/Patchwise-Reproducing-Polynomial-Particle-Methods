MODULE INTERPOLATION

  USE NEWTYPE

  IMPLICIT NONE

    INTERFACE PHYSF
        MODULE PROCEDURE PHYSFVEC, PHYSF2D
    END INTERFACE PHYSF

CONTAINS

SUBROUTINE PHYSFVEC(VECVAL, VEC, NVEC, INDX, IORDER, INTSUPP)

    INTEGER, INTENT(IN) :: NVEC, IORDER
    TYPE(FVALUE), INTENT(OUT) :: VECVAL(NVEC)
    TYPE(POINT2D), INTENT(IN) :: VEC(NVEC)
    TYPE(RECPATCH), INTENT(IN) :: INTSUPP
    TYPE(INT2D), INTENT(IN) :: INDX
    TYPE(POINT2D) :: TMP_PT
    REAL*8 :: X_FACTOR, Y_FACTOR
    TYPE(FVALUE) :: TMP_SF
    INTEGER :: I

    X_FACTOR = IORDER/(INTSUPP%PT2%X-INTSUPP%PT1%X)
    Y_FACTOR = IORDER/(INTSUPP%PT3%Y-INTSUPP%PT1%Y)

    DO I=1, NVEC
        TMP_PT = POINT2D(X_FACTOR*(VEC(I)%X-INTSUPP%PT1%X), Y_FACTOR*(VEC(I)%Y-INTSUPP%PT1%Y))

        TMP_SF = REFSF(TMP_PT,INDX,IORDER)

        VECVAL(I)%D00 = TMP_SF%D00
        VECVAL(I)%D10 = X_FACTOR*TMP_SF%D10
        VECVAL(I)%D01 = Y_FACTOR*TMP_SF%D01
    ENDDO

END SUBROUTINE PHYSFVEC

TYPE(FVALUE) FUNCTION PHYSF2D(PT,INDX,IORDER,INTSUPP)

    TYPE(POINT2D), INTENT(IN) :: PT
    TYPE(INT2D), INTENT(IN) :: INDX
    TYPE(RECPATCH), INTENT(IN) :: INTSUPP
    INTEGER, INTENT(IN) :: IORDER
    TYPE(POINT2D) :: TMP_PT
    REAL(8) :: X_FACTOR, Y_FACTOR
    TYPE(FVALUE) :: TMP_SF

    X_FACTOR = IORDER/(INTSUPP%PT2%X-INTSUPP%PT1%X)
    Y_FACTOR = IORDER/(INTSUPP%PT3%Y-INTSUPP%PT1%Y)

    TMP_PT = POINT2D(X_FACTOR*(PT%X-INTSUPP%PT1%X), Y_FACTOR*(PT%Y-INTSUPP%PT1%Y))

    TMP_SF = REFSF(TMP_PT, INDX, IORDER)

    PHYSF2D%D00 = TMP_SF%D00
    PHYSF2D%D10 = X_FACTOR*TMP_SF%D10
    PHYSF2D%D01 = Y_FACTOR*TMP_SF%D01

END FUNCTION PHYSF2D


TYPE(FVALUE) FUNCTION REFSF(PT,INDX,IORDER)

      TYPE(POINT2D), INTENT(IN) :: PT
      TYPE(INT2D), INTENT(IN) :: INDX
      INTEGER, INTENT(IN) :: IORDER

      SELECT CASE (IORDER)
      CASE (2)
        REFSF%D00 = REF_1D_K2(PT%X,INDX%A,0)*REF_1D_K2(PT%Y,INDX%B,0)
        REFSF%D10 = REF_1D_K2(PT%X,INDX%A,1)*REF_1D_K2(PT%Y,INDX%B,0)
        REFSF%D01 = REF_1D_K2(PT%X,INDX%A,0)*REF_1D_K2(PT%Y,INDX%B,1)
      CASE(4)
        REFSF%D00 = REF_1D_K4(PT%X,INDX%A,0)*REF_1D_K4(PT%Y,INDX%B,0)
        REFSF%D10 = REF_1D_K4(PT%X,INDX%A,1)*REF_1D_K4(PT%Y,INDX%B,0)
        REFSF%D01 = REF_1D_K4(PT%X,INDX%A,0)*REF_1D_K4(PT%Y,INDX%B,1)
      CASE(6)
        REFSF%D00 = REF_1D_K6(PT%X,INDX%A,0)*REF_1D_K6(PT%Y,INDX%B,0)
        REFSF%D10 = REF_1D_K6(PT%X,INDX%A,1)*REF_1D_K6(PT%Y,INDX%B,0)
        REFSF%D01 = REF_1D_K6(PT%X,INDX%A,0)*REF_1D_K6(PT%Y,INDX%B,1)
      CASE(8)
        REFSF%D00 = REF_1D_K8(PT%X,INDX%A,0)*REF_1D_K8(PT%Y,INDX%B,0)
        REFSF%D10 = REF_1D_K8(PT%X,INDX%A,1)*REF_1D_K8(PT%Y,INDX%B,0)
        REFSF%D01 = REF_1D_K8(PT%X,INDX%A,0)*REF_1D_K8(PT%Y,INDX%B,1)
      END SELECT

END FUNCTION REFSF


  REAL(8) FUNCTION REF_1D_K2(X,INDX,DIFF)

    REAL(8), INTENT(IN) :: X
    INTEGER, INTENT(IN) :: INDX, DIFF

    REF_1D_K2 = 0.0D0

    IF (DIFF .EQ. 0) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K2 = ((-2+X)*(-1+X))*(1/2.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K2 = -((-2+X)*X)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K2 = ((-1+X)*X)*(1/2.D0)
      END IF
    END IF
    IF (DIFF .EQ. 1) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K2 = (-3+2*X)*(1/2.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K2 = -2.0D0*(-1+X)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K2 = (-1+2*X)*(1/2.D0)
      END IF
    END IF
    IF (DIFF .EQ. 2) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K2 = 1.0D0
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K2 = -1.0D0
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K2 = 1.0D0
      END IF
    END IF

  END FUNCTION REF_1D_K2


  REAL(8) FUNCTION REF_1D_K4(X,INDX,DIFF)

    REAL(8), INTENT(IN) :: X
    INTEGER, INTENT(IN) :: INDX, DIFF


    REF_1D_K4 = 0.0D0

    IF (DIFF .EQ. 0) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K4 = ((-4+X)*(-3+X)*(-2+X)*(-1+X))*(1/24.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K4 = -((-4+X)*(-3+X)*(-2+X)*X)*(1/6.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K4 = ((-4+X)*(-3+X)*(-1+X)*X)*(1/4.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K4 = -((-4+X)*(-2+X)*(-1+X)*X)*(1/6.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K4 = ((-3+X)*(-2+X)*(-1+X)*X)*(1/24.D0)
      END IF
    END IF
    IF (DIFF .EQ. 1) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K4 = ((-5+2*X)*(5-5*X+X**2))*(1/12.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K4 = (24-52*X+27*X**2-4*X**3)*(1/6.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K4 = ((-2+X)*(3-8*X+2*X**2))*(1/2.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K4 = (8-28*X+21*X**2-4*X**3)*(1/6.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K4 = ((-3+2*X)*(1-3*X+X**2))*(1/12.D0)
      END IF
    END IF
    IF (DIFF .EQ. 2) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K4 = (35-30*X+6*X**2)*(1.0D0/12.0D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K4 = (-26+27*X-6*X**2)*(1.0D0/3.0D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K4 = (19-24*X+6*X**2)*(1.0D0/2.0D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K4 = (-14+21*X-6*X**2)*(1.0D0/3.0D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K4 = (11-18*X+6*X**2)*(1.0D0/12.0D0)
      END IF
    END IF

  END FUNCTION REF_1D_K4


  REAL(8) FUNCTION REF_1D_K6(X,INDX,DIFF)

    REAL(8), INTENT(IN) :: X
    INTEGER, INTENT(IN) :: INDX, DIFF


    REF_1D_K6 = 0.0D0

    IF (DIFF .EQ. 0) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K6 = ((-6+X)*(-5+X)*(-4+X)*(-3+X)*(-2+X)*(-1+X))*(1/720.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K6 = -((-6+X)*(-5+X)*(-4+X)*(-3+X)*(-2+X)*X)*(1/120.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K6 = ((-6+X)*(-5+X)*(-4+X)*(-3+X)*(-1+X)*X)*(1/48.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K6 = -((-6+X)*(-5+X)*(-4+X)*(-2+X)*(-1+X)*X)*(1/36.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K6 = ((-6+X)*(-5+X)*(-3+X)*(-2+X)*(-1+X)*X)*(1/48.D0)
      ELSE IF (INDX .EQ. 5) THEN
        REF_1D_K6 = -((-6+X)*(-4+X)*(-3+X)*(-2+X)*(-1+X)*X)*(1/120.D0)
      ELSE IF (INDX .EQ. 6) THEN
        REF_1D_K6 = ((-5+X)*(-4+X)*(-3+X)*(-2+X)*(-1+X)*X)*(1/720.D0)
      END IF
    END IF
    IF (DIFF .EQ. 1) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K6 = ((-7+2*X)*(252-392*X+203*X**2-42*X**3+3*X**4))*(1/720.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K6 = (360-1044*X+870*X**2-310*X**3+50*X**4-3*X**5)*(1/60.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K6 = (-360+1404*X-1383*X**2+548*X**3-95*X**4+6*X**5)*(1/48.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K6 = -((-3+X)*(2-6*X+X**2)*(20-18*X+3*X**2))*(1/18.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K6 = (-180+792*X-921*X**2+428*X**3-85*X**4+6*X**5)*(1/48.D0)
      ELSE IF (INDX .EQ. 5) THEN
        REF_1D_K6 = (72-324*X+390*X**2-190*X**3+40*X**4-3*X**5)*(1/60.D0)
      ELSE IF (INDX .EQ. 6) THEN
        REF_1D_K6 = ((-5+2*X)*(24-100*X+95*X**2-30*X**3+3*X**4))*(1/720.D0)
      END IF
    END IF
    IF (DIFF .EQ. 2) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K6 = (1624-2205*X+1050*X**2-210*X**3+15*X**4)*(1.0D0/360.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K6 = (-1044+1740*X-930*X**2+200*X**3-15*X**4)*(1.0D0/60.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K6 = (702-1383*X+822*X**2-190*X**3+15*X**4)*(1.0D0/24.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K6 = (-508+1116*X-726*X**2+180*X**3-15*X**4)*(1.0D0/18.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K6 = (396-921*X+642*X**2-170*X**3+15*X**4)*(1.0D0/24.D0)
      ELSE IF (INDX .EQ. 5) THEN
        REF_1D_K6 = (-324+780*X-570*X**2+160*X**3-15*X**4)*(1.0D0/60.D0)
      ELSE IF (INDX .EQ. 6) THEN
        REF_1D_K6 = (274-675*X+510*X**2-150*X**3+15*X**4)*(1.0D0/360.D0)
      END IF
    END IF

  END FUNCTION REF_1D_K6


  REAL(8) FUNCTION REF_1D_K8(X,INDX,DIFF)

    REAL(8), INTENT(IN) :: X
    INTEGER, INTENT(IN) :: INDX, DIFF


    REF_1D_K8 = 0.0D0

    IF (DIFF .EQ. 0) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K8 = ((X-1)*(X-2)*(X-3)*(X-4)*(X-5)*(X-6)*(X-7)*(X-8))*(1/40320.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K8 = -((X)*(X-2)*(X-3)*(X-4)*(X-5)*(X-6)*(X-7)*(X-8))*(1/5040.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K8 = ((X)*(X-1)*(X-3)*(X-4)*(X-5)*(X-6)*(X-7)*(X-8))*(1/1440.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K8 = -((X)*(X-1)*(X-2)*(X-4)*(X-5)*(X-6)*(X-7)*(X-8))*(1/720.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K8 = ((X)*(X-1)*(X-2)*(X-3)*(X-5)*(X-6)*(X-7)*(X-8))*(1/576.D0)
      ELSE IF (INDX .EQ. 5) THEN
        REF_1D_K8 = -((X)*(X-1)*(X-2)*(X-3)*(X-4)*(X-6)*(X-7)*(X-8))*(1/720.D0)
      ELSE IF (INDX .EQ. 6) THEN
        REF_1D_K8 = ((X)*(X-1)*(X-2)*(X-3)*(X-4)*(X-5)*(X-7)*(X-8))*(1/1440.D0)
      ELSE IF (INDX .EQ. 7) THEN
        REF_1D_K8 = -((X)*(X-1)*(X-2)*(X-3)*(X-4)*(X-5)*(X-6)*(X-8))*(1/5040.D0)
      ELSE IF (INDX .EQ. 8) THEN
        REF_1D_K8 = ((X)*(X-1)*(X-2)*(X-3)*(X-4)*(X-5)*(X-6)*(X-7))*(1/40320.D0)
      END IF
    END IF
    IF (DIFF .EQ. 1) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K8 = ((-9+2*X)*(3044-5886*X+4299*X**2-1539*X**3+288*X**4-27*X**5+X**6))*(1/10080.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K8 = (40320-138528*X+146580*X**2-73696*X**3+20125*X**4-3066*X**5+245*X**6-8*X**7)*(1/5040.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K8 = (-10080+44712*X-55059*X**2+30578*X**3-8950*X**4+1434*X**5-119*X**6+4*X**7)*(1/720.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K8 = (13440-64096*X+86076*X**2-51456*X**3+15975*X**4-2682*X**5+231*X**6-8*X**7)*(1/720.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K8 = ((-4+X)*(630-2952*X+3633*X**2-1840*X**3+435*X**4-48*X**5+2*X**6))*(1/144.D0)
      ELSE IF (INDX .EQ. 5) THEN
        REF_1D_K8 = (8064-40608*X+58692*X**2-38176*X**3+12905*X**4-2346*X**5+217*X**6-8*X**7)*(1/720.D0)
      ELSE IF (INDX .EQ. 6) THEN
        REF_1D_K8 = (-3360+17144*X-25245*X**2+16818*X**3-5850*X**4+1098*X**5-105*X**6+4*X**7)*(1/720.D0)
      ELSE IF (INDX .EQ. 7) THEN
        REF_1D_K8 = (5760-29664*X+44268*X**2-30016*X**3+10675*X**4-2058*X**5+203*X**6-8*X**7)*(1/5040.D0)
      ELSE IF (INDX .EQ. 8) THEN
        REF_1D_K8 = ((-7+2*X)*(180-882*X+1155*X**2-637*X**3+168*X**4-21*X**5+X**6))*(1/10080.D0)
      END IF
    END IF
    IF (DIFF .EQ. 2) THEN
      IF (INDX .EQ. 0) THEN
        REF_1D_K8 = (59062-100926*X+67347*X**2-22680*X**3+4095*X**4-378*X**5+14*X**6)*(1.0D0/10080.D0)
      ELSE IF (INDX .EQ. 1) THEN
        REF_1D_K8 = (-69264+146580*X-110544*X**2+40250*X**3-7665*X**4+735*X**5-28*X**6)*(1.0D0/2520.D0)
      ELSE IF (INDX .EQ. 2) THEN
        REF_1D_K8 = (22356-55059*X+45867*X**2-17900*X**3+3585*X**4-357*X**5+14*X**6)*(1.0D0/360.D0)
      ELSE IF (INDX .EQ. 3) THEN
        REF_1D_K8 = (-32048+86076*X-77184*X**2+31950*X**3-6705*X**4+693*X**5-28*X**6)*(1.0D0/360.D0)
      ELSE IF (INDX .EQ. 4) THEN
        REF_1D_K8 = (12438-34968*X+32979*X**2-14320*X**3+3135*X**4-336*X**5+14*X**6)*(1.0D0/144.D0)
      ELSE IF (INDX .EQ. 5) THEN
        REF_1D_K8 = (-20304+58692*X-57264*X**2+25810*X**3-5865*X**4+651*X**5-28*X**6)*(1.0D0/360.D0)
      ELSE IF (INDX .EQ. 6) THEN
        REF_1D_K8 = (8572-25245*X+25227*X**2-11700*X**3+2745*X**4-315*X**5+14*X**6)*(1.0D0/360.D0)
      ELSE IF (INDX .EQ. 7) THEN
        REF_1D_K8 = (-14832+44268*X-45024*X**2+21350*X**3-5145*X**4+609*X**5-28*X**6)*(1.0D0/2520.D0)
      ELSE IF (INDX .EQ. 8) THEN
        REF_1D_K8 = (6534-19698*X+20307*X**2-9800*X**3+2415*X**4-294*X**5+14*X**6)*(1.0D0/10080.D0)
      END IF
    END IF

  END FUNCTION REF_1D_K8

END MODULE INTERPOLATION
