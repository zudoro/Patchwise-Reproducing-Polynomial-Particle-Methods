MODULE INTERPOLATION

	USE GLBVAR

  IMPLICIT NONE

    INTERFACE PHYSF
        MODULE PROCEDURE PHYSFVEC, PHYSF2D
    END INTERFACE PHYSF

CONTAINS

SUBROUTINE PHYSFVEC(FT, PT, INDX, INTSUPP)

	TYPE(POINT2D), INTENT(IN) :: PT(:)
  TYPE(FVALUE), INTENT(OUT) :: FT(UBOUND(PT,1))
  TYPE(RECPATCH), INTENT(IN) :: INTSUPP
  TYPE(INT2D), INTENT(IN) :: INDX
  TYPE(TRANSFORM2D) :: TMP_PT
  TYPE(FVALUE) :: TMP_SF
  INTEGER :: I, N

	N = UBOUND(PT,1)

	DO I=1, N
		IF ((PT(I).IN.INTSUPP).EQV..TRUE.) THEN
			TMP_PT = PHYREC_TO_MSTREC(PT(I), INTSUPP, MSTREC)
			TMP_SF = REFSF(TMP_PT%PT, INDX)
			FT(I)%D00 = TMP_SF%D00
			FT(I)%D10 = TMP_SF%D10*TMP_PT%DXDX + TMP_SF%D01*TMP_PT%DYDX
			FT(I)%D01 = TMP_SF%D10*TMP_PT%DXDY + TMP_SF%D01*TMP_PT%DYDY
		ELSE
			FT(I) = FVALUE(0.D0,0.D0,0.D0)
		ENDIF
	ENDDO

END SUBROUTINE PHYSFVEC

TYPE(FVALUE) FUNCTION PHYSF2D(PT, INDX, INTSUPP)

	TYPE(POINT2D), INTENT(IN) :: PT
	TYPE(INT2D), INTENT(IN) :: INDX
	TYPE(RECPATCH), INTENT(IN) :: INTSUPP
	TYPE(TRANSFORM2D) :: TMP_PT
	TYPE(FVALUE) :: TMP_SF

	IF ((PT.IN.INTSUPP).EQV..TRUE.) THEN
		TMP_PT = PHYREC_TO_MSTREC(PT, INTSUPP, MSTREC)
		TMP_SF = REFSF(TMP_PT%PT, INDX)

    PHYSF2D%D00 = TMP_SF%D00
    PHYSF2D%D10 = TMP_SF%D10*TMP_PT%DXDX + TMP_SF%D01*TMP_PT%DYDX
    PHYSF2D%D01 = TMP_SF%D10*TMP_PT%DXDY + TMP_SF%D01*TMP_PT%DYDY
	ELSE
		PHYSF2D = FVALUE(0.D0,0.D0,0.D0)
	ENDIF

END FUNCTION PHYSF2D


TYPE(FVALUE) FUNCTION REFSF(PT,INDX)
  TYPE(POINT2D), INTENT(IN) :: PT
	TYPE(INT2D), INTENT(IN) :: INDX

  REFSF%D00 = REF1D(PT%X,INDX%A,0)*REF1D(PT%Y,INDX%B,0)
	REFSF%D10 = REF1D(PT%X,INDX%A,1)*REF1D(PT%Y,INDX%B,0)
	REFSF%D01 = REF1D(PT%X,INDX%A,0)*REF1D(PT%Y,INDX%B,1)

END FUNCTION REFSF


REAL*8 FUNCTION REF1D(X,INDX,DIFF)

	REAL*8, INTENT(IN)::X
	INTEGER, INTENT(IN)::INDX,DIFF

	IF(IORDER.EQ.2)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=((-1+X)*X)/2.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=1.0D0-X**2
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(X*(1+X))/2.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=-0.50D0+X
			ELSEIF(INDX.EQ.1)THEN
				REF1D=-2.0D0*X
			ELSEIF(INDX.EQ.2)THEN
				REF1D=0.50D0+X
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.3)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=-((-1+X)*(-1+3*X)*(1+3*X))/16.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(9*(-1+X)*(1+X)*(-1+3*X))/16.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(-9*(-1+X)*(1+X)*(1+3*X))/16.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=((1+X)*(-1+3*X)*(1+3*X))/16.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(1.0D0+18*X-27*X**2)/16.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(9*(-3.0D0-2*X+9*X**2))/16.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(-9*(-3.0D0+2*X+9*X**2))/16.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-1.0D0+18*X+27*X**2)/16.0D0
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.4)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(2.0D0*(-1.0D0+X)*(-0.50D0+X)*X*(0.50D0+X))/3.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-8.0D0*(-1.0D0+X)*(-0.50D0+X)*X*(1.0D0+X))/3.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=4*(-1.0D0+X)*(-0.50D0+X)*(0.50D0+X)*(1.0D0+X)
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-8.0D0*(-1.0D0+X)*X*(0.50D0+X)*(1.0D0+X))/3.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(2.0D0*(-0.50D0+X)*X*(0.50D0+X)*(1.0D0+X))/3.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(1/6.0D0)-(X/3)-(2*X**2)+(8*X**3/3)
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-4*(1-4*X-3*X**2+8*X**3))/3.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=2*X*(-5+8*X**2)
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-4*(-1-4*X+3*X**2+8*X**3))/3.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=-(1/6.0D0)-(X/3)+(2*X**2)+(8*X**3/3)
			ENDIF
		ELSEIF (DIFF.EQ.2) THEN
			IF(INDX.EQ.0)THEN
				REF1D = -1.D0/3.D0 - 4.D0*X + 8.D0*X**2
			ELSEIF(INDX.EQ.1)THEN
				REF1D = -(4.D0/3.D0)*(-4.D0 - 6.D0*X + 24.D0*X**2)
			ELSEIF(INDX.EQ.2)THEN
				REF1D = 32.D0*X**2 + 2.D0*(-5.D0 + 8.D0*X**2)
			ELSEIF(INDX.EQ.3)THEN
				REF1D = -(4.D0/3.D0)*(-4.D0 + 6.D0*X + 24.D0*X**2)
			ELSEIF(INDX.EQ.4)THEN
				REF1D = -1.D0/3.D0 + 4.D0*X + 8.D0*X**2
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.5)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=-((-1+X)*(-3+5*X)*(-1+5*X)*(1+5*X)*(3+5*X))/768.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(25*(-1+X)*(1+X)*(-3+5*X)*(-1+5*X)*(1+5*X))/768.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(-25*(-1+X)*(1+X)*(-3+5*X)*(-1+5*X)*(3+5*X))/384.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(25*(-1+X)*(1+X)*(-3+5*X)*(1+5*X)*(3+5*X))/384.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(-25*(-1+X)*(1+X)*(-1+5*X)*(1+5*X)*(3+5*X))/768.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=((1+X)*(-3+5*X)*(-1+5*X)*(1+5*X)*(3+5*X))/768.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(-9-500*X+750*X**2+2500*X**3-3125*X**4)/768.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(25*(5+156*X-390*X**2-300*X**3+625*X**4))/768.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(-25*(45+68*X-510*X**2-100*X**3+625*X**4))/384.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(25*(45-68*X-510*X**2+100*X**3+625*X**4))/384.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(-25*(5-156*X-390*X**2+300*X**3+625*X**4))/768.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(9-500*X-750*X**2+2500*X**3+3125*X**4)/768.0D0
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.6)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=((-1+X)*X*(-2+3*X)*(-1+3*X)*(1+3*X)*(2+3*X))/80.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-9*(-1+X)*X*(1+X)*(-2+3*X)*(-1+3*X)*(1+3*X))/40.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(9*(-1+X)*X*(1+X)*(-2+3*X)*(-1+3*X)*(2+3*X))/16.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(4-49*X**2+126*X**4-81*X**6)/4.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(9*(-1+X)*X*(1+X)*(-2+3*X)*(1+3*X)*(2+3*X))/16.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(-9*(-1+X)*X*(1+X)*(-1+3*X)*(1+3*X)*(2+3*X))/40.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(X*(1+X)*(-2+3*X)*(-1+3*X)*(1+3*X)*(2+3*X))/80.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(-4+8*X+135*X**2-180*X**3-405*X**4+486*X**5)/80.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-9*(-1+3*X+30*X**2-60*X**3-45*X**4+81*X**5))/20.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(9*(-4+24*X+39*X**2-156*X**3-45*X**4+162*X**5))/16.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-49*X+252*X**3-243*X**5)/2.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(9*(4+24*X-39*X**2-156*X**3+45*X**4+162*X**5))/16.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(-9*(1+3*X-30*X**2-60*X**3+45*X**4+81*X**5))/20.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(4+8*X-135*X**2-180*X**3+405*X**4+486*X**5)/80.0D0
			ENDIF
		ELSEIF(DIFF.EQ.2)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(8 + 270*X - 540*X**2 - 1620*X**3 + 2430*X**4)/80.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-9*(3 + 60*X - 180*X**2 - 180*X**3 + 405*X**4))/20.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(9*(24 + 78*X - 468*X**2 - 180*X**3 + 810*X**4))/16.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-49 + 756*X**2 - 1215*X**4)/2.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(9*(24 - 78*X - 468*X**2 + 180*X**3 + 810*X**4))/16.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(-9*(3 - 60*X - 180*X**2 + 180*X**3 + 405*X**4))/20.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(8 - 270*X - 540*X**2 + 1620*X**3 + 2430*X**4)/80.0D0
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.7)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=-((-1+X)*(-5+7*X)*(-3+7*X)*(-1+7*X)*(1+7*X)*(3+7*X)*(5+7*X))/92160.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(49*(-1+X)*(1+X)*(-5+7*X)*(-3+7*X)*(-1+7*X)*(1+7*X)*(3+7*X))/92160.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(-49*(-1+X)*(1+X)*(-5+7*X)*(-3+7*X)*(-1+7*X)*(1+7*X)*(5+7*X))/30720.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(49*(-1+X)*(1+X)*(-5+7*X)*(-3+7*X)*(-1+7*X)*(3+7*X)*(5+7*X))/18432.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(-49*(-1+X)*(1+X)*(-5+7*X)*(-3+7*X)*(1+7*X)*(3+7*X)*(5+7*X))/18432.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(49*(-1+X)*(1+X)*(-5+7*X)*(-1+7*X)*(1+7*X)*(3+7*X)*(5+7*X))/30720.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(-49*(-1+X)*(1+X)*(-3+7*X)*(-1+7*X)*(1+7*X)*(3+7*X)*(5+7*X))/92160.0D0
			ELSEIF(INDX.EQ.7)THEN
				REF1D=((1+X)*(-5+7*X)*(-3+7*X)*(-1+7*X)*(1+7*X)*(3+7*X)*(5+7*X))/92160.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(225+25382*X-38073*X**2-336140*X**3+420175*X**4+705894*X**5-823543*X**6)/92160.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(49*(-63-4990*X+10479*X**2+57820*X**3-101185*X**4-72030*X**5+117649*X**6))/92160.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(-49*(-175-7794*X+27279*X**2+44100*X**3-128625*X**4-43218*X**5+117649*X**6))/30720.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(49*(-1575-3782*X+39711*X**2+16268*X**3-142345*X**4-14406*X**5+117649*X**6))/18432.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(-49*(-1575+3782*X+39711*X**2-16268*X**3-142345*X**4+14406*X**5+117649*X**6))/18432.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(49*(-175+7794*X+27279*X**2-44100*X**3-128625*X**4+43218*X**5+117649*X**6))/30720.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(-49*(-63+4990*X+10479*X**2-57820*X**3-101185*X**4+72030*X**5+117649*X**6))/92160.0D0
			ELSEIF(INDX.EQ.7)THEN
				REF1D=(-225+25382*X+38073*X**2-336140*X**3-420175*X**4+705894*X**5+823543*X**6)/92160.0D0
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.8)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D =((-1+X)*X*(-1+2*X)*(1+2*X)*(-3+4*X)*(-1+4*X)*(1+4*X)*(3+4*X))/630.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-16*(-1+X)*X*(1+X)*(-1+2*X)*(1+2*X)*(-3+4*X)*(-1+4*X)*(1+4*X))/315.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(4*(-1+X)*X*(1+X)*(-1+2*X)*(-3+4*X)*(-1+4*X)*(1+4*X)*(3+4*X))/45.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-16*(-1+X)*X*(1+X)*(-1+2*X)*(1+2*X)*(-3+4*X)*(-1+4*X)*(3+4*X))/45.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(9.0D0-205*X**2+1092*X**4-1920*X**6+1024*X**8)/9.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(-16*(-1+X)*X*(1+X)*(-1+2*X)*(1+2*X)*(-3+4*X)*(1+4*X)*(3+4*X))/45.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(4*(-1+X)*X*(1+X)*(1+2*X)*(-3+4*X)*(-1+4*X)*(1+4*X)*(3+4*X))/45.0D0
			ELSEIF(INDX.EQ.7)THEN
				REF1D=(-16*(-1+X)*X*(1+X)*(-1+2*X)*(1+2*X)*(-1+4*X)*(1+4*X)*(3+4*X))/315.0D0
			ELSEIF(INDX.EQ.8)THEN
				REF1D=(X*(1+X)*(-1+2*X)*(1+2*X)*(-3+4*X)*(-1+4*X)*(1+4*X)*(3+4*X))/630.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(9-18*X-588*X**2+784*X**3+4480*X**4-5376*X**5-7168*X**6+8192*X**7)/630.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-16*(3-8*X-189*X**2+336*X**3+1260*X**4-2016*X**5-1344*X**6+2048*X**7))/315.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(4*(9-36*X-507*X**2+1352*X**3+2080*X**4-4992*X**5-1792*X**6+4096*X**7))/45.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-16*(9-72*X-183*X**2+976*X**3+580*X**4-2784*X**5-448*X**6+2048*X**7))/45.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(2*X*(-205+2184*X**2-5760*X**4+4096*X**6))/9.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(-16*(-9-72*X+183*X**2+976*X**3-580*X**4-2784*X**5+448*X**6+2048*X**7))/45.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(4*(-9-36*X+507*X**2+1352*X**3-2080*X**4-4992*X**5+1792*X**6+4096*X**7))/45.0D0
			ELSEIF(INDX.EQ.7)THEN
				REF1D=(-16*(-3-8*X+189*X**2+336*X**3-1260*X**4-2016*X**5+1344*X**6+2048*X**7))/315.0D0
			ELSEIF(INDX.EQ.8)THEN
				REF1D=(-9-18*X+588*X**2+784*X**3-4480*X**4-5376*X**5+7168*X**6+8192*X**7)/630.0D0
			ENDIF
		ELSEIF(DIFF.EQ.2) THEN
			IF(INDX.EQ.0) THEN
				REF1D = (-9d0 - 588*x + 1176*x**2 + 8960*x**3 - 13440*x**4 - 21504*x**5 + 28672*x**6)/315.d0
			ELSEIF(INDX.EQ.1)THEN
				REF1D = (-32*(-4d0 - 189*x + 504*x**2 + 2520*x**3 - 5040*x**4 - 4032*x**5 + 7168*x**6))/315.d0
			ELSEIF(INDX.EQ.2)THEN
				REF1D = (8*(-18 - 507*x + 2028*x**2 + 4160*x**3 - 12480*x**4 - 5376*x**5 + 14336*x**6))/45.d0
			ELSEIF(INDX.EQ.3)THEN
				REF1D = (-32*(-36d0 - 183*x + 1464*x**2 + 1160*x**3 - 6960*x**4 - 1344*x**5 + 7168*x**6))/45.d0
			ELSEIF(INDX.EQ.4)THEN
				REF1D = -45.55555555555555550d0 + 1456*x**2 - 6400*x**4 + (57344*x**6)/9.d0
			ELSEIF(INDX.EQ.5)THEN
				REF1D = (-32*(-36d0 + 183*x + 1464*x**2 - 1160*x**3 - 6960*x**4 + 1344*x**5 + 7168*x**6))/45.d0
			ELSEIF(INDX.EQ.6)THEN
				REF1D = (8*(-18d0 + 507*x + 2028*x**2 - 4160*x**3 - 12480*x**4 + 5376*x**5 + 14336*x**6))/45.d0
			ELSEIF(INDX.EQ.7)THEN
				REF1D = (-32*(-4d0 + 189*x + 504*x**2 - 2520*x**3 - 5040*x**4 + 4032*x**5 + 7168*x**6))/315.d0
			ELSEIF(INDX.EQ.8)THEN
				REF1D = (-9d0 + 588*x + 1176*x**2 - 8960*x**3 - 13440*x**4 + 21504*x**5 + 28672*x**6)/315.d0
			ENDIF
		ENDIF
	ELSEIF(IORDER.EQ.10)THEN
		IF(DIFF.EQ.0)THEN
			IF(INDX.EQ.0)THEN
				REF1D=((-1+X)*X*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/145152.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-25*(-1+X)*X*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)&
				*(1+5*X)*(2+5*X)*(3+5*X))/72576.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(25*(-1+X)*X*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)&
				*(1+5*X)*(2+5*X)*(4+5*X))/16128.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-25*(-1+X)*X*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)&
				*(1+5*X)*(3+5*X)*(4+5*X))/6048.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(25*(-1+X)*X*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/3456.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=-((-1+X)*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/576.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(25*(-1+X)*X*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/3456.0D0
			ELSEIF(INDX.EQ.7)THEN
				REF1D=(-25*(-1+X)*X*(1+X)*(-4+5*X)*(-3+5*X)*(-1+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/6048.0D0
			ELSEIF(INDX.EQ.8)THEN
				REF1D=(25*(-1+X)*X*(1+X)*(-4+5*X)*(-2+5*X)*(-1+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/16128.0D0
			ELSEIF(INDX.EQ.9)THEN
				REF1D=(-25*(-1+X)*X*(1+X)*(-3+5*X)*(-2+5*X)*(-1+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/72576.0D0
			ELSEIF(INDX.EQ.10)THEN
				REF1D=(X*(1+X)*(-4+5*X)*(-3+5*X)*(-2+5*X)*(-1+5*X)*(1+5*X)&
				*(2+5*X)*(3+5*X)*(4+5*X))/145152.0D0
			ENDIF
		ELSEIF(DIFF.EQ.1)THEN
			IF(INDX.EQ.0)THEN
				REF1D=(-576+1152*X+61500*X**2-82000*X**3-853125*X**4+1023750*X**5 &
				+3281250*X**6-3750000*X**7-3515625*X**8+3906250*X**9)/145152.0D0
			ELSEIF(INDX.EQ.1)THEN
				REF1D=(-25*(-72+180*X+7566*X**2-12610*X**3-99750*X**4+149625*X**5&
				+341250*X**6-487500*X**7-281250*X**8+390625*X**9))/36288.0D0
			ELSEIF(INDX.EQ.2)THEN
				REF1D=(25*(-192+640*X+19476*X**2-43280*X**3-228375*X**4+456750*X**5&
				+603750*X**6-1150000*X**7-421875*X**8+781250*X**9))/16128.0D0
			ELSEIF(INDX.EQ.3)THEN
				REF1D=(-25*(-144+720*X+13107*X**2-43690*X**3-102375*X**4+307125*X**5&
				+223125*X**6-637500*X**7-140625*X**8+390625*X**9))/3024.0D0
			ELSEIF(INDX.EQ.4)THEN
				REF1D=(25*(-576+5760*X+20028*X**2-133520*X**3-121125*X**4+726750*X**5&
				+236250*X**6-1350000*X**7 &
				-140625*X**8+781250*X**9))/3456.0D0
			ELSEIF(INDX.EQ.5)THEN
				REF1D=(-21076*X+382250*X**3-1918125*X**5+3437500*X**7-1953125*X**9)/288.0D0
			ELSEIF(INDX.EQ.6)THEN
				REF1D=(25*(576+5760*X-20028*X**2-133520*X**3+121125*X**4+726750*X**5&
				-236250*X**6-1350000*X**7 &
				+140625*X**8+781250*X**9))/3456.0D0
			ELSEIF(INDX.EQ.7)THEN
				REF1D=(-25*(144+720*X-13107*X**2-43690*X**3+102375*X**4+307125*X**5&
				-223125*X**6-637500*X**7 &
				+140625*X**8+390625*X**9))/3024.0D0
			ELSEIF(INDX.EQ.8)THEN
				REF1D=(25*(192+640*X-19476*X**2-43280*X**3+228375*X**4+456750*X**5&
				-603750*X**6-1150000*X**7 &
				+421875*X**8+781250*X**9))/16128.0D0
			ELSEIF(INDX.EQ.9)THEN
				REF1D=(-25*(72+180*X-7566*X**2-12610*X**3+99750*X**4+149625*X**5&
				-341250*X**6-487500*X**7 &
				+281250*X**8+390625*X**9))/36288.0D0
			ELSEIF(INDX.EQ.10)THEN
				REF1D=(576+1152*X-61500*X**2-82000*X**3+853125*X**4+1023750*X**5&
				-3281250*X**6-3750000*X**7 &
				+3515625*X**8+3906250*X**9)/145152.0D0
			ENDIF
		ELSEIF(DIFF.EQ.2) THEN
			IF(INDX.EQ.0) THEN
				REF1D = (192.d0 + 20500*x - 41000*x**2 - 568750*x**3 + 853125*x**4 + 3281250*x**5 - &
									4375000*x**6 - 4687500*x**7 + 5859375*x**8)/24192.d0
			ELSEIF(INDX.EQ.1) THEN
				REF1D = (-25.d0*(60.d0 + 5044*x - 12610*x**2 - 133000*x**3 + 249375*x**4 + 682500*x**5 - &
									1137500*x**6 - 750000*x**7 + 1171875*x**8))/12096.d0
			ELSEIF(INDX.EQ.2) THEN
				REF1D = (25.d0*(320.d0 + 19476*x - 64920*x**2 - 456750*x**3 + 1141875*x**4 + 1811250*x**5 - &
									4025000*x**6 - 1687500*x**7 + 3515625*x**8))/8064.d0
			ELSEIF(INDX.EQ.3) THEN
				REF1D = (-25.d0*(240.d0 + 8738*x - 43690*x**2 - 136500*x**3 + 511875*x**4 + 446250*x**5 - &
									1487500*x**6 - 375000*x**7 + 1171875*x**8))/1008.d0
			ELSEIF(INDX.EQ.4) THEN
				REF1D = (25.d0*(960.d0 + 6676*x - 66760*x**2 - 80750*x**3 + 605625*x**4 + 236250*x**5 - &
									1575000*x**6 - 187500*x**7 + 1171875*x**8))/576.d0
			ELSEIF(INDX.EQ.5) THEN
				REF1D = (-21076.d0 + 1146750*x**2 - 9590625*x**4 + 24062500*x**6 - 17578125*x**8)/288.d0
			ELSEIF(INDX.EQ.6) THEN
				REF1D = (25.d0*(960.d0 - 6676*x - 66760*x**2 + 80750*x**3 + 605625*x**4 - 236250*x**5 - &
									1575000*x**6 + 187500*x**7 + 1171875*x**8))/576.d0
			ELSEIF(INDX.EQ.7) THEN
				REF1D = (-25.d0*(240.d0 - 8738*x - 43690*x**2 + 136500*x**3 + 511875*x**4 - 446250*x**5 - &
									1487500*x**6 + 375000*x**7 + 1171875*x**8))/1008.d0
			ELSEIF(INDX.EQ.8) THEN
				REF1D = (25.d0*(320.d0 - 19476*x - 64920*x**2 + 456750*x**3 + 1141875*x**4 - 1811250*x**5 - &
									4025000*x**6 + 1687500*x**7 + 3515625*x**8))/8064.d0
			ELSEIF(INDX.EQ.9) THEN
				REF1D = (-25.d0*(60.d0 - 5044*x - 12610*x**2 + 133000*x**3 + 249375*x**4 - 682500*x**5 - &
									1137500*x**6 + 750000*x**7 + 1171875*x**8))/12096.d0
			ELSEIF(INDX.EQ.10) THEN
				REF1D = (192.d0 - 20500*x - 41000*x**2 + 568750*x**3 + 853125*x**4 - 3281250*x**5 - &
									4375000*x**6 + 4687500*x**7 + 5859375*x**8)/24192.d0
			ENDIF
		ENDIF
	ELSEIF (IORDER.EQ.12) THEN
		IF (DIFF.EQ.0) THEN
			IF (INDX.EQ.0) THEN
				REF1D = ((-1 + x)*x*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/92400.d0
			ELSEIF (INDX.EQ.1) THEN
				REF1D = (-3*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*&
									(1 + 3*x)*(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x))/3850.d0
			ELSEIF (INDX.EQ.2) THEN
				REF1D = (3*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*&
									(1 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/1400.d0
			ELSEIF (INDX.EQ.3) THEN
				REF1D = -((-1 + x)*x*(1 + x)*(-1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/210.d0
			ELSEIF (INDX.EQ.4) THEN
				REF1D = (9*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/560.d0
			ELSEIF (INDX.EQ.5) THEN
				REF1D = (-9*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*&
									(1 + 3*x)*(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(5 + 6*x))/175.d0
			ELSEIF (INDX.EQ.6) THEN
				REF1D = ((-1 + x)*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/100.d0
			ELSEIF (INDX.EQ.7) THEN
				REF1D = (-9*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*&
									(1 + 3*x)*(2 + 3*x)*(-5 + 6*x)*(1 + 6*x)*(5 + 6*x))/175.d0
			ELSEIF (INDX.EQ.8) THEN
				REF1D = (9*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/560.d0
			ELSEIF (INDX.EQ.9) THEN
				REF1D = -((-1 + x)*x*(1 + x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/210.d0
			ELSEIF (INDX.EQ.10) THEN
				REF1D = (3*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-1 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/1400.d0
			ELSEIF (INDX.EQ.11) THEN
				REF1D = (-3*(-1 + x)*x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*&
									(1 + 3*x)*(2 + 3*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/3850.d0
			ELSEIF (INDX.EQ.12) THEN
				REF1D = (x*(1 + x)*(-1 + 2*x)*(1 + 2*x)*(-2 + 3*x)*(-1 + 3*x)*(1 + 3*x)*&
									(2 + 3*x)*(-5 + 6*x)*(-1 + 6*x)*(1 + 6*x)*(5 + 6*x))/92400.d0
			ENDIF
		ELSEIF (DIFF.EQ.1) THEN
			IF (INDX.EQ.0) THEN
				REF1D = (100.d0 - 200*x - 15807*x**2 + 21076*x**3 + 344025*x**4 - 412830*x**5 - 2320164*x**6 + &
									2651616*x**7 + 5773680*x**8 - 6415200*x**9 - 4618944*x**10 + 5038848*x**11)/92400.d0
			ELSEIF (INDX.EQ.1) THEN
				REF1D = (-3.d0*(20.d0 - 48*x - 3135*x**2 + 5016*x**3 + 66550*x**4 - 95832*x**5 - 426195*x**6 + &
									584496*x**7 + 962280*x**8 - 1283040*x**9 - 641520*x**10 + 839808*x**11))/3850.d0
			ELSEIF (INDX.EQ.2) THEN
				REF1D = (3.d0*(25.d0 - 75*x - 3858*x**2 + 7716*x**3 + 78125*x**4 - 140625*x**5 - 454356*x**6 + &
									778896*x**7 + 874800*x**8 - 1458000*x**9 - 513216*x**10 + 839808*x**11))/700.d0
			ELSEIF (INDX.EQ.3) THEN
				REF1D = (-100.d0 + 400*x + 14907*x**2 - 39752*x**3 - 270990*x**4 + 650376*x**5 + 1284255*x**6 - &
									2935440*x**7 - 2152008*x**8 + 4782240*x**9 + 1154736*x**10 - 2519424*x**11)/210.d0
			ELSEIF (INDX.EQ.4) THEN
				REF1D = (9.d0*(100.d0 - 600*x - 13407*x**2 + 53628*x**3 + 169265*x**4 - 609354*x**5 - 669060*x**6 + &
									2293920*x**7 + 1014768*x**8 - 3382560*x**9 - 513216*x**10 + 1679616*x**11))/560.d0
			ELSEIF (INDX.EQ.5) THEN
				REF1D = (-9.d0*(100.d0 - 1200*x - 5307*x**2 + 42456*x**3 + 51950*x**4 - 374040*x**5 - 183519*x**6 + &
									1258416*x**7 + 262440*x**8 - 1749600*x**9 - 128304*x**10 + 839808*x**11))/175.d0
			ELSEIF (INDX.EQ.6) THEN
				REF1D = (x*(-5369.d0 + 148148*x**2 - 1200771*x**4 + 3891888*x**6 - 5307120*x**8 + 2519424*x**10))/50.d0
			ELSEIF (INDX.EQ.7) THEN
				REF1D = (-9.d0*(-100.d0 - 1200*x + 5307*x**2 + 42456*x**3 - 51950*x**4 - 374040*x**5 + 183519*x**6 + &
									1258416*x**7 - 262440*x**8 - 1749600*x**9 + 128304*x**10 + 839808*x**11))/175.d0
			ELSEIF (INDX.EQ.8) THEN
				REF1D = (9.d0*(-100.d0 - 600*x + 13407*x**2 + 53628*x**3 - 169265*x**4 - 609354*x**5 + 669060*x**6 + &
									2293920*x**7 - 1014768*x**8 - 3382560*x**9 + 513216*x**10 + 1679616*x**11))/560.d0
			ELSEIF (INDX.EQ.9) THEN
				REF1D = (100.d0 + 400*x - 14907*x**2 - 39752*x**3 + 270990*x**4 + 650376*x**5 - 1284255*x**6 - &
									2935440*x**7 + 2152008*x**8 + 4782240*x**9 - 1154736*x**10 - 2519424*x**11)/210.d0
			ELSEIF (INDX.EQ.10) THEN
				REF1D = (3.d0*(-25.d0 - 75*x + 3858*x**2 + 7716*x**3 - 78125*x**4 - 140625*x**5 + 454356*x**6 + &
									778896*x**7 - 874800*x**8 - 1458000*x**9 + 513216*x**10 + 839808*x**11))/700.d0
			ELSEIF (INDX.EQ.11) THEN
				REF1D = (-3.d0*(-20.d0 - 48*x + 3135*x**2 + 5016*x**3 - 66550*x**4 - 95832*x**5 + 426195*x**6 + &
									584496*x**7 - 962280*x**8 - 1283040*x**9 + 641520*x**10 + 839808*x**11))/3850.d0
			ELSEIF (INDX.EQ.12) THEN
				REF1D = (-100.d0 - 200*x + 15807*x**2 + 21076*x**3 - 344025*x**4 - 412830*x**5 + 2320164*x**6 + &
									2651616*x**7 - 5773680*x**8 - 6415200*x**9 + 4618944*x**10 + 5038848*x**11)/92400.d0
			ENDIF
		ELSEIF (DIFF.EQ.2) THEN
			IF (INDX.EQ.0) THEN
				REF1D = -0.00216450216450216450d0 - (479*x)/1400.d0 + (479*x**2)/700.d0 + (417*x**3)/28.d0 - (1251*x**4)/56.d0 - &
								(7533*x**5)/50.d0 + (5022*x**6)/25.d0 + (17496*x**7)/35.d0 - (4374*x**8)/7.d0 - (17496*x**9)/35.d0 + (104976*x**10)/175.d0
			ELSEIF (INDX.EQ.1) THEN
				REF1D = (-3.d0*(-24.d0 - 3135*x + 7524*x**2 + 133100*x**3 - 239580*x**4 - 1278585*x**5 + 2045736*x**6 + &
								3849120*x**7 - 5773680*x**8 - 3207600*x**9 + 4618944*x**10))/1925.d0
			ELSEIF (INDX.EQ.2) THEN
				REF1D = (3.d0*(-75.d0 - 7716*x + 23148*x**2 + 312500*x**3 - 703125*x**4 - 2726136*x**5 + 5452272*x**6 + &
								6998400*x**7 - 13122000*x**8 - 5132160*x**9 + 9237888*x**10))/700.d0
			ELSEIF (INDX.EQ.3) THEN
				REF1D = (200.d0 + 14907*x - 59628*x**2 - 541980*x**3 + 1625940*x**4 + 3852765*x**5 - 10274040*x**6 - &
								8608032*x**7 + 21520080*x**8 + 5773680*x**9 - 13856832*x**10)/105.d0
			ELSEIF (INDX.EQ.4) THEN
				REF1D = (9.d0*(-300.d0 - 13407*x + 80442*x**2 + 338530*x**3 - 1523385*x**4 - 2007180*x**5 + 8028720*x**6 + &
								4059072*x**7 - 15221520*x**8 - 2566080*x**9 + 9237888*x**10))/280.d0
			ELSEIF (INDX.EQ.5) THEN
				REF1D = (-18.d0*(-600.d0 - 5307*x + 63684*x**2 + 103900*x**3 - 935100*x**4 - 550557*x**5 + 4404456*x**6 + &
								1049760*x**7 - 7873200*x**8 - 641520*x**9 + 4618944*x**10))/175.d0
			ELSEIF (INDX.EQ.6) THEN
				REF1D = (-5369.d0 + 444444.d0*x**2 - 6003855*x**4 + 27243216*x**6 - 47764080*x**8 + 27713664*x**10)/50.d0
			ELSEIF (INDX.EQ.7) THEN
				REF1D = (-18.d0*(-600.d0 + 5307*x + 63684*x**2 - 103900*x**3 - 935100*x**4 + 550557*x**5 + 4404456*x**6 - &
								1049760*x**7 - 7873200*x**8 + 641520*x**9 + 4618944*x**10))/175.d0
			ELSEIF (INDX.EQ.8) THEN
				REF1D = (9.d0*(-300.d0 + 13407*x + 80442*x**2 - 338530*x**3 - 1523385*x**4 + 2007180*x**5 + 8028720*x**6 - &
								4059072*x**7 - 15221520*x**8 + 2566080*x**9 + 9237888*x**10))/280.d0
			ELSEIF (INDX.EQ.9) THEN
				REF1D = (200.d0 - 14907*x - 59628*x**2 + 541980*x**3 + 1625940*x**4 - 3852765*x**5 - 10274040*x**6 + &
								8608032*x**7 + 21520080*x**8 - 5773680*x**9 - 13856832*x**10)/105.d0
			ELSEIF (INDX.EQ.10) THEN
				REF1D = (3.d0*(-75.d0 + 7716*x + 23148*x**2 - 312500*x**3 - 703125*x**4 + 2726136*x**5 + 5452272*x**6 - &
								6998400*x**7 - 13122000*x**8 + 5132160*x**9 + 9237888*x**10))/700.d0
			ELSEIF (INDX.EQ.11) THEN
				REF1D = (-3.d0*(-24.d0 + 3135*x + 7524*x**2 - 133100*x**3 - 239580*x**4 + 1278585*x**5 + 2045736*x**6 - &
								3849120*x**7 - 5773680*x**8 + 3207600*x**9 + 4618944*x**10))/1925.d0
			ELSEIF (INDX.EQ.12) THEN
				REF1D = -0.00216450216450216450d0 + (479*x)/1400.d0 + (479*x**2)/700.d0 - (417*x**3)/28.d0 - (1251*x**4)/56.d0 + &
								(7533*x**5)/50.d0 + (5022*x**6)/25.d0 - (17496*x**7)/35.d0 - (4374*x**8)/7.d0 + &
								(17496*x**9)/35.d0 + (104976*x**10)/175.d0
			ENDIF
		ENDIF
! 	ELSEIF (IORDER.EQ.14) THEN
! 		IF (DIFF.EQ.0) THEN
! 			IF (INDX.EQ.0) THEN
! 				REF1D = ((-1 + x)*x*(-6 + 7*x)*(-5 + 7*x)*(-4 + 7*x)*(-3 + 7*x)*(-2 + 7*x)*(-1 + 7*x)*(1 + 7*x)*(2 + 7*x)*(3 + 7*x)*(4 + 7*x)*(5 + 7*x)*(6 + 7*x))/1.7791488e9
! 			ELSEIF (INDX.EQ.1) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.2) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.3) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.4) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.5) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.6) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.7) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.8) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.9) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.10) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.11) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.12) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.13) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.14) THEN
! 				REF1D = 
! 			ENDIF
! 		ELSEIF (DIFF.EQ.1) THEN
! 			IF (INDX.EQ.0) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.1) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.2) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.3) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.4) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.5) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.6) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.7) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.8) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.9) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.10) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.11) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.12) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.13) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.14) THEN
! 				REF1D = 
! 			ENDIF
! 		ELSEIF (DIFF.EQ.2) THEN
! 			IF (INDX.EQ.0) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.1) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.2) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.3) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.4) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.5) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.6) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.7) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.8) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.9) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.10) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.11) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.12) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.13) THEN
! 				REF1D = 
! 			ELSEIF (INDX.EQ.14) THEN
! 				REF1D = 
! 			ENDIF
! 		ENDIF
	ENDIF

END FUNCTION REF1D

END MODULE INTERPOLATION
