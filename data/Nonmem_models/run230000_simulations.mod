$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING RATE
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=9

$MODEL
COMP=(CENTRAL) ; 1
COMP=(PERIPH) ; 2
COMP=(PERIP) ; 3
COMP=(AUC) ; 4

$ABBREVIATED COMRES = 2

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0

TVCL   =  THETA(1) ; 4.0
TVVC   =  THETA(2) ; 70.0
TVVP   =  THETA(3) ; 50.0
TVQ = THETA(4) ; 4.0
TVV3 = THETA(5) ; 30
TVQ3 = THETA(6) ; 4


IF(ICALL == 4) THEN
	CL    = TVCL  * EXP(ETA(1))
	VC    = TVVC  * EXP(ETA(2))

	VP    = TVVP  * EXP(ETA(3))
	Q     = TVQ

	V3 = TVV3 * EXP(ETA(4))
	Q3 = TVQ3

	K23 = Q / VC
	K32 = Q / VP

	K24 = Q3 / VC
	K42 = Q3 / V3

	K20 = CL / VC

	a0 = K20 * K32 * K42
	a1 = K20 * K42 + K32 * K42 + K32 * K24 + K20 * K32 + K42 * K23
	a2 = K20 + K23 + K24 + K32 + K42

	p1 = a1 - (a2 * a2/3)
	q1 = (2 * a2 * a2 * a2 /27) - (a1 * a2/3) + a0
	r_1 = SQRT(-(p1*p1*p1)/27)
	phi = ACOS((-q1/2)/r_1)/3
	r_2 = 2* EXP(LOG(r_1)/3)

	root1 = -(COS(phi)*r_2-a2/3)
	root2 = -(COS(phi+2*3.141593/3)*r_2-a2/3)
	root3 = -(COS(phi+4*3.141593/3)*r_2-a2/3)

	l1 = root1
	IF(root2.GT.root1.AND.root2.GE.root3) l1=root2
	IF(root3.GT.root2.AND.root3.GT.root1) l1=root3

	l3 = root1
	IF(root2.LT.root1.AND.root2.LE.root3) l3 = root2
	IF(root3.LT.root2.AND.root3.LT.root1) l3 = root3

	HALFLIFE1 = LOG(2)/l1

	HALFLIFE3 = LOG(2)/l3

	
	DOWHILE (HALFLIFE1.GT.15.OR.HALFLIFE3.GT.150.OR.HALFLIFE1.LT.1.OR.HALFLIFE3.LT.6) ; 5x "true" half life
		CALL SIMETA(ETA) 

		CL    = TVCL  * EXP(ETA(1))
		VC    = TVVC  * EXP(ETA(2))

		VP    = TVVP  * EXP(ETA(3))
		Q     = TVQ

		V3 = TVV3 * EXP(ETA(4))
		Q3 = TVQ3

		K23 = Q / VC
		K32 = Q / VP

		K24 = Q3 / VC
		K42 = Q3 / V3

		K20 = CL / VC

		a0 = K20 * K32 * K42
		a1 = K20 * K42 + K32 * K42 + K32 * K24 + K20 * K32 + K42 * K23
		a2 = K20 + K23 + K24 + K32 + K42

		p1 = a1 - (a2 * a2/3)
		q1 = (2 * a2 * a2 * a2 /27) - (a1 * a2/3) + a0
		r_1 = SQRT(-(p1*p1*p1)/27)
		phi = ACOS((-q1/2)/r_1)/3
		r_2 = 2* EXP(LOG(r_1)/3)

		root1 = -(COS(phi)*r_2-a2/3)
		root2 = -(COS(phi+2*3.141593/3)*r_2-a2/3)
		root3 = -(COS(phi+4*3.141593/3)*r_2-a2/3)

		l1 = root1
		IF(root2.GT.root1.AND.root2.GE.root3) l1=root2
		IF(root3.GT.root2.AND.root3.GT.root1) l1=root3

		l3 = root1
		IF(root2.LT.root1.AND.root2.LE.root3) l3 = root2
		IF(root3.LT.root2.AND.root3.LT.root1) l3 = root3

		HALFLIFE1 = LOG(2)/l1
		HALFLIFE3 = LOG(2)/l3

	ENDDO
ENDIF



S1 = VC/1000 ; dv in µg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax 
  COM(2)=-1 ; holder of Tmax 
ENDIF

A_0(4) = 0


$DES
DADT(1) =  - K20*A(1) - K23 * A(1) + K32 * A(2) - K24 * A(1) + K42 * A(3)
DADT(2) = K23 * A(1) - K32 * A(2)
DADT(3) = K24 * A(1) - K42 * A(3)
DADT(4) = A(1)/S1


CC = A(1)/S1

IF (S2_SAMPLING.GE.30) THEN
   TDOSE = 216
ELSE
   TDOSE = 0
ENDIF

IF(CC.GT.COM(1).AND.T.GE.TDOSE) THEN 
   COM(1)=CC 
   COM(2)=T
ENDIF


$ERROR
CMAX	= COM(1)
IF (S2_SAMPLING.GE.30) THEN
   TMAX = COM(2) - 216
ELSE
   TMAX = COM(2) 
ENDIF

AUC = A(4)

IPRED = A(1)/S1

DEL=0
IF (IPRED.EQ.0) DEL=1
W=(IPRED+DEL)*THETA(7)
Y=IPRED+W*ETA(5)

IRES=DV-IPRED
IWRES=IRES/W

$THETA
4 ; CL
70 ; V
50 ; V2
4 ; Q
30 ; V3
4 ; Q3
0.22 ; err prop 

$OMEGA STANDARD CORRELATION BLOCK(4)
0.45 ; CL
0 0.45 ; VC
0 0 0.45 ; VP
0 0 0 0.45 ; V3

$OMEGA 1 FIX
 
$SIGMA
0 FIX    ; prop


$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1


$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VC VP Q Q3 V3 ETA1 ETA2 ETA3 ETA4 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run104.tab FORMAT=,
