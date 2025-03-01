$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING bodyweight sex RATE
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=9

$MODEL
COMP=(CENTRAL) ; 1
COMP=(PERIPH) ; 2
COMP=(AUC) ; 3

$ABBREVIATED COMRES = 2

$PK 

TVCL   = THETA(1) ; 4.0
TVVC   = THETA(2) ; 70.0
TVVP   =  THETA(3) ;50.0
TVQ = THETA(4) ; 4


IF(ICALL == 4) THEN
	CL    = TVCL  * EXP(ETA(1)) * ((bodyweight/75)**1.5) * (0.75**sex)
	VC    = TVVC  * EXP(ETA(2))

	VP    = TVVP  * EXP(ETA(3))
	Q     = TVQ

	K23 = Q / VC
	K32 = Q / VP
	K20 = CL / VC

	A0 = K20 * K32
	A1 = -(K23 + K32 + K20)

	L1 = (- A1 + SQRT(A1*A1 - 4*A0))/2
	L2 = (- A1 - SQRT(A1*A1-4*A0))/2

	HALFLIFE1 = LOG(2)/L1
	HALFLIFE2 = LOG(2)/L2
	
	DOWHILE (HALFLIFE1.GT.25.OR.HALFLIFE2.GT.125.OR.HALFLIFE1.LT.1.OR.HALFLIFE2.LT.5) ; 5x "true" half life
		CALL SIMETA(ETA) 

		CL    = TVCL  * EXP(ETA(1)) * ((bodyweight/75)**1.5) * (0.75**sex)
		VC    = TVVC  * EXP(ETA(2))

		VP    = TVVP  * EXP(ETA(3))
		Q     = TVQ

		K23 = Q / VC
		K32 = Q / VP
		K20 = CL / VC

		A0 = K20 * K32
		A1 = -(K23 + K32 + K20)

		L1 = (- A1 + SQRT(A1*A1 - 4*A0))/2
		L2 = (- A1 - SQRT(A1*A1-4*A0))/2

		HALFLIFE1 = LOG(2)/L1
		HALFLIFE2 = LOG(2)/L2

	ENDDO
ENDIF


S1 = VC/1000 ; dv in µg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax 
  COM(2)=-1 ; holder of Tmax 
ENDIF

A_0(3) = 0

$DES
DADT(1) = - K20*A(1) - K23 * A(1) + K32 * A(2)
DADT(2) = K23 * A(1) - K32 * A(2)
DADT(3) = A(1)/S1

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

AUC = A(3)

DEL=0
IF (F.EQ.0) DEL=1
W=(F+DEL)*THETA(5)
Y=F+W*ETA(4)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W

$THETA
4 ; CL
70 ; VC
50 ; VP
4 ; Q
0.22 ; err prop 

$OMEGA STANDARD CORRELATION BLOCK(3)
0.45 ; CL
0 0.45 ; VC
0 0 0.45 ; VP 

$OMEGA 1 FIX
 
$SIGMA
0 FIX    ; prop

$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1


$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VC VP Q ETA1 ETA2 ETA3 ETA4 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run102.tab FORMAT=,

