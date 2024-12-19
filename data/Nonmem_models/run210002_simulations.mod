$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING RATE
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=9

$MODEL
COMP=(CENTRAL) ; 1
COMP=(AUC) ; 2


$ABBREVIATED COMRES = 2


$PK 


TVCL   = THETA(1) ; 4.0
TVVC   = THETA(2) ; 70.0



IF(ICALL == 4) THEN
	CL    = TVCL  * EXP(ETA(1))
	VC    = TVVC  * EXP(ETA(2) )

	HALFLIFE = LOG(2)/(CL/VC)
	
	DOWHILE (HALFLIFE.GT.60.OR.HALFLIFE.LT.2)
		CALL SIMETA(ETA) 
		CL    = TVCL  * EXP(ETA(1) )
		VC    = TVVC  * EXP(ETA(2) )

		HALFLIFE = LOG(2)/(CL/VC)
	ENDDO
ENDIF


K20 = CL / VC

S1 = VC/1000 ; dv in µg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax 
  COM(2)=-1 ; holder of Tmax
ENDIF

A_0(2) = 0


$DES
DADT(1) = - K20*A(1)
DADT(2) = A(1)/S1

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

AUC = A(2)

DEL=0
IF (F.EQ.0) DEL=1
W=(F+DEL)*THETA(3)
Y=F+W*ETA(3)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W

$THETA
4 ; CL
70 ; VC
0.22 ; err prop 

$OMEGA STANDARD CORRELATION BLOCK(2)
0.77 ; CL
0 0.77 ; VC

$OMEGA 1 FIX
 
$SIGMA
0 FIX    ; prop


$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1


$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VC ETA1 ETA2 ETA3 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run002.tab FORMAT=,

