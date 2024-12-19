$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=8

$MODEL
COMP=(DEPOT) ; 1
COMP=(CENTRAL) ; 2
COMP=(AUC) ; 3

$ABBREVIATED COMRES = 2

$PK 

TVVMAX = THETA(1) ;10000
TVVC   = THETA(2) ;70.0
TVKA   = THETA(3) ;1.0
TVKM   = THETA(4) ;2500
TVCL = THETA(5) ; 4

IF(ICALL == 4) THEN
	VMAX  = TVVMAX* EXP(ETA(1))
	VC    = TVVC  * EXP(ETA(2))
	KA    = TVKA  * EXP(ETA(3))
	KM    = TVKM * EXP(ETA(4))
	CL = TVCL * EXP(ETA(5))
	
	CLEARANCE = CL + (VMAX*1000/KM)
	HALFLIFE = LOG(2)/(CLEARANCE/VC)

	DOWHILE (HALFLIFE.GT.75.OR.HALFLIFE.LT.3)
		CALL SIMETA(ETA) 
		VMAX  = TVVMAX* EXP(ETA(1))
		VC    = TVVC  * EXP(ETA(2))
		KA    = TVKA  * EXP(ETA(3))
		KM    = TVKM * EXP(ETA(4))
		CL = TVCL * EXP(ETA(5))
	
		CLEARANCE = CL + (VMAX*1000/KM)
		HALFLIFE = LOG(2)/(CLEARANCE/VC)

	ENDDO
ENDIF


K20 = CL/VC
S2 = VC/1000 ; dv in mg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax 
  COM(2)=-1 ; holder of Tmax
ENDIF

A_0(3) = 0

$DES
DADT(1) = -KA * A(1)
DADT(2) =  KA * A(1) - VMAX*(A(2))/(KM+(A(2)/S2)) - K20*A(2)
DADT(3) = A(2)/S2

CC = A(2)/S2


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
W=(F+DEL)*THETA(6)
Y=F+W*ETA(6)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W


$THETA
5 ; VMAX
70     ; VC
1      ; KA
3500   ; KM
2 ; CL
0.22 ; err prop

$OMEGA STANDARD CORRELATION BLOCK(5)
0.45 ; VMAX
0 0.45 ; VC
0 0 0.45 ; KA 
0 0 0 0.45 ; KM
0 0 0 0 0.45 ; CL

$OMEGA 1 FIX
 
$SIGMA
0 FIX    ; prop

$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1

$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VMAX VC KM KA ETA1 ETA2 ETA3 ETA4 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run201.tab FORMAT=,

