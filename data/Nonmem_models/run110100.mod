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
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0
ETVMAX  = ETA(1) + iETA1 ; VMAX
ETVC  = ETA(2) + iETA2 ; VC
ETKA  = ETA(3) + iETA3 ; KA
ETKM  = ETA(4) + iETA4 ; KM

TVVMAX = THETA(1) ;10000
TVVC   = THETA(2) ;70.0
TVKA   = THETA(3) ;1.0
TVKM   = THETA(4) ;2500

VMAX  = TVVMAX* EXP(ETVMAX)
VC    = TVVC  * EXP(ETVC)
KA    = TVKA  * EXP(ETKA)
KM    = TVKM * EXP(ETKM)

S2 = VC/1000 ; dv in mg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax 
  COM(2)=-1 ; holder of Tmax
ENDIF

A_0(3) = 0

$DES
DADT(1) = -KA * A(1)
DADT(2) =  KA * A(1) - VMAX*(A(2))/(KM+(A(2)/S2))
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

IPRED = A(2)/S2

DEL=0
IF (IPRED.EQ.0) DEL=1
W=(IPRED+DEL)*THETA(5)
Y=IPRED+W*EPS(1)

IRES=DV-IPRED
IWRES=IRES/W



$THETA
10  ; VMAX
70     ; VC
1     ; KA
3500  ; KM
0.22 ; err prop

$OMEGA STANDARD CORRELATION BLOCK(4)
0.45 ; VMAX
0 0.45 ; VC
0 0 0.45 ; KA 
0 0 0 0.45 ; KM

$SIGMA 1 FIX   

;$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1
$EST METHOD=1 INTER MAXEVAL=0 NOABORT SIG=3 PRINT=1 POSTHOC


$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED VMAX VC KM KA ETA1 ETA2 ETA3 ETA4 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run201.tab FORMAT=,

