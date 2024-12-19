$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=8

$MODEL
COMP=(DEPOT) ; 1
COMP=(CENTRAL) ; 2
COMP=(PERIPH) ; 3 
COMP=(AUC) ; 4

$ABBREVIATED COMRES = 2

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0
iETA5 = 0
ETCL  = ETA(1) + iETA1 ; CL
ETVC  = ETA(2) + iETA2 ; VC
ETKA  = ETA(3) + iETA3 ; KA
ETVP  = ETA(4) + iETA4 ; VP
ETQ   = ETA(5) + iETA5 ; Q

TVCL   =  THETA(1) ; 4.0
TVVC   =  THETA(2) ; 70.0
TVKA   =  THETA(3) ; 1.0
TVVP   =  THETA(4) ;50.0
TVQ = THETA(5) ; 4

CL    = TVCL  * EXP(ETCL)
VC    = TVVC  * EXP(ETVC)
KA    = TVKA  * EXP(ETKA)
VP    = TVVP  * EXP(ETVP)
Q     = TVQ  * EXP(ETQ)

K23 = Q / VC
K32 = Q / VP
K20 = CL / VC

S2 = VC/1000 ; dv in µg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax  
  COM(2)=-1 ; holder of Tmax
ENDIF

A_0(4) = 0


$DES
DADT(1) = -KA * A(1)
DADT(2) =  KA * A(1) - K20*A(2) - K23 * A(2) + K32 * A(3)
DADT(3) = K23 * A(2) - K32 * A(3)
DADT(4) = A(2)/S2

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

AUC = A(4)

DEL=0
IF (F.EQ.0) DEL=1
W=(F+DEL) * THETA(6)
Y=F+W*EPS(1)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W

$THETA 
4 ; CL
70 ; VC
1 ; KA
50 ; V2
4 ; Q
0.22 ; err prop

$OMEGA STANDARD CORRELATION BLOCK(5)
0.77 ; CL
0 0.77 ; VC
0 0 0.77 ; KA 
0 0 0 0.77 ; VP
0 0 0 0 0.77 ; Q
 
$SIGMA 1 FIX 


;$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1
$EST METHOD=1 INTER MAXEVAL=0 NOABORT SIG=3 PRINT=1 POSTHOC

$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VC KA VP Q ETA1 ETA2 ETA3 ETA4 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run101.tab FORMAT=,
