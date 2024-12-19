$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING RATE
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN6 TOL=9

$MODEL
COMP=(CENTRAL) ; 1
COMP=(AUC) ; 2

$ABBREVIATED COMRES = 2


$PK 
iETA1 = 0
iETA2 = 0
ETCL  = ETA(1) + iETA1 ; CL
ETVC  = ETA(2) + iETA2 ; VC

TVCL   =  THETA(1) ; 4.0
TVVC   =  THETA(2) ; 70.0

CL    = TVCL  * EXP(ETCL )
VC    = TVVC  * EXP(ETVC )
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

W = THETA(3)
Y=F+W*EPS(1)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W

$THETA 
4 ; CL
70 ; VC
50 ; err add 

$OMEGA STANDARD CORRELATION BLOCK(2)
0.45 ; CL
0 0.45 ; VC

$SIGMA 1 FIX ; prop


;$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1
$EST METHOD=1 INTER MAXEVAL=0 NOABORT SIGL=9 NSIG=3 PRINT=1 POSTHOC

$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VC ETA1 ETA2 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run404.tab FORMAT=,

