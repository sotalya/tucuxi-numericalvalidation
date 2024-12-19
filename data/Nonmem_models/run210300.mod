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
iETA3 = 0
iETA4 = 0
ETVMAX  = ETA(1) + iETA1 ; VMAX
ETVC  = ETA(2) + iETA2 ; VC
ETKM  = ETA(3) + iETA3 ; KM
ETCL = ETA(4) + iETA4 ; CL

TVVMAX = THETA(1) ;10000
TVVC   = THETA(2) ;70.0
TVKM   = THETA(3) ;2500
TVCL = THETA(4) ; 4

VMAX  = TVVMAX* EXP(ETVMAX)
VC    = TVVC  * EXP(ETVC)
KM    = TVKM * EXP(ETKM)
CL = TVCL * EXP(ETCL)

K20 = CL/VC

S1 = VC/1000 ; dv in mg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax 
  COM(2)=-1 ; holder of Tmax
ENDIF

A_0(2) = 0

$DES
DADT(1) = - VMAX*(A(1))/(KM+(A(1)/S1)) - K20 * A(1)
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
W=(F+DEL) * THETA(5)
Y=F+W*EPS(1)
IPRED=F
IRES=DV-IPRED
IWRES=IRES/W


$THETA
5 ; VMAX
70  ; VC
3500  ; KM
2 ; CL
0.22  ; err prop

$OMEGA STANDARD CORRELATION BLOCK(4)
0.45 ; VMAX
0 0.45 ; VC
0 0 0.45 ; KM 
0 0 0 0.45  ; CL

$SIGMA 1 FIX


$EST METHOD=1 INTER MAXEVAL=0 NOABORT SIG=3 PRINT=1 POSTHOC


$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VMAX VC KM ETA1 ETA2 ETA3 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run202.tab FORMAT=,
