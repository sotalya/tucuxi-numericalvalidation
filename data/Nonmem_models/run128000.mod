$PROBLEM EXAMPLINIB1
$INPUT ID TIME EVID AMT CMT SS II DV S2_SAMPLING
$DATA data.csv IGNORE=@
$SUBROUTINES ADVAN13 TOL=8

$MODEL
COMP=(ERLANG1) ; 1
COMP=(CENTRAL) ; 2
COMP=(PERIPH) ; 3 
COMP=(ERLANG2) ; 4
COMP=(ERLANG3) ; 5
COMP=(ERLANG4) ; 6
COMP=(ERLANG5) ; 7
COMP=(ERLANG6) ; 8
COMP=(ERLANG7) ; 9
COMP=(AUC) ; 10

$ABBREVIATED COMRES = 2

$PK 
iETA1 = 0
iETA2 = 0
iETA3 = 0
iETA4 = 0

ETCL  = ETA(1) + iETA1 ; CL
ETVC  = ETA(2) + iETA2 ; VC
ETVP  = ETA(3) + iETA4 ; VP


ETKTR  = ETA(4) + iETA3 ; KTR

TVCL   =  THETA(1) 
TVVC   =  THETA(2) 
TVVP   =  THETA(3) 
TVQ = THETA(4) 


TVKTR = THETA(5) 


CL    = TVCL  * EXP(ETCL)
VC    = TVVC  * EXP(ETVC)
VP    = TVVP  * EXP(ETVP)
Q     = TVQ


KTR    = TVKTR  * EXP(ETKTR)

K23 = Q / VC
K32 = Q / VP
K20 = CL / VC

S2 = VC/1000 ; dv in µg/l ; amt in mg


IF(NEWIND.LE.1) THEN ; assign negative Cmax Tmax for the new subject 
  COM(1)=-1 ; holder of Cmax  
  COM(2)=-1 ; holder of Tmax
ENDIF

A_0(10) = 0


$DES
DADT(1) = -KTR * A(1)
DADT(4) = KTR * A(1) -KTR * A(4)
DADT(5) = KTR * A(4) -KTR * A(5)
DADT(6) = KTR * A(5) -KTR * A(6)
DADT(7) = KTR * A(6) -KTR * A(7)
DADT(8) = KTR * A(7) -KTR * A(8)
DADT(9) = KTR * A(8) -KTR * A(9)

DADT(2) =  KTR * A(9) - K20*A(2) - K23 * A(2) + K32 * A(3)
DADT(3) = K23 * A(2) - K32 * A(3)

DADT(10) = A(2)/S2

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

AUC = A(10)

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
50 ; V2
4 ; Q
5 ; KTR

0.22 ; err prop

$OMEGA STANDARD CORRELATION BLOCK(4)
0.45 ; CL
0 0.45 ; VC
0 0 0.45 ; VP 
0 0 0 0.45 ; KTR
 
$SIGMA 1 FIX 


;$SIM ONLYSIM (12345 NEW) SUBPROBLEM=1
$EST METHOD=1 INTER MAXEVAL=0 NOABORT SIG=3 PRINT=1 POSTHOC

$TABLE ID TIME EVID DV S2_SAMPLING PRED IPRED CL VC VP Q KTR ETA1 ETA2 ETA3 ETA4 AUC CMAX TMAX NOPRINT NOAPPEND ONEHEADER FILE=run101.tab FORMAT=,
