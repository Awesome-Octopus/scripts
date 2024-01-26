#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1500
XMIN	3345.000000
XWID	149.900000
*
*	Documentational Text:
*
TITL	'TEMPO200uM_h20-45mm3'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    08/09/18
TIME    19:58:17
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.342
A1SW    0.015
MWFQ    9.610525e+09
MWPW    0.01
AVGS    10
SPTP    0.04115
RCAG    35
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    1500
RCTC    0.04096
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     ezcw, 1.0

MetalIonsDef       False
OrgRadicalsDef     False

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3420.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.002300
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3345.000 G
SweepDirection     Up
SweepWidth         150.0 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.610525 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        128.8
Power              10.00 mW
PowerAtten         10.00 dB
QValue             2568

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         10
NbScansDone        10
NbScansToDo        10
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       0

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           1.21e-08
ConvTime           41.15 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               35 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             1.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           10.000
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.993287,0,0,0.00240032,0,0,-0.00984044,0,0,8.8985e-05,0,0,-4.75288e-06,0,0,-1.05279e-06,0,0,2.9393e-07,0,0,2.70257e-09,0,0,-7.36288e-10,0
QuadMode           False
ResCenter          125 mm
ResLength          25 mm
Resolution         1500
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          61.72 s
TimeConst          40.96 ms
TimeExp            False
TuneCaps           1

*
************************************************************
