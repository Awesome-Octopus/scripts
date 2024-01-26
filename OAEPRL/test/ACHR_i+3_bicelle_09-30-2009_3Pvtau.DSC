#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	IDX
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	512
XMIN	-0.000000
XWID	4088.000000
YPTS	70
YMIN	-0.000000
YWID	276.000000
*
*	Documentational Text:
*
TITL	'ACHR_i+3_bicelle_09-30-2009_3Pvtau'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'Time'
YNAM	'Time'
IRUNI	''
IIUNI	''
XUNI	'ns'
YUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    10/01/09
TIME    07:53:41
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    
OXS1    
AXS1    
AXS2    
AXS3    
MWPW    0.0006332
A1CT    0.348
B0VL    0.3295186
A1SW    0.1
MWFQ    9.265087e+09
AVGS    4
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
Power              0.6332 mW
PowerAtten         25 dB

.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        3480.00 G
Delay              0.0 s
FieldFlyback       On
FieldResol         10.00000
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3480.000 G
SweepDirection     Up
SweepWidth         1000.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.265087 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        15.00 dB
ELDORAtt           30 dB
FrequencyA         9.500000 GHz
VideoBW            200 MHz
VideoGain          33 dB

.DVC     ftEpr, 1.0

AutoTimeOut        True
AveragesPerScan    100
FTEzB0SettlDel     0 s
FTEzB0Step         0 G
FTEzCurrELDOR      9.100000 GHz
FTEzD11deer        0 ns
FTEzD20deer        0 ns
FTEzD3deer         300 ns
FTEzDXdeer         0 ns
FTEzDelay0         368 ns
FTEzDelay1         120 ns
FTEzDelay2         60 ns
FTEzDeltaX         8 ns
FTEzDeltaY         4 ns
FTEzELDORStart     9.100000 GHz
FTEzELDORWidth     0.000000 GHz
FTEzENDORDetDel    0 ns
FTEzExtTrgDel      200 ns
FTEzIntgBase       1 ns
FTEzIntgGate       72 ns
FTEzMWPiHalf       16 ns
FTEzP20deer        0 ns
FTEzP3             40 ns
FTEzRF1Pi          0 ns
FTEzRF1PlsSep      4 ns
FTEzRF1Pump        0.500 MHz
FTEzRF1Start       0.500 MHz
FTEzRF1Width       0.000 MHz
FTEzSRELDORExp     False
FTEzSXdeer         1024
FTEzSXdeerExp      False
FTEzSleep          0 s
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     3-Pls vs 2-Pls
PlsSPELLISTSlct    4-step
PlsSPELPhPrgEx     Normal
QuadDetect         True
ReplaceMode        Off
ShotRepTime        999.60 us
ShotsPLoop         100
SptProgress        100 %
StochMode          False
SweepsPExp         4
TriggerTimeOut     66 s
XSpecRes           512
YSpecRes           70

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         4
NbScansDone        3
NbScansToDo        4
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
