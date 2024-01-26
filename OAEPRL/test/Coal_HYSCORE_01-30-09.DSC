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
XPTS	128
XMIN	-0.000000
XWID	2032.000000
YPTS	128
YMIN	-0.000000
YWID	2032.000000
*
*	Documentational Text:
*
TITL	'Coal_HYSCORE_01-30-09'
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
DATE    01/30/09
TIME    19:27:34
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    
OXS1    
AXS1    
AXS2    
AXS3    
MWPW    1.998e-07
A1CT    0.3482
B0VL    0.3482
A1SW    0.01
MWFQ    9.778462e+09
AVGS    1
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
Power              0.0001998 mW
PowerAtten         60 dB

.DVC     fieldCtrl, 1.0

AllegroMode        False
CenterField        3482.00 G
Delay              0.0 s
FieldFlyback       On
FieldResol         10.00000
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3480.000 G
SweepDirection     Up
SweepWidth         100.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.778462 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        8.00 dB
ELDORAtt           30 dB
FrequencyA         9.500000 GHz
VideoBW            200 MHz
VideoGain          6 dB

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
FTEzDelay0         400 ns
FTEzDelay1         128 ns
FTEzDelay2         300 ns
FTEzDeltaX         16 ns
FTEzDeltaY         16 ns
FTEzELDORStart     9.100000 GHz
FTEzELDORWidth     0.000000 GHz
FTEzENDORDetDel    0 ns
FTEzExtTrgDel      200 ns
FTEzIntgBase       1 ns
FTEzIntgGate       4 ns
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
PlsSPELEXPSlct     HYSCORE
PlsSPELLISTSlct    4-step
PlsSPELPhPrgEx     Normal
QuadDetect         True
ReplaceMode        Off
ShotRepTime        1999.20 us
ShotsPLoop         100
SptProgress        100 %
StochMode          False
SweepsPExp         1
TriggerTimeOut     38 s
XSpecRes           128
YSpecRes           128

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       1

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************
