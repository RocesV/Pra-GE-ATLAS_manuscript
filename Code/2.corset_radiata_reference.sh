#!bin/bash

CORSET=/home/rocesv/corset-1.09-linux64/
INPUT=$(ls /mnt/g/Pra-GE-ATLAS/ATLAS/Pinus_radiata/1.Assemblies/3.Corset/Corset_okayaltokayReduced/0.mapping_bams/salmon/mapping/*/*/*.out/aux_info/eq_classes.txt)

$CORSET/corset -D 9999999 -g 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10 \
 -n BU1,BU2,BU3,BU4,BU5,BU6,BU7,BU8,BU9,BU10,BU11,BU12,BU13,BU14,BU15,BU16,ME1,NE1,NE2,NE3,NE4,NE5,NE6,NE7,NE8,NE9,NE10,NE11,NE12,NE13,NE14,NE15,NE16,NE17,NE18,NE19,NE20,NE21,NE22,NE23,NE24,NE25,NE26,NE27,NE28,NE29,PH1,PH2,XY1,XY2,XY3,XY4,XY5,XY6,HE1,HE2,HE3,HE4,HE5,HE6,HE7,HE8,DO1,DO2,DO3,DO4,DO5,DO6,DO7,DO8,DO9,DO10,PY1,PY2,PY3,PY4,PY5,PY6,FU1,FU2,FU3,FU4,FU5,FU6,FU7,FU8,FU9,FU10,FU11,FU12,FU13,FU14,FU15,FU16,FU17,FU18,FU19,FU20,FU21,FU22,FU23,FU24,FU25,FU26,FU27,FU28,FU29,FU30,FU31,FU32,FU33,FU34,FU35,FU36,FU37,FU38,FU39,FU40,FU41,FU42,FU43,FU44,FU45,FU46,FU47,FU48,BD1,BD2,BD3,BD4,BD5,BD6,BD7,BD8,BD9,BD10,BD11,BD12,BD13,BD14,BD15,BD16,BD17,BD18,BD19,BD20,BD21 \
 -i salmon_eq_classes $INPUT
