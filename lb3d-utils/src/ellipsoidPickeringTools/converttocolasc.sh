#!/bin/bash
id=""
nummer="07"
outname="test"
prognr="09"
gr_out_file=$outname$nummer
prog="ma2tdropletCalc"
progstart="./"$prog
prog=$prog$prognr".cc"
#aus="OP"$nummer

#echo $id > id.conf
#echo $outname > outname.conf


#ls colour_adsorb01_t00006400-0199872626.h5
#h5dump colour_absorb01_t00006400-0199872626.h5 > colour_absorb01_t00006400-0199872626.asc
echo "datei auslesen"
#echo "reading file"
for (( i=0; $i <= 9; i++ )) ; do
#for (( i=0; $i <= 2; i++ )) ; do
    col="colour_"$gr_out_file"_t00000"$i"00$id.h5"
    colasc="colour_"$gr_out_file"_t00000"$i"00$id.asc"
    h5dump $col > $colasc
    echo "$nummer: $i"
done
for (( i=10; $i <= 99; i++ )) ; do
    col="colour_"$gr_out_file"_t0000"$i"00$id.h5"
    colasc="colour_"$gr_out_file"_t0000"$i"00$id.asc"
    h5dump $col > $colasc
    echo "$nummer: $i"
done
for (( i=100; $i <= 999; i++ )) ; do
    col="colour_"$gr_out_file"_t000"$i"00$id.h5"
    colasc="colour_"$gr_out_file"_t000"$i"00$id.asc"
    h5dump $col > $colasc
    echo "$nummer: $i"
done
for (( i=1000; $i <= 1000; i++ )) ; do
    col="colour_"$gr_out_file"_t00"$i"00$id.h5"
    colasc="colour_"$gr_out_file"_t00"$i"00$id.asc"
    h5dump $col > $colasc
    echo "$nummer: $i"
done
#g++ -o $prog $progcode
#echo "Berechnung starten"
#$progstart
ls colour*asc
#rm colour*asc
rm id.conf
rm outname.conf
