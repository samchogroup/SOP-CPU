#!/bin/bash

# directories
SOPDIR=/home/anans14/SOP-CPU
INPUTDIR=/home/anans14/Benchmark/70s/INPUT
OUTPUTDIR=/home/anans14/Benchmark/70s/OUTPUT

method=NAIVE

curr=125
 # set up input file for restart
sed "s/CURR/$curr/g" $INPUTDIR/input_first_$method.temp > $INPUTDIR/input_$curr

# run first trajectory
$SOPDIR/sop.x $INPUTDIR/input_$curr > temp.txt
time=`(cat temp.txt | grep Total | awk '{print $5}')`
energy=$(awk '{print $7}' OUTPUT/70s.0.60.2718.$method.$curr.out )

totalTime=$time
echo $curr $time $totalTime $energy > timings$method.txt

prev=$curr
rm INPUT/input_$curr

for ((i=250; i<=2000; i=i*2))
do
    echo Current step is $curr.
    curr=$i
    timeDiff=$((curr-prev))

    # set up input file for restart
    sed "s/START/$prev/g" $INPUTDIR/input_restart_$method.temp | sed "s/UPDATE/$timeDiff/g" | sed "s/STOP/$curr/g" | sed "s/CURR/$curr/g" > $INPUTDIR/input_$curr

    # copy files from previous run for current run
    cp $OUTPUTDIR/70s.0.60.2718.$method.$prev.coords.out             $OUTPUTDIR/70s.0.60.2718.$method.$curr.coords.out
    cp $OUTPUTDIR/70s.0.60.2718.$method.$prev.coords_uncorrected.out $OUTPUTDIR/70s.0.60.2718.$method.$curr.coords_uncorrected.out
    cp $OUTPUTDIR/70s.0.60.2718.$method.$prev.velocs.out             $OUTPUTDIR/70s.0.60.2718.$method.$curr.velocs.out

    # run restart trajectory
    $SOPDIR/sop.x $INPUTDIR/input_$curr > temp.txt
    time=`(cat temp.txt | grep Total | awk '{print $5}')`
    energy=$(awk '{print $7}' OUTPUT/70s.0.60.2718.$method.$curr.out )
    
    totalTime=`(echo $time + $totalTime | bc -l )`
    echo $curr $time $totalTime $energy >> timings$method.txt

    prev=$curr
    rm INPUT/input_$curr
done

for ((i=3000; i<=10000; i=i+1000))
do
    echo Current step is $curr.
    curr=$i
    timeDiff=$((curr-prev))

    # set up input file for restart
    sed "s/START/$prev/g" $INPUTDIR/input_restart_$method.temp | sed "s/UPDATE/$timeDiff/g" | sed "s/STOP/$curr/g" | sed "s/CURR/$curr/g" > $INPUTDIR/input_$method_$curr

    # copy files from previous run for current run
    cp $OUTPUTDIR/70s.0.60.2718.$method.$prev.coords.out             $OUTPUTDIR/70s.0.60.2718.$method.$curr.coords.out
    cp $OUTPUTDIR/70s.0.60.2718.$method.$prev.coords_uncorrected.out $OUTPUTDIR/70s.0.60.2718.$method.$curr.coords_uncorrected.out
    cp $OUTPUTDIR/70s.0.60.2718.$method.$prev.velocs.out             $OUTPUTDIR/70s.0.60.2718.$method.$curr.velocs.out

    # run restart trajectory
    $SOPDIR/sop.x $INPUTDIR/input_$curr > temp.txt
    time=`(cat temp.txt | grep Total | awk '{print $5}')`
    energy=$(awk '{print $7}' OUTPUT/70s.0.60.2718.$method.$curr.out )
    
    totalTime=`(echo $time + $totalTime | bc -l )`
    echo $curr $time $totalTime $energy >> timings$method.txt

    prev=$curr
    rm INPUT/input_$curr
done

rm temp.txt

