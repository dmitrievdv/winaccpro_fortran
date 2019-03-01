MODEL=$1
LINES=$3
ACCARG=$2

PROGDIR=/home/denis/progs/fortran/profile/bin
CURRDIR=`pwd`

cd $PROGDIR

[ -z `echo $ACCARG | grep -o p` ]; PLANE=$?
[ -z `echo $ACCARG | grep -o z` ]; Z=$?
[ -z `echo $ACCARG | grep -o l` ]; LINE=$?

echo $PLANE,$Z,$LINE
ACCARGCHECK=$(($PLANE+$Z+$LINE))
if [ $ACCARGCHECK -eq 0 ]; then
  echo "No accuracy argument"
  exit
fi

if [ -z $LINES ]; then
  LINES="Ha Hb"
fi

Ha=3\ 2
Hb=4\ 2
Hg=5\ 2
Pa=4\ 3
Pb=5\ 3
Pg=6\ 3
Bra=5\ 4
Brb=6\ 4
Brg=7\ 4

DATADIR=$MODEL

if [ ! -d $DATADIR ]; then
  echo "Nothing to upgrade!"
  exit
fi

cd logs
if [ ! -f $MODEL.log ]; then
  echo "No accuracy data"
  exit
else
  ACCDATA=`head -6 $MODEL.log | tail -1`
  ACCPLANE=`echo $ACCDATA | cut -d ' ' -f1`
  ACCZ=`echo $ACCDATA | cut -d ' ' -f2`
  ACCLINE=`echo $ACCDATA | cut -d ' ' -f3`
  echo $ACCPLANE,$ACCZ,$ACCLINE
fi
cd ..

if [ $PLANE -ne 0 ]; then
  ACCPLANE=$(( $ACCPLANE * 2 ))
fi

if [ $LINE -ne 0 ]; then
  ACCLINE=$(( $ACCLINE * 2 ))
fi

if [ $Z -ne 0 ]; then
  ACCZ=$(( $ACCZ * 2 ))
fi

echo $ACCPLANE,$ACCZ,$ACCLINE

if [ -f "$PROGDIR"/in_data/"$DATADIR"_popul.dat ]; then
  cp "$PROGDIR"/in_data/"$DATADIR"_popul.dat "$PROGDIR"/in_data/level_data
elif [ "$2" != "nomodel" ]; then
  echo "Can't find $DATADIR model"
  echo "Can't find $DATADIR model" >> script.log
  # echo "Can't find $DATADIR model" | mail -s 'errors' dmitrievdv242@gmail.com
  exit 2
fi

HOTSPOT=hotspot

MDOT=`echo $DATADIR | cut -d '_' -f2`
echo $MDOT
SIZE=${#MDOT}

if [ $SIZE == '2' ]; then
  i=$(($SIZE-1))
  echo $i
  MDOT=${MDOT::-1}'.'${MDOT:$i:1}
fi

echo $MDOT

for i in `seq 15 15 75`
do 
  for line in "Hb" "Ha" # "Hg" "Brg"
  do
    levels=${!line}
    if [ "$HOTSPOT" == '' ]; then
      call="MAGPROF :accr 3.0 :axis $i 0 0 :dipl 2.2 3 :rad 3 :acc $ACCPLANE $ACCZ $ACCLINE :star 4000 2 0.8 0 :temp 8000 mag :grid pol :line 2d $levels"
    else
      call="MAGPROF :accr 3.0 :axis $i 0 0 :dipl 2.2 3 :rad 3 :acc $ACCPLANE $ACCZ $ACCLINE :star 4000 2 0.8 0 :temp 8000 mag :grid pol :line 2d $levels +hotspot -$MDOT"
    fi
    echo $call
    echo $call >> script.log
    ./$call > logs/"$DATADIR".log

    # if [ -f "$PROGDIR"/output/prof.dat ]; then
      mv "$PROGDIR"/output/prof.dat "$DATADIR"/raw_data/"$i"_"$line"_prof.dat
    # else 
      # cp "$DATADIR"_nohotspot/raw_data/"$i"_"$line"_prof.dat "$DATADIR"/raw_data
    # fi
    # if [ -f "$PROGDIR"/output/emm.dat ]; then
      mv "$PROGDIR"/output/emm.dat "$DATADIR"/raw_data/"$i"_"$line"_emm.dat
    # else 
      # cp "$DATADIR"_nohotspot/raw_data/"$i"_"$line"_emm.dat "$DATADIR"/raw_data
    # fi
    # if [ -f "$PROGDIR"/output/abs.dat ]; then
      mv "$PROGDIR"/output/absorb.dat "$DATADIR"/raw_data/"$i"_"$line"_abs.dat
    # else 
      # cp "$DATADIR"_nohotspot/raw_data/"$i"_"$line"_abs.dat "$DATADIR"/raw_data
    # fi
    # exit
  done
done 


cd $CURRDIR
