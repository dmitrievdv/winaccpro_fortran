CURRDIR=`pwd`
PROGDIR="/home/denis/progs/fortran/profile/bin"
cd $PROGDIR
DATADIR=$1
HOTSPOT=$2

echo "Model: $DATADIR"
echo "Model: $DATADIR" > script.log


if [ -d $DATADIR ]; then
  rm -rf $DATADIR
fi

if [ -z $DATADIR ]; then
  echo 'Empty DATADIR argument'
  echo 'Empty DATADIR argument' >> script.log
  exit 1
fi

mkdir $DATADIR
cd $DATADIR
mkdir png
mkdir eps
mkdir raw_data
cd $PROGDIR
cp plot.gnu $DATADIR


if [ -f "$PROGDIR"/in_data/"$DATADIR"_popul.dat ]; then
  cp "$PROGDIR"/in_data/"$DATADIR"_popul.dat "$PROGDIR"/in_data/level_data
elif [ "$2" != "nomodel" ]; then
  echo "Can't find $DATADIR model"
  echo "Can't find $DATADIR model" >> script.log
  # echo "Can't find $DATADIR model" | mail -s 'errors' dmitrievdv242@gmail.com
  exit 2
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
      call="MAGPROF :accr 3.0 :axis $i 0 0 :dipl 2.2 3 :rad 3 :acc 200 100 200 :star 4000 2 0.8 0 :temp 8000 mag :grid pol :line 2d $levels"
    else
      call="MAGPROF :accr 3.0 :axis $i 0 0 :dipl 2.2 3 :rad 3 :acc 200 100 200 :star 4000 2 0.8 0 :temp 8000 mag :grid pol :line 2d $levels +hotspot -$MDOT"
    fi
    echo $call
    echo $call >> script.log
    ./$call >> logs/"$DATADIR".log

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
# i=10
# line="Hb"
# levels=${!line}

# call="MAGPROF :accr 3.0 :axis $i 0 0 :dipl 2.2 3 :rad 3 :acc 200 100 200 :star 4000 2 0.8 0 :temp 8000 mag :grid dec :line 2d $levels"
# echo $call
# ./$call

# mv "$PROGDIR"/output/prof.dat "$DATADIR"/raw_data/"$i"_"$line"_prof.dat
# mv "$PROGDIR"/output/emm.dat "$DATADIR"/raw_data/"$i"_"$line"_emm.dat
# mv "$PROGDIR"/output/absorb.dat "$DATADIR"/raw_data/"$i"_"$line"_abs.dat

cd $DATADIR
gnuplot plot.gnu
cd $PROGDIR
# zip -r $DATADIR"_png.zip" $DATADIR/png
# zip -r $DATADIR"_eps.zip" $DATADIR/eps

# echo "" | mail -a $DATADIR"_png.zip" -a $DATADIR"_eps.zip" -s "results" dmitrievdv242@gmail.com
# echo "" | mail -a $DATADIR"_png.zip" -a $DATADIR"_eps.zip" -s "results" vgcrao@mail.ru 
# rm -rf *.zip



