CURRDIR=`pwd`
PROGDIR="/home/denis/progs/fortran/profile/bin"
cd $PROGDIR
DATADIR=$1"_rot"

if [ -d $DATADIR ]; then
  rm -rf $DATADIR
fi

if [ -z $DATADIR ]; then
  echo 'Empty DATADIR argument' | mail -s 'errors' dmitrievdv242@gmail.com 
  exit 1
fi

mkdir $DATADIR
cd $DATADIR
mkdir png
mkdir eps
mkdir raw_data
cd $PROGDIR
cp plot.gnu $DATADIR

if [ -f "$PROGDIR"/in_data/"$1"_popul.dat ]; then
  cp "$PROGDIR"/in_data/"$1"_popul.dat "$PROGDIR"/in_data/level_data
elif [ "$2" != "nomodel" ]; then
  echo "Can't find $1 model"
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

for i in `seq 15 30 75`
do 
  for line in "Ha" "Hb"
  do
    levels=${!line}
    for psi in `seq 0 10 180`
    do
      # call="MAGPROF :accr 3.0 :axis $i $psi 15 :dipl 2.2 3 :rad 3 :acc 200 100 200 :star 4000 2 0.8 15 :temp 8000 mag :grid dec :line 2d $levels"
      call="MAGPROF :accr 3.0 :axis $i $psi 15 :dipl 2.2 3 :rad 3 :acc 200 100 200 :star 4000 2 0.8 15 :temp 8000 mag :grid pol :line 2d $levels +hotspot -$MDOT"
      echo $call 
      ./$call >> logs/rotlog.log

      mv "$PROGDIR"/output/prof.dat "$DATADIR"/raw_data/"$line"_"$i"_"$psi"_prof.dat
      mv "$PROGDIR"/output/emm.dat "$DATADIR"/raw_data/"$line"_"$i"_"$psi"_emm.dat
      mv "$PROGDIR"/output/absorb.dat "$DATADIR"/raw_data/"$line"_"$i"_"$psi"_abs.dat
    done
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

# cd $DATADIR
# gnuplot plot.gnu
# cd $PROGDIR
# zip -r $DATADIR"_png.zip" $DATADIR/png
# zip -r $DATADIR"_eps.zip" $DATADIR/eps

# echo "rot for $DATADIR done" | mail -s "results" dmitrievdv242@gmail.com
# echo "" | mail -a $DATADIR"_png.zip" -a $DATADIR"_eps.zip" -s "results" vgcrao@mail.ru 
# rm -rf *.zip



