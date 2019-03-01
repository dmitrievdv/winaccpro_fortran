
cd in_data
MODEL_FILES=`cat chosen`
cd ..
# echo $MODEL_FILES
# exit
# pwd
if [ ! -d "logs" ]; then
	mkdir logs
fi
for MODEL_FILE in $MODEL_FILES; do
	# echo $MODEL_FILE

	MODEL=`echo $MODEL_FILE | rev | cut -d '_' -f2- | rev`
	# if [ ! -d "$MODEL""_nohotspot" ]; then
	if [ ! -d  "$MODEL" ]; then
		echo $MODEL 
		./script.sh $MODEL hotspot
	fi
	# fi
done
# ./script.sh 