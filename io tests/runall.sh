cd ~/Documents/GitHub/engine
make
begin1=$(date +%s)
for folder in io*; do
	echo $folder
	cd "$folder"
	rm *.bmp > /dev/null 2>&1
	begin=$(date +%s)
	../engine *.ini
	end=$(date +%s)
	echo "$(expr $end - $begin)"
	cd ..
done
end1=$(date +%s)
echo "$(expr $end1 - $begin1)"
#
#for folder in io*; do
#	echo $folder
#	cd "$folder"
#	../compare.sh . 1
#	cd ..
#done
#
