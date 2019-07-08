#!/usr/bin/env/ python3
#!/bin/bash

#RECORDS TIME TAKEN
# starttime='date +%s'

#GENERATES RECOMBINATION ABUNDANCE
cd "./Code/eucode_recomb"
rm *.o lirec
make
printf "\nModeling Recombination...\n"
./lirec
cd ..
cp ./eucode_recomb/halfabund.data ./eucode_background/initabundance.data
cp ./eucode_recomb/halfabund.data ./eucode_star/initabundance.data
cp ./eucode_recomb/halfabund.data ./eucode_grid/initabundance.data
printf "done\n\n"

#DISTRIBUTES MODEL.DATA
cd ..
cp ./Code/model.data ./Code/model.grid.data
cp ./Code/model.data ./Code/eucode_background/model.data
cp ./Code/model.data ./Code/eucode_star/model.data
cp ./Code/model.data ./Code/eucode_grid/model.data

#GENERATES BACKGROUND REIONIZATION ABUNDANCE
cd "./Code/eucode_background"
rm *.o lirec
make
printf "\nModeling Reionization Background...\n"
./lirec
cd ..
printf "done\n\n"


#REMOVES OLD FILES
# cd "./eucode_grid/"
# rm zrd_*
rm ./eucode_grid/zrd_*
# cd ..
cd ..
if [ -e "./Results" ]; then
  rm -rf ./Results
fi


#PRODUCES STARS
printf "\nCreating stars...\n"
python ./Code/stars.py
rm ./Code/model.grid.data
printf "done\n\n"

#RECOMPILE EUCODE_STAR
cd "./Code/eucode_star"
rm *.o lirec
make
cd ..
cd ..

#MODELS STROMGREN SPHERES
# #UBUNTU CODE
#declare -a stars
# readarray stars < ./Results/Stars/stardata.txt

#MAC OSX CODE
declare -a stars
let i=0
while IFS=$'\n' read -r linedata; do
    stars[i]="${linedata}"
    ((++i))
done < ./Results/Stars/stardata.txt

printf ${stars[@]}

for i in $(seq 0 $((${#stars[@]}-1))); do
  printf "\n\nModeling star %s: \n" "$((${i}+1))"

  #LOADS STAR DATA
  star="${stars[$i]}"

  #MODELS STROMGREN SPHERE
  cd "./Code/eucode_star"
  printf "%s\n" "${star[@]}" > stardata.txt
  tr -s ' '  '\n' < stardata.txt > star.data
  rm stardata.txt zrd.txt
  ./lirec
  cd ..
  cd ..

  #MANIPULATES DATA
  front=$(<./Code/eucode_star/rfront.txt)
  printf "%s\n" "$front" >> "./Results/Stars/frontdata.txt"
  cp "./Code/eucode_star/zrd.txt" "./Code/eucode_grid/zrd_"`printf "%03d" $i`".txt"
done


#PREPS XSTIFF FILE
cp "./Code/xstiff.grid.f" "./Code/xstiff.f"
cp "./Code/xstiff.dead.f" "./Code/xstiff.d.f"

#DEFINES GRID AND DISTANCES
printf "\nDefining grid...\n"
python ./Code/grid.py
printf "done\n\n"

#UPDATES FILES FOR EUCODE_GRID
# cd "./Code/eucode_grid/"
# rm  stars.data xstiff.f
rm ./Code/eucode_grid/stars.data ./Code/eucode_grid/xstiff.f
# cd ..
# cd ..
cp "./Results/Stars/stardata.txt" "./Code/eucode_grid/stars.data"
mv "./Code/xstiff.f" "./Code/eucode_grid/xstiff.f"
mv "./Code/xstiff.d.f" "./Code/eucode_dead/xstiff.f"


#COUNTS GRID POINTS
cd ./Results/Grid/
num_points=$(ls | wc -l)
cd ..
cd ..

#RECOMPILES EUCODE_GRID
cd "./Code/eucode_grid/"
rm  *.o lirec abund.*
make
cd ..
cd "./eucode_dead/"
rm *.o lirec
make
cd ..
cd ..

#MODELS GRID POINTS
# if ! [ "$(ls ./Code/eucode_grid/sabund)" ]; then
#   mkdir ./Code/eucode_grid/sabund/
#   mkdir ./Code/eucode_grid/model
#   printf "\nYEET\n\nYEET\n"
# fi

for i in $(seq 0 $((${num_points}-2))); do
  printf "\n\nModeling grid point %s: \n\n" "$((${i}+1))"

  #COPIES DISTANCES INTO EUCODE_GRID
  cp "./Results/Grid/point_"`printf "%03d" $i`"/point.data" "./Code/eucode_grid/point.data"
  cp "./Results/Grid/point_"`printf "%03d" $i`"/point.data" "./Code/eucode_dead/point.data"

  #MODELS GRID POINT
  # printf "\nWhile alive: \n"
  cd ./Code/eucode_grid
  ./lirec
  cd ..

  # if ! [ `ls -1A ./eucode_grid/sabund/ | wc -l` -eq 0 ] ; then
  count=`ls -1 ./eucode_grid/sabund/*.txt 2>/dev/null | wc -l`
  if [ $count != 0 ]
  then
    mv "./eucode_grid/sabund/"*".txt" "./eucode_dead/sabund/"
    mv "./eucode_grid/model/"*".txt" "./eucode_dead/model/"
  fi 

  # if [[ -e "./eucode_grid/sabund/000.txt" ]]; then
  #   mv "./eucode_grid/sabund/"*".txt" "./eucode_dead/sabund/"
  #   mv "./eucode_grid/model/"*".txt" "./eucode_dead/model/"
  # fi


  # printf "\n\nWhile dead: \n"
  cd ./eucode_dead
  ./lirec
  cd ..
  cd ..

  #SAVES ABUNDANCE DATA
  for j in $(seq 0 $((${#stars[@]}-1))); do
    if [[ -e "./Code/eucode_dead/sabund/"`printf "%03d" $j`".txt" ]]; then
      printf "%s\n" "${i[@]}" > ./Code/eucode_grid/sloop.data
      printf "%s\n" "${j[@]}" >> ./Code/eucode_grid/sloop.data
      python ./Code/append.py
      rm "./Code/eucode_dead/sabund/"`printf "%03d" $j`".txt"
      rm "./Code/eucode_dead/model/"`printf "%03d" $j`".txt"
      rm "./Code/eucode_grid/abund."`printf "%03d" $j`".txt"
      rm "./Code/eucode_grid/sloop.data"
    fi
  done
  rm "./Code/eucode_grid/point.data"
  rm "./Code/eucode_dead/point.data"
done
printf "done\n\n"


#AVERAGES OVER GRID POINTS
printf "\nAveraging over grid points...\n"
python ./Code/grid.avg.py
printf "done\n"

#PLOTS RESULTS
printf "\nPlotting results...\n"
cp ./Code/plot/combined.plt ./Results/combined.plt
cd ./Results
gnuplot combined.plt
rm combined.plt
printf "done\n"

# cp ./Code/plot/separate.plt ./Results/separate.plt
# cd ./Results
# gnuplot separate.plt
#rm separate.plt

cd ..

#OUTPUTS RUNTIMES
# endtime=$'date +%s'
# runtime=$((endtime-starttime))
#
# printf "\nFinished at: \n"
# printf "%s\n" "$(end)"
#
# printf "\nTime Taken: \n"
# printf "%s" "$(runtime)"

say Done!
