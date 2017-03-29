#!/bin/bash

#Modify the line below to specify the location of the Ash3d directory,
#and the wind file directory
ASH3DBIN="./"
ASH3DDATA="/Users/yuhsuan/WorkDir/ash3d"
WINDFILEDIR="/Users/yuhsuan/WorkDir/windfiles"
folder="MSH_1980_cloud"                             #name of zip file

echo "removing old .dat, .kmz, .zip, and .txt files"
rm -f *.kmz *.dat AshArrivalTimes.txt *.zip *.gif *.nc *.gif


#The graphic files copied below are legends used in the kmz files
# echo "copying legend graphics"
# cp ${ASH3DDATA}/src/autorun_scripts/USGS_warning3.png .
# cp ${ASH3DDATA}/src/autorun_scripts/concentration_legend.png .
# cp ${ASH3DDATA}/src/autorun_scripts/CloudHeight_hsv.jpg .
# cp ${ASH3DDATA}/src/autorun_scripts/CloudLoad_hsv.png .
# cp ${ASH3DDATA}/src/autorun_scripts/cloud_arrival_time.png .

echo "copying airports file"
cp ${ASH3DDATA}/src/input_files/GlobalAirports_ewert.txt .

echo "creating soft link to wind file directory"
ln -s ${WINDFILEDIR} Wind_nc

echo "running Ash3d model"
${ASH3DBIN}/cloud ash3d_input.inp

#check for errors
rc=$((rc + $?))
echo "rc=$rc"
if [[ "$rc" -gt 0 ]] ; then
   echo "error.  Exiting."
   exit 1
fi

# echo "zipping up kml files"
# zip cloud_arrivaltimes_airports.kmz AshArrivalTimes.kml    USGS_warning3.png
# zip cloud_arrivaltimes_hours.kmz    CloudArrivalTime.kml   USGS_warning3.png cloud_arrival_time.png
# zip CloudConcentration.kmz          CloudConcentration.kml USGS_warning3.png concentration_legend.png
# zip CloudHeight.kmz                 CloudHeight.kml        USGS_warning3.png CloudHeight_hsv.jpg
# zip CloudLoad.kmz                   CloudLoad.kml          USGS_warning3.png CloudLoad_hsv.png

#Uncomment these lines if you want to create a gif map of the output.  You will need GMT and 
#other utilities to run the script GFSVolc_to_gif_tvar.sh
#cp ${ASH3DBIN}/src/autorun_scripts/Ash3d_*.cpt .
#ln -s ${ASH3DBIN}/src/input_files/world_cities.txt .
#${ASH3DBIN}/src/autorun_scripts/GFSVolc_to_gif_tvar.sh 3
#rm -f current_time.txt *.xy Temp.epsi world_cities.txt *.cpt caption.txt

#check for errors
rc=$?
echo "errors=$rc"
if [[ "$rc" -gt 0 ]] ; then
   echo "error.  Exiting."
   exit 1
fi

echo "removing extra files"
rm -f *.kml *.cpt *.grd *.png *.jpg Wind_nc GlobalAirports_ewert.txt

#zip up results
zip ${folder}.zip *.kmz *.dat Ash3d.lst ash3d_input.inp \
                  AshArrivalTimes.txt ${folder}

#write out time of simulation
time2=`date`
echo "${folder}   start: ${time1}  end: ${time2}"
echo "Done with ${folder}."

echo "All done"
