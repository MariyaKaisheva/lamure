#echo "Test $# $1 $2"

#!/bin/bash

if [ $# -le "1" ]; 
    then echo "illegal number of parameters"
    exit -1
fi

depth=10
input_folder=${1}
input_rot_mat_folder=${2}

#for xyz_file in ${input_folder}*.xyz 
#do
#  ./lamure/install/bin/lamure_preprocessing ${xyz_file}
#done

cd ../install/bin

LOG_FILE_NAME=test_log

touch ${LOG_FILE_NAME}

for bvh_file in `ls -Sr ${input_folder}*.bvh`;
do
	bvh_filename_wo_extensions=$(basename "$bvh_file")

	mkdir ${input_folder}processed/${bvh_filename_wo_extensions}
	mkdir ${input_folder}processed/${bvh_filename_wo_extensions}/point_clouds
	mkdir ${input_folder}processed/${bvh_filename_wo_extensions}/line_objs
	mkdir ${input_folder}processed/${bvh_filename_wo_extensions}/LOGS

	for rot_file in ${input_rot_mat_folder}*.rot
	do

	  rot_filename_wo_extensions=$(basename "$rot_file")

	  #empty the tmp_log_file
	  echo "" > ${LOG_FILE_NAME}

	  #echo "Test log input: " >> ${LOG_FILE_NAME}
	  time ./lamure_points_to_lines_obj_exporter -f $bvh_file -t $rot_file -d $depth --apply_nurbs_fitting --apply_alpha_shapes > ${LOG_FILE_NAME}

	  cat ${LOG_FILE_NAME} > ${bvh_filename_wo_extensions}"."${rot_filename_wo_extensions}".LOG"

	  time ./lamure_points_to_lines_obj_exporter -f $bvh_file -t $rot_file -d $depth --write_xyz_points

	  for xyz_all_file in *.xyz_all
	  do
		./lamure_preprocessing ${xyz_all_file}
	  done
	done


	  mv *.obj ${input_folder}processed/${bvh_filename_wo_extensions}/line_objs/.
	  mv ./*.bvh ${input_folder}processed/${bvh_filename_wo_extensions}/point_clouds/.
	  mv ./*.lod ${input_folder}processed/${bvh_filename_wo_extensions}/point_clouds/.
	  mv ./*.xyz_all ${input_folder}processed/${bvh_filename_wo_extensions}/point_clouds/.
	  mv ./*.LOG ${input_folder}processed/${bvh_filename_wo_extensions}/LOGS/.
done

#
#LOG_FILE_NAME_2=xyz_log
#touch ${LOG_FILE_NAME_2}
#for xyz_file in ${input_rot_mat_folder}*.xyz_all
#do
#   ./lamure_preprocessing ${xyz_all_file}.xyz_all >> ${LOG_FILE_NAME_2}
#done

	


