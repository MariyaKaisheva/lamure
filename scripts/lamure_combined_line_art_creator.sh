#!/bin/bash

lamure_app_parameter_array=( "$@" )

echo "Number of Parameters $#"

base_command_to_execute=""

rotation_array_size=0
for arg in "${lamure_app_parameter_array[@]}"; do
  echo "$arg"

  #check if argument ends on '.rot', where the point needs to be escaped
  #if the point would not be escaped, it would be a regex placeholder 
  if [[ $arg == *\.rot ]]; then
    rotation_file_array[$rotation_array_size]=$arg
    rotation_array_size=$(($rotation_array_size+1))
   
    echo "FOUND ROT FILE"
  else
    if [ $arg == "-t" ]; then
      echo "IGNORING '-t' FLAG"
    else
      base_command_to_execute="$base_command_to_execute $arg"
    fi
  fi
done

echo "rot array size: $rotation_array_size"

echo "Base command to execute: $base_command_to_execute"


mkdir combined_files
cd combined_files
echo "" > combined_line_art_file.lob

processing_iteration_counter=0
for rotation_string in "${rotation_file_array[@]}"; do
 
  final_command_to_execute="$base_command_to_execute -t $rotation_string"
  echo "EXECUTING COMMAND $final_command_to_execute"


  command_output=`eval "$final_command_to_execute"`


  output_lob_file=`echo "$command_output" | grep "lob"`
  echo "OUTPUT LOB FILENAME IN ITERATION: $output_lob_file"

  cat $output_lob_file >> combined_line_art_file.lob
 
  echo ""
  echo ""
  echo ""  
done

cd ..
