#!/bin/bash
# run_cellchat_contrast.sh
source activate R_envs
while IFS= read -r line; do
  trace=$(echo "$line" | awk -F'[[:space:]]+' '{gsub(/"/, "", $1); print $1}')
  case=$(echo "$line" | awk -F'[[:space:]]+' '{gsub(/"/, "", $2); print $2}')
  echo "  trace: $trace"
  echo "  case: $case"
  echo "Running: $trace vs $case"
  Rscript /mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/cellchat/contrast/cellchat_contrast.R \
    --pos.dataset "$trace" \
    --other.dataset "$case" \
    --pos_rds /mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/cellchat/"$trace"/"$trace"_cellchat_netP.rds \
    --other_rds /mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/cellchat/"$case"/"$case"_cellchat_netP.rds \
    --output /mnt/NexData_file/home/wyt/S203_project/project/xenium_FWSP20255467_human_skin/anno_A_NC/cellchat/contrast/
    
  done < skin_contrace.txt

