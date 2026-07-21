cd Xenium_Align

#28
python /mnt/NexData_file/home/hhq/Xenium-Align/Xenium_Align/data_preprocess_check.py -sample 28 -preservation_method ffpe -data_file_path /mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/FWSP20255020-02/output-XETG00305__0074658__28__20251017__072721

python /mnt/NexData_file/home/hhq/Xenium-Align/Xenium_Align/cellpose_image_segmentation.py -sample 28 -preservation_method ffpe -data_file_path /mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/FWSP20255020-02/output-XETG00305__0074658__28__20251017__072721 -channel_cellpose 1 -min_size 15 -flow_threshold 0.8 -use_gpu True


python /mnt/NexData_file/home/hhq/Xenium-Align/Xenium_Align/xenium_alignment_for_keypoints.py -sample 28 -preservation_method ffpe -data_file_path /mnt/NexData_file/home/hhq/project/Xenium/FWSP20255020-02_mouse_lung_xenium/FWSP20255020-02/output-XETG00305__0074658__28__20251017__072721 -crop_radius_pixel 400 -center_move_pixel 300 -check_cell_num 100 -mip_ome_extract_ratio 0.125 -mip_ome_extract_min 50 -segment_method cellpose -overlap_type overlap_ave -overlap_threshold_ave 0.9 -keypoints_min_num 15 -epoch_num 30