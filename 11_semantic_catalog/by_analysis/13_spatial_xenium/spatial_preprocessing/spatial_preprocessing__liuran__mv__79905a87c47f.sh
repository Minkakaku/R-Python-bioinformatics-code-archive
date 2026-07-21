#!/bin/bash
# mv.sh

# 查找所有匹配的目录
for dir in ./output-XETG00305__0076917__*__20260313__054531; do
    if [ -d "$dir" ]; then
        echo "Processing: $dir"
        
        # 检查源文件是否存在
        if [ -f "$dir/morphology_focus/ch0000_dapi.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0000_dapi.ome.tif" \
               "$dir/morphology_focus/morphology_focus_0000.ome.tif"
        fi
        
        if [ -f "$dir/morphology_focus/ch0001_atp1a1_cd45_e-cadherin.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0001_atp1a1_cd45_e-cadherin.ome.tif" \
               "$dir/morphology_focus/morphology_focus_0001.ome.tif"
        fi
        
        if [ -f "$dir/morphology_focus/ch0002_18s.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0002_18s.ome.tif" \
               "$dir/morphology_focus/morphology_focus_0002.ome.tif"
        fi
        
        if [ -f "$dir/morphology_focus/ch0003_alphasma_vimentin.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0003_alphasma_vimentin.ome.tif" \
               "$dir/morphology_focus/morphology_focus_0003.ome.tif"
        fi
    fi
done