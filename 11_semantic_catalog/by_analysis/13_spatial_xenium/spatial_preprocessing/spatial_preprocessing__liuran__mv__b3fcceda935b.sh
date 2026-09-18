#!/bin/bash
# mv.sh

# 使用循环处理每个匹配的目录
for dir in ./output-XETG00305__0064340__TMA_hum__20260121__052026; do
    if [ -d "$dir" ]; then
        echo "处理目录: $dir"
        
        # 检查源文件是否存在
        if [ -f "$dir/morphology_focus/ch0000_dapi.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0000_dapi.ome.tif" "$dir/morphology_focus/morphology_focus_0000.ome.tif"
        fi
        
        if [ -f "$dir/morphology_focus/ch0001_atp1a1_cd45_e-cadherin.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0001_atp1a1_cd45_e-cadherin.ome.tif" "$dir/morphology_focus/morphology_focus_0001.ome.tif"
        fi
        
        if [ -f "$dir/morphology_focus/ch0002_18s.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0002_18s.ome.tif" "$dir/morphology_focus/morphology_focus_0002.ome.tif"
        fi
        
        if [ -f "$dir/morphology_focus/ch0003_alphasma_vimentin.ome.tif" ]; then
            mv "$dir/morphology_focus/ch0003_alphasma_vimentin.ome.tif" "$dir/morphology_focus/morphology_focus_0003.ome.tif"
        fi
    fi
done