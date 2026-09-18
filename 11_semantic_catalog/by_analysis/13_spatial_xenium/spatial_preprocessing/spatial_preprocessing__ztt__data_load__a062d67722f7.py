import os
import threading
import scanpy as sc
import spatialdata as sd
from spatialdata_io import xenium


# 多线程读取Xenium下机数据读取
def xenium_data_load_multithreaded(data_dir, sample_info):
    def sd_read_xenium(sample_data, sample_name, sdata_dict):
        sdata = xenium(path=sample_data, cells_boundaries=True, n_jobs=6)
        sdata_dict[sample_name] = sdata

    threads = []
    sdata_dict = {}
    sample_2_group = {}
    with open(sample_info, 'r') as f:
        for line in f:
            raw_name, sample_name, group_name = line.strip().split('\t')[:3]
            sample_2_group[sample_name] = group_name
            thread = threading.Thread(target=sd_read_xenium, args=(os.path.join(data_dir, raw_name),sample_name, sdata_dict,))
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()
    
    sdata = sd.concatenate(
                sdata_dict,
                concatenate_tables=True,
                obs_names_make_unique=True
            )

    sdata.tables['table'].obs["sample"] = sdata.tables['table'].obs["region"].str.replace('cell_circles-', '')
    sdata.tables['table'].obs["group"] = sdata.tables['table'].obs["sample"].apply(lambda x: sample_2_group[x])

    sdata.tables['table'].obs["cell_boundaries"] = sdata.tables['table'].obs["region"].str.replace('cell_circles', 'cell_boundaries')
    sdata.set_table_annotates_spatialelement(table_name='table', region=[i for i in sdata.shapes.keys() if i.startswith('cell_boundaries-')], region_key='cell_boundaries')

    return sdata


if __name__ == '__main__':
    data_dir = '/mnt/NexData_file/home/ztt/test_data'
    sample_info = './sample.txt'
    sdata, adata = xenium_data_load_multithreaded(data_dir, sample_info)






