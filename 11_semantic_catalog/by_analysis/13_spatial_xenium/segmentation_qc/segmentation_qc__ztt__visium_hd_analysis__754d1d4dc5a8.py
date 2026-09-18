import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import json
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
import seaborn as sns
import geopandas as gpd
from PIL import Image
from spatialdata_io import visium_hd, geojson
from spatialdata.models import TableModel,ShapesModel
from spatialdata_io._constants._constants import VisiumHDKeys
from spatialdata.transformations.transformations import Identity
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context

import warnings
warnings.filterwarnings("ignore")

plt.style.use('default')
plt.rcParams['figure.facecolor'] = 'white'


def read_visiumhd_segmented(data_dir, sample_info, bin_size=8, boundary='cell_segmentations'):
    sample2group = {}
    with open(sample_info, 'r') as f:
        for line in f.readlines():
            line = line.strip().split('\t')
            if len(line) > 1:
                sample2group[line[0]] = line[1]
            else:
                print('sample.txt need two columns')
    
    adata_list = []
    for sample in os.listdir(data_dir):
        sdata = visium_hd(
            path=os.path.join(data_dir, sample),
            dataset_id=sample,
            filtered_counts_file=True,
            var_names_make_unique=True,
            bin_size=bin_size
        )
        with open(os.path.join(data_dir, sample, 'segmented_outputs/spatial/scalefactors_json.json')) as file:
            scalefactors = json.load(file)
        # Estimate spot diameter
        hires_img = np.asarray(Image.open(os.path.join(data_dir, sample, 'segmented_outputs/spatial/tissue_hires_image.png')))
        lowres_img = np.asarray(Image.open(os.path.join(data_dir, sample, 'segmented_outputs/spatial/tissue_lowres_image.png')))
        shapes = geojson(os.path.join(data_dir, sample, "segmented_outputs/cell_segmentations.geojson"), coordinate_system=sample)
        areas = shapes.geometry.area
        mean_area = np.mean(areas)
        diameter_pixels = 2 * np.sqrt(mean_area / np.pi)
        scalefactors["spot_diameter_fullres"] = diameter_pixels
        uns_spatial = {
                sample: {
                    "images": {
                        "hires": hires_img,
                        "lowres": lowres_img,
                    },
                    "scalefactors": scalefactors
                }
            }

        if boundary=='cell_segmentations' and os.path.exists(os.path.join(data_dir, sample, "segmented_outputs")):
            # add segmentated scRNA data
            adata = sc.read_10x_h5(os.path.join(data_dir, sample, "segmented_outputs/filtered_feature_cell_matrix.h5"), gex_only=False)
            adata.var_names_make_unique()
            adata.obs['sample'] = sample
            adata.obs['group'] = sample2group[sample]
            shapes_name = "cell_segmentation"
            adata.obs[VisiumHDKeys.INSTANCE_KEY] = adata.obs.index.str.replace('cellid_', '').str.replace('-1','').astype(int)
            adata.obs[VisiumHDKeys.REGION_KEY] = shapes_name
            adata.obs[VisiumHDKeys.REGION_KEY] = adata.obs[VisiumHDKeys.REGION_KEY].astype("category")
            adata = TableModel.parse(
                adata=adata,
                region=shapes_name,
                region_key=str(VisiumHDKeys.REGION_KEY),
                instance_key=str(VisiumHDKeys.INSTANCE_KEY)
            )
            sdata.tables['cells'] = adata.copy()

            # add cell_segmentation boundary
            gdf_seg = gpd.read_file(os.path.join(data_dir, sample, f'segmented_outputs/{boundary}.geojson'))
            gdf_seg.rename(columns={'cell_id':'location_id'}, inplace=True)
            gdf_seg = gdf_seg[gdf_seg['location_id'].isin(adata.obs['location_id'])]
            gdf_seg.set_index('location_id', inplace=True)
            sdata.shapes['cell_segmentation'] = ShapesModel.parse(gdf_seg, transformations={sample: Identity()})

            # add spatial_coords into adata.obsm['spatial]
            gdf = sdata.shapes['cell_segmentation']
            gdf.crs = None
            spatial_coords = np.column_stack([
                gdf.centroid.x.values, 
                gdf.centroid.y.values
            ])
            sdata.tables['cells'].obsm['spatial'] = spatial_coords.copy()

            # add cell area in adata.obs, cell area unit is um2
            gdf_seg['cell_area'] = gdf_seg.area * scalefactors['microns_per_pixel'] * scalefactors['microns_per_pixel']
            sdata.tables['cells'].obs = pd.merge(sdata.tables['cells'].obs, gdf_seg['cell_area'], on='location_id', how='inner')

            sdata.tables['cells'].uns["spatial"] = uns_spatial
            adata_list.append(sdata.tables['cells'])
        else:
            tb_name = list(sdata.tables.keys())[0]
            tmp_adata = sdata.tables[tb_name]
            tmp_adata.obs['sample'] = sample
            tmp_adata.obs['group'] = sample2group[sample]
            tmp_adata.var_names_make_unique()

            tmp_adata.uns["spatial"] = uns_spatial
            adata_list.append(tmp_adata)

    adata = ad.concat(adata_list, join='inner', merge='same', uns_merge='first', index_unique=None)
    adata.obs.reset_index(inplace=True, drop=True)
    adata.obs.index = adata.obs.index.astype(str)

    # Data QC
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-') 
    adata.var['ribo'] = adata.var_names.str.upper().str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes_by_counts']) / np.log10(adata.obs['total_counts'])

    return adata











