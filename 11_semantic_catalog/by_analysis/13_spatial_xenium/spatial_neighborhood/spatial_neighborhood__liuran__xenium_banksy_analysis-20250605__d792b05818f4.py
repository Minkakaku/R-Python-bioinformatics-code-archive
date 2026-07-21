import os
import sys
import time
import logging
import warnings
import psutil  # 用于监控内存使用
import threading

# Suppress warnings
warnings.filterwarnings("ignore")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("xenium_banksy_script.log"),
        logging.StreamHandler(sys.stdout)
    ]
)

# Add custom directory to system path
sys.path.append('/home/liuran/software/Banksy_py')

# Import required libraries
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import umap
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.cluster_methods import run_Leiden_partition

# Global variable to store memory usage data
memory_usage_data = []
remaining_memory_data = []

def monitor_memory_usage(interval=10):
    """
    Monitor memory usage at specified intervals and log the data.
    
    Args:
        interval (int): Interval in seconds between memory checks.
    """
    process = psutil.Process(os.getpid())
    while True:
        # 获取当前进程的内存使用情况（驻留集大小，RSS）
        mem_info = process.memory_info()
        memory_usage = mem_info.rss / (1024 ** 2)  # Convert to MB
        
        # 获取服务器总内存和剩余内存
        mem = psutil.virtual_memory()
        total_memory = mem.total / (1024 ** 2)  # Convert to MB
        remaining_memory = mem.available / (1024 ** 2)  # Convert to MB
        
        # 记录内存使用数据
        memory_usage_data.append(memory_usage)
        remaining_memory_data.append(remaining_memory)
        
        # 记录日志
        logging.info(f"Current memory usage: {memory_usage:.2f} MB | Remaining memory: {remaining_memory:.2f} MB")
        
        time.sleep(interval)

def main():
    try:
        start_time = time.time()
        logging.info("Starting Xenium Banksy Analysis Script")
        
        # Start memory monitoring in a separate thread
        memory_thread = threading.Thread(target=monitor_memory_usage, args=(10,))
        memory_thread.daemon = True  # Ensure the thread exits when the main program exits
        memory_thread.start()
        
        # Read the input data
        logging.info("Reading input data from xenium_for_banksy.h5ad...")
        adata = sc.read_h5ad('xenium_for_banksy.h5ad')
        
        # Banksy parameters
        logging.info("Setting Banksy parameters...")
        nbr_weight_decay = 'scaled_gaussian'
        k_geom = 15
        lambda_list = [0.2]  # list of lambda parameters
        m = 1
        
        # Include spatial coordinates information
        logging.info("Including spatial coordinates information...")
        coord_keys = ('x_pixel', 'y_pixel', 'coord_xy')
        x_coord, y_coord, xy_coord = coord_keys[0], coord_keys[1], coord_keys[2]
        raw_y, raw_x = adata.obs[y_coord], adata.obs[x_coord]
        adata.obsm[xy_coord] = np.vstack((adata.obs[x_coord].values, adata.obs[y_coord].values)).T
        
        # Run Banksy
        logging.info("Initializing Banksy...")
        banksy_dict = initialize_banksy(
            adata,
            coord_keys,
            k_geom,
            nbr_weight_decay=nbr_weight_decay,
            max_m=m,
            plt_edge_hist=False,
            plt_nbr_weights=False,
            plt_agf_angles=False,
            plt_theta=False
        )
        
        logging.info("Generating Banksy matrix...")
        banksy_dict, banksy_matrix = generate_banksy_matrix(
            adata,
            banksy_dict,
            lambda_list,
            max_m=m
        )
        
        # PCA and Harmony integration
        logging.info("Performing PCA and Harmony integration...")
        sc.tl.pca(banksy_matrix)
        sce.pp.harmony_integrate(banksy_matrix, "batch")
        
        # Neighbors and UMAP
        logging.info("Computing neighbors and UMAP...")
        sc.pp.neighbors(banksy_matrix, use_rep='X_pca_harmony')
        sc.tl.umap(banksy_matrix)
        
        # Leiden clustering
        logging.info("Performing Leiden clustering...")
        sc.tl.leiden(banksy_matrix)
        
        # Save results
        logging.info("Saving processed data to xenium_for_banksy_0.2.h5ad...")
        banksy_matrix.write('xenium_for_banksy_0.2.h5ad')
        
        logging.info("Saving observation data to xenium_for_banksy_0.2_obs.csv...")
        banksy_matrix.obs.to_csv("xenium_for_banksy_0.2_obs.csv", index=True)
        
        end_time = time.time()
        logging.info("Script completed successfully!")
        logging.info(f"Total execution time: {end_time - start_time:.2f} seconds")
        
        # Save memory usage data to a file
        with open("memory_usage.log", "w") as f:
            f.write("Memory Usage (MB),Remaining Memory (MB)\n")
            for mem, rem in zip(memory_usage_data, remaining_memory_data):
                f.write(f"{mem:.2f},{rem:.2f}\n")
        logging.info("Memory usage data saved to memory_usage.log")
        
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise e

if __name__ == "__main__":
    import time  # Import the time module here
    main()