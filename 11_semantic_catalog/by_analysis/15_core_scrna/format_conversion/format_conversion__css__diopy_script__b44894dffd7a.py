{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import anndata as ad\n",
    "import diopy\n",
    "import argparse\n",
    "import warnings\n",
    "import os"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "usage: diopy.py [-h] -i H5AD -o OUT_DIR -n NAME\n",
      "diopy.py: error: the following arguments are required: -i, -o, -n\n"
     ]
    },
    {
     "ename": "SystemExit",
     "evalue": "2",
     "output_type": "error",
     "traceback": [
      "An exception has occurred, use %tb to see the full traceback.\n",
      "\u001b[0;31mSystemExit\u001b[0m\u001b[0;31m:\u001b[0m 2\n"
     ]
    }
   ],
   "source": [
    "warnings.filterwarnings(\"ignore\")\n",
    "parser = argparse.ArgumentParser(prog='diopy.py', description=\"h5ad2h5\")\n",
    "parser.add_argument('-i', action='store', dest='h5ad', required=True, help='adata.h5ad')\n",
    "parser.add_argument('-o', action='store', dest='out_dir',required=True, help='Result output directory.')\n",
    "parser.add_argument('-n', action='store', dest='name',required=True, help='Result output h5name.')\n",
    "parser = parser.parse_args()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "adata = ad.read_h5ad(parser.h5ad)\n",
    "os.chdir(parser.out_dir)\n",
    "adata.obs_names_make_unique()\n",
    "# 这里有时候会报错\n",
    "# diopy.output.write_h5(adata,file = parser.name,save_X=True,save_graph = True)\n",
    "diopy.output.write_h5(adata,file = parser.name,save_X=False,save_graph = True)\n",
    "os.system(f'sh /home/wzb/pipline/dio/docker_run.sh {os.path.join(parser.out_dir,parser.name)} {parser.out_dir} \"h5_trans.rds\"')"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "scRNA",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.16"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
