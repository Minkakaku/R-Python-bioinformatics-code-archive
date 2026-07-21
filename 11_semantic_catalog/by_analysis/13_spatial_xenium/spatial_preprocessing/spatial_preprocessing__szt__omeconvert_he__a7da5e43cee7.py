import tifffile as tf
import sys
import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.ndimage import affine_transform

def image_affine_transform(img, affine_mat, type = "HE"):
    res_img = None
    if len(img.shape) == 3:
        if not (type == "HE" or type == "qptiff"):
            raise Exception("type should only be HE or qptiff")
        
        img_ls = []
        if type == "qptiff":
            for i in range(img.shape[0]):
                shape = img[i, :, :].shape
                img_one = affine_transform(img[i, :, :], affine_mat, output_shape = (max(shape), max(shape)))
                img_ls.append(img_one)
            res_img = np.stack(img_ls, axis=0)
            
        elif type == "HE":
            for i in range(img.shape[2]):
                shape = img[:, :, i].shape
                img_one = affine_transform(img[:, :, i], affine_mat, output_shape = (max(shape), max(shape)))
                img_ls.append(img_one)
            res_img = np.stack(img_ls, axis=2)
    else:
        res_img = affine_transform(img, affine_mat, output_shape = (max(img.shape, img.shape)))
    
    return res_img


def img_resize(img,scale_factor):
    width = int(np.floor(img.shape[1] * scale_factor))
    height = int(np.floor(img.shape[0] * scale_factor))
    return cv2.resize(img, (width, height), interpolation = cv2.INTER_AREA)

def img_resize_full(img, scale_factor, type = "HE"):
    res_img = None
    if not (type == "HE" or type == "qptiff"):
            raise Exception("type should only be HE or qptiff")
    
    if len(img.shape) == 3 and type == "qptiff":
        img_ls = []
        for i in range(img.shape[0]):
            img_one = img_resize(img[i, :, :], scale_factor)
            img_ls.append(img_one)
            
        res_img = np.stack(img_ls, axis=0)

    else:
        res_img = img_resize(img, scale_factor)
    
    return res_img

def write_ome_tif(filename, image, channel_names, photometric_interp, metadata, subresolutions):
    subresolutions = subresolutions

    fn = filename + ".ome.tif"

    with tf.TiffWriter(fn,  bigtiff=True, ) as tif:
        px_size_x = metadata['PhysicalSizeX']
        px_size_y = metadata['PhysicalSizeY']

        options = dict(
            photometric=photometric_interp,
            tile=(1024, 1024),
            maxworkers=4,
            compression='jpeg2000',
            compressionargs={'level':85},
            resolutionunit='CENTIMETER',
        )

        print("Writing pyramid level 0")
        tif.write(
            image,
            subifds=subresolutions,
            resolution=(1e4 / px_size_x, 1e4 / px_size_y),
            metadata=metadata,
            **options
        )

        scale = 1
        for i in range(subresolutions):
            scale /= 2
            if photometric_interp == 'minisblack':
                if image.shape[0] < image.shape[-1]:
                    image = np.moveaxis(image,0,-1)
                    image = img_resize_full(image,0.5,type = "HE")
                    image = np.moveaxis(image,-1,0)
            else:
                image = img_resize_full(image,0.5,type = "HE")

            print("Writing pyramid level {}".format(i+1))
            tif.write(
                image,
                subfiletype=1,
                resolution=(1e4 / scale / px_size_x, 1e4 / scale / px_size_y),
                **options
            )


with tf.TiffFile(sys.argv[1]) as tif:
    image = tif.asarray()

    channel_names=None

    if tif.pages[0].samplesperpixel == 3:
        photometric_interp='rgb'
    elif tif.pages[0].samplesperpixel == 1:
        photometric_interp='minisblack'

    if tif.pages[0].get_resolution():
        res = tif.pages[0].get_resolution()
        res = [1/r for r in res]
        if tif.pages[0].resolutionunit == 2:
            unit = "in"
        elif tif.pages[0].resolutionunit == 3:
            unit = "cm"
        elif tif.pages[0].resolutionunit == 4:
            unit = "mm"
        elif tif.pages[0].resolutionunit == 5:
            unit = "um"

    if tif.is_ome:
        meta = tf.xml2dict(tif.ome_metadata)
        res = (meta['OME']['Image']['Pixels']['PhysicalSizeX'],meta['OME']['Image']['Pixels']['PhysicalSizeY'])
        unit = meta['OME']['Image']['Pixels']['PhysicalSizeXUnit']

        try:
            channel_names=[]
            for i, element in enumerate(meta['OME']['Image']['Pixels']['Channel']):
                channel_names.append(meta['OME']['Image']['Pixels']['Channel'][i]['Name'])
        except KeyError:
            channel_names=None

    metadata={
        'PhysicalSizeX': res[0],
        'PhysicalSizeXUnit': unit,
        'PhysicalSizeY': res[1],
        'PhysicalSizeYUnit': unit,
        'Channel': {'Name': channel_names}
        }

filename = sys.argv[1].rsplit('.',1)[0]

# crop image
# (w1, h1, w2, h2)
if len(sys.argv) > 2: 
    # parse coords
    coords = eval(sys.argv[2]) 
    print(f"coords: {coords}")

    if len(coords) ==  1 and coords == -1:
        w1, h1, w2, h2 = 0, 0, image.shape[0], image.shape[1]
    elif len(coords) != 4:
        raise Exception("coordinates should have length of 4.")
    w1, h1, w2, h2 = coords

    # crop
    if len(image.shape) == 3: 
        image = image[h1:h2, w1:w2, :] 
    elif len(image.shape) == 2: 
        image = image[h1:h2, w1:w2]  

# sample id
try:
    sample_id = sys.argv[3]
except Exception:
    sample_id = None


# if have affine matrix
if len(sys.argv) > 4:
    affine_matrix = eval(sys.argv[4]) 
    affine_matrix = np.array(affine_matrix).reshape((2, -1))
    print(f"affine_matrix: {affine_matrix}")

    image = image_affine_transform(image, affine_matrix, type = "HE")


# export thumbnail
plt.rcParams['figure.figsize'] = (8,8)
plt.imshow(img_resize_full(image, scale_factor=0.1, type = "HE"))
thumbnail_filename = f"{filename} - {sample_id}.png" if sample_id is not None else f"{filename}.png"
plt.savefig(thumbnail_filename)

# export image
filename = f"{filename} - {sample_id}" if sample_id is not None else filename
write_ome_tif(filename, image, channel_names, photometric_interp, metadata, subresolutions = 5)

# conda create --name tiff tifffile opencv
# conda activate tiff
# Usage: python omecovert.py <path to image file being converted>
# python omeconvert.py Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif

# parameters list
# @1：path of current script file
# @2：path of image to be processed
# @3：coords of ROI; -1 means not crop image
# @4：sample id to rename the new image
# @5：affine matrix

# affine matrix
# sample1: 0.002465042,1.138533541,1002.355568,-1.138533541,0.002465042,33820.7169
# sample2: 0.002716359,1.138724738,659.508548,-1.138724738,0.002716359,29063.62312
# sample3: 0.002743924,1.138161928,-848.0258569,-1.138161928,0.002743924,32916.96413


# cd data
# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.tif 796,7860,29407,48106 sample1 "0.002465042,1.138533541,1002.355568,-1.138533541,0.002465042,33820.7169"
# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.tif 31798,5230,58216,46960 sample2 "0.002716359,1.138724738,659.508548,-1.138724738,0.002716359,29063.62312"
# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.tif 58232,1455,86604,49114 sample3 "0.002743924,1.138161928,-848.0258569,-1.138161928,0.002743924,32916.96413"