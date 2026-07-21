import tifffile as tf
import sys
import numpy as np
import cv2
import matplotlib.pyplot as plt

def img_resize(img,scale_factor):
    width = int(np.floor(img.shape[1] * scale_factor))
    height = int(np.floor(img.shape[0] * scale_factor))
    return cv2.resize(img, (width, height), interpolation = cv2.INTER_AREA)

def img_resize_full(img, scale_factor):
    res_img = None
    if len(img.shape) == 3:
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
            image = img_resize_full(image,0.5)

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
        # res = tif.pages[0].get_resolution()
        x_res = tif.pages[0].tags['XResolution'].value
        y_res = tif.pages[0].tags['XResolution'].value

        # x_res: (n, d) -> n：pixel, d: phsical length
        if len(x_res) == 2: x_res = x_res[0]/x_res[1]
        if len(y_res) == 2: y_res = y_res[0]/y_res[1]

        res = (1/x_res, 1/y_res)

        if tif.pages[0].resolutionunit == 2:
            unit = "in"
        elif tif.pages[0].resolutionunit == 3:
            unit = "cm"
        elif tif.pages[0].resolutionunit == 4:
            unit = "mm"
        elif tif.pages[0].resolutionunit == 5:
            unit = "um"
        # else:
        #     res = (1, 1)
        #     unit = None    

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
    else:
        channel_names=[]
        for i in range(image.shape[0]):
            meta = tf.xml2dict(tif.pages[i].description)
            key = list(meta.keys())[0]
            channel_names.append(meta[key]['Biomarker'])


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
    if len(coords) != 4:
        raise Exception("coordinates should have length of 4.")
    w1, h1, w2, h2 = coords

    # crop
    if len(image.shape) == 3: 
        image = image[:, h1:h2, w1:w2] 
    elif len(image.shape) == 2: 
        image = image[h1:h2, w1:w2]  

# sample id
try:
    sample_id = sys.argv[3]
except Exception:
    sample_id = None

# export thumbnail
plt.rcParams['figure.figsize'] = (8,8)
plt.imshow(img_resize_full(image[0, :, :], scale_factor=0.1))
thumbnail_filename = f"{filename} - {sample_id}.png" if sample_id is not None else f"{filename}.png"
plt.savefig(thumbnail_filename)

# export image
filename = f"{filename} - {sample_id}" if sample_id is not None else filename
write_ome_tif(filename, image, channel_names, photometric_interp, metadata, subresolutions = 5)

# conda create --name tiff tifffile opencv
# conda activate tiff
# Usage: python omecovert.py <path to image file being converted>
# python omeconvert.py Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif


# Qupath acquire coordinates
# GsonTools.getInstance(true).toJson(getAnnotationObjects())

#   "coordinates": [
#     [
#       [4120, 3764],
#       [23111, 3764],
#       [23111, 17303],
#       [4120, 17303],
#       [4120, 3764]
#     ]
#   ]
# point1 and point3 -> 4120,3764,23111,17303

#   [3803, 18012],
#   [24014, 18012],
#   [24014, 30932],
#   [3803, 30932],
#   [3803, 18012]
# point1 and point3 -> 3803,18012,24014,30932

#   [3446, 30993],
#   [25783, 30993],
#   [25783, 44986],
#   [3446, 44986],
#   [3446, 30993]
# point1 and point3 -> 3446,30993,25783,44986

# run python code: 
# python ../omeconvert_qptiff.py FWSP20241617\ yuan_Scan1.er.qptiff 4120,3764,23111,17303 sample1
# python ../omeconvert_qptiff.py FWSP20241617\ yuan_Scan1.er.qptiff 3803,18012,24014,30932 sample2
# python ../omeconvert_qptiff.py FWSP20241617\ yuan_Scan1.er.qptiff 3446,30993,25783,44986 sample3