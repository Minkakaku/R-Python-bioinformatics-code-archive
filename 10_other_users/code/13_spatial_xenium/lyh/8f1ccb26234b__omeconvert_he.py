import tifffile as tf
import sys
import numpy as np
import cv2
import matplotlib.pyplot as plt

def img_resize(img,scale_factor):
    width = int(np.floor(img.shape[1] * scale_factor))
    height = int(np.floor(img.shape[0] * scale_factor))
    return cv2.resize(img, (width, height), interpolation = cv2.INTER_AREA)

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
                    image = img_resize(image,0.5)
                    image = np.moveaxis(image,-1,0)
            else:
                image = img_resize(image,0.5)

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
    if len(coords) != 4:
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

# export thumbnail
plt.rcParams['figure.figsize'] = (8,8)
plt.imshow(img_resize(image, scale_factor=0.1))
thumbnail_filename = f"{filename} - {sample_id}.png" if sample_id is not None else f"{filename}.png"
plt.savefig(thumbnail_filename)

# export image
filename = f"{filename} - {sample_id}" if sample_id is not None else filename
write_ome_tif(filename, image, channel_names, photometric_interp, metadata, subresolutions = 5)

# conda create --name tiff tifffile opencv matplotlib
# conda activate tiff
# Usage: python omecovert.py <path to image file being converted>
# python omeconvert.py Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif
