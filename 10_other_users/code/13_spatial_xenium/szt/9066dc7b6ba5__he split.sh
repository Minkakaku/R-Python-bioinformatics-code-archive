# Qupath acquire coordinates
# GsonTools.getInstance(true).toJson(getAnnotationObjects())

# INFO: Result: [
#   {
#     "type": "Feature",
#     "id": "4ed7db6d-8182-449c-b7bb-b61ade8d26cf",
#     "geometry": {
#       "type": "Polygon",
#       "coordinates": [
#         [
#           [58232, 1455],
#           [86604, 1455],
#           [86604, 49114],
#           [58232, 49114],
#           [58232, 1455]
#         ]
#       ]
#     },
#     "properties": {
#       "objectType": "annotation",
#       "name": "ID100"
#     }
#   },
#   {
#     "type": "Feature",
#     "id": "10965d3f-6d77-4d45-b605-54f063e9cde9",
#     "geometry": {
#       "type": "Polygon",
#       "coordinates": [
#         [
#           [796, 7860],
#           [29407, 7860],
#           [29407, 48106],
#           [796, 48106],
#           [796, 7860]
#         ]
#       ]
#     },
#     "properties": {
#       "objectType": "annotation",
#       "name": "ID1"
#     }
#   },
#   {
#     "type": "Feature",
#     "id": "adceb5f3-8bb6-44a1-a31f-ee9fb1914557",
#     "geometry": {
#       "type": "Polygon",
#       "coordinates": [
#         [
#           [31798, 5230],
#           [58216, 5230],
#           [58216, 46960],
#           [31798, 46960],
#           [31798, 5230]
#         ]
#       ]
#     },
#     "properties": {
#       "objectType": "annotation",
#       "name": "ID10"
#     }
#   }
# ]

# ID1
#           [796, 7860],
#           [29407, 7860],
#           [29407, 48106],
#           [796, 48106],
#           [796, 7860]
# point1 and point3 -> 796,7860,29407,48106

# ID10
#           [31798, 5230],
#           [58216, 5230],
#           [58216, 46960],
#           [31798, 46960],
#           [31798, 5230]
# point1 and point3 -> 31798,5230,58216,46960

# ID100
#           [58232, 1455],
#           [86604, 1455],
#           [86604, 49114],
#           [58232, 49114],
#           [58232, 1455]
# point1 and point3 -> 58232,1455,86604,49114

# run python code: 
# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.tif  

# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.ome.tif 796,7860,29407,48106 sample1
# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.ome.tif 31798,5230,58216,46960 sample2
# python ../omeconvert_he.py 原老师xenium+PCF下机-HE.ome.tif 58232,1455,86604,49114 sample3