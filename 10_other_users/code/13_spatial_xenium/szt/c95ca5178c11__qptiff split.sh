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
# python ../omeconvert.py FWSP20241617\ yuan_Scan1.er.qptiff 4120,3764,23111,17303 sample1
# python ../omeconvert.py FWSP20241617\ yuan_Scan1.er.qptiff 3803,18012,24014,30932 sample2
# python ../omeconvert.py FWSP20241617\ yuan_Scan1.er.qptiff 3446,30993,25783,44986 sample3