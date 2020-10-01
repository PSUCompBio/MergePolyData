import sys
import json
import linecache

jsonOutputFile = sys.argv[1]

outputJson = json.loads(open(jsonOutputFile).read())
meshType = outputJson["output-file"].split('_')[0]
cellDataFile = meshType'_cellcentres.txt'

# Part list 
partList = ["msc", "stem", "cerebellum", "frontal", "parietal", "occipital",
        "temporal"]

# Populate region from cell centres file
maxInjuryMetrics = ["principal-max-strain", "principal-min-strain",
        "maximum-shear-strain", "maximum-PSxSR"];
threshInjuryMetrics = ["CSDM-10", "CSDM-15", "CSDM-30", "CSDM-5", 
        "MPS-95"];

# Loop over metric with single element output
for metric in maxInjuryMetrics:
    maxLocation = int(outputJson[metric]["global-element-id"])
    maxLine = linecache.getline('cellcentres.txt', maxLocation).split()
    outputJson[metric]["brain-region"] = maxLine[4]
    outputJson[metric]["location-femtech-reference"] = \
            [float(maxLine[1]), float(maxLine[2]), float(maxLine[3])]
    outputJson[metric].pop("global-element-id", None)
    outputJson[metric].pop("location", None)

# Loop over metric with multiple element output
for metric in threshInjuryMetrics:
    maxLocation = outputJson[metric]["global-element-id"]
    for part in partList:
        outputJson[metric][part] = []
    for loc in maxLocation:
        maxLine = linecache.getline('cellcentres.txt', int(loc)).split()
        outputJson[metric][maxLine[4]].append([float(maxLine[1]), \
                float(maxLine[2]), float(maxLine[3])])
    outputJson[metric].pop("global-element-id", None)

# jstr = json.dumps(outputJson, indent = 2, sort_keys=False)
# print(jstr)
with open(jsonOutputFile, 'w') as outfile:
    json.dump(outputJson, outfile, indent = 2, sort_keys=False)
