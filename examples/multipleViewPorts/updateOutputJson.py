import sys
import json
import linecache

jsonOutputFile = sys.argv[1]

outputJson = json.loads(open(jsonOutputFile).read())

# Populate region from cell centres file
injuryMetrics = ["principal-max-strain", "principal-min-strain",
        "maximum-shear-strain", "maximum-PSxSR"];
for metric in injuryMetrics:
    maxLocation = int(outputJson[metric]["global-element-id"])
    maxLine = linecache.getline('cellcentres.txt', maxLocation).split()
    outputJson[metric]["brain-region"] = maxLine[4]
    outputJson[metric]["location-femtech-reference"] = \
            [float(maxLine[1]), float(maxLine[2]), float(maxLine[3])]

# jstr = json.dumps(outputJson, indent = 2, sort_keys=False)
# print(jstr)
with open(jsonOutputFile, 'w') as outfile:
    json.dump(outputJson, outfile, indent = 2, sort_keys=False)
