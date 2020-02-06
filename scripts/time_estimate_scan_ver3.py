import json

def TimeEstimate_for_Xray():
    return

f = open('source_data/xray_data.json', 'w')

f.close()

Am-241 = 'E=60.0keV'
XrayTube = ['E=' + str(i) +'0.0keV' for i in range(1,6)]
print(XrayTube)

print('Done!')