import json

Material = ['Si','Cu','Al','ABSresin','Polyester','Polystyrene']
Title = ['Photon Energy',
         'Coherent Scattering',
         'Incoherent Scattering',
         'Photoelectric Absorption',
         'Total with Coherent',
         'Total without Coherent']
#Energy = [
#    'E=' + str(i) + 'keV' for i in range(1,101)
#]

data = dict()
for material in Material:
    print('--- open file : ', str(material))

    f_txt = open('xray_mu_' + str(material) + '.txt', 'r')
    list = f_txt.readlines()


    for i in range(3):
        list.pop(0)
    dict_energy = dict()
    for i in range(len(list)):
        s = list[i]
        list_split = s.split()
        dict_data = dict(zip(Title,list_split))
        e_str = dict_data[Title[0]]
        e_float = float(e_str) * 10**3
        energy = '%.1f' % e_float
        Energy = 'E=' + str(energy) + 'keV'
        dict_energy[Energy] = dict_data

    data[material] = dict_energy

    f_txt.close()

f_json = open('xray_data.json', 'w')
json.dump(data,f_json,indent=4)
f_json.close()
print('Done!')