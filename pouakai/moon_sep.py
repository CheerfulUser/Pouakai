from core import pouakai
import pandas as pd

test = pd.read_csv('../dev/moon_seperation_test.csv')

files = test['file'].values

for i in range(len(files)-4):
    i += 4
    print(files[i])
    try:
        pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/moon_test/',time_tolerence=1000)
    except:
        print('failed image')
