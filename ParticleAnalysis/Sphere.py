#200nm
import math
import random
import numpy as np
import pandas as pd

def sphere_packing_ratios(monomer_diameter) -> pd.Series:
    return pd.read_csv('sphere_packing',sep='\t')['ratio']*monomer_diameter
def sphere_volume(radius):
    return 4*math.pi*radius**3/3


def find_sphere_packing_index(packing_ratios:pd.Series,target_diameter):
    for index,diameter in packing_ratios.items():
        if target_diameter<diameter:
            return index-1


class Sphere:
    def __init__(self,radius):
        self.radius=radius
        self.coordinates=np.array((random.randint(-1000,1000),random.randint(-1000,1000),random.randint(-1000,1000)))

class SphereContainer:
    def __init__(self):
        self.spheres=[]
    def add_sphere(self,sphere:Sphere):
        print(sphere)
        self.spheres.append(sphere)
    def find_center(self):
        pass
# if __name__=='__main__':
    # contian=SphereContainer()
    # for x in range(3):
    #     contian.add_sphere(Sphere(50))
    # print(contian.spheres)
