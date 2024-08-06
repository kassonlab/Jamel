#200nm
import math
import random
import numpy as np
import pandas as pd


def sphere_volume(radius):
    return 4*math.pi*radius**3/3

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
if __name__=='__main__':
    print(pd.read_csv('sphere_packing',sep='\t')['ratio']*100)
    # contian=SphereContainer()
    # for x in range(3):
    #     contian.add_sphere(Sphere(50))
    # print(contian.spheres)
