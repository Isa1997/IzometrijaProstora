
import numpy as np
import math


def Euler2A(xa,ya,za):
    rx = [[1,0,0],[0,math.cos(xa),-math.sin(xa)],[0,math.sin(xa),math.cos(xa)]]
    Rx = np.array(rx)

    ry = [[math.cos(ya),0,math.sin(ya)],[0,1,0],[-math.sin(ya),0,math.cos(ya)]]
    Ry = np.array(ry)

    rz = [[math.cos(za),-math.sin(za),0],[math.sin(za),math.cos(za),0],[0,0,1]]
    Rz = np.array(rz)


    matrixA = Rz.dot(Ry.dot(Rx))
    print("Matrix A:")
    print(matrixA)
    return matrixA


def main():
    print("ss")
    Euler2A(-math.atan(1/4),-math.asin(8/9),math.atan(4))



if __name__ == "__main__":
    main()