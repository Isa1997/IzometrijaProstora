import sys

import numpy as np
import math
import  copy



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

def mulpliplyVectors(A,B):
    x = A[1] * B[2] - A[2] * B[1]
    y = -A[0] * B[2] + A[2] * B[0]
    z = A[0] * B[1] - A[1] * B[0]
    return np.array([x,y,z])

def AxisAngle(A):

    if round(np.linalg.det(A),0) != 1:
        print(np.linalg.det(A))
        print("Determinant not equal to 1!")
        sys.exit(1)

    Atr = np.transpose(A)
    if np.array_equal(Atr.dot(A), np.identity(3)):
        print("Matrix is not orthogonal!")
        sys.exit(1)
    Apom = copy.deepcopy(A)
    Apom[0][0] -= 1
    Apom[1][1] -= 1
    Apom[2][2] -= 1

    ii = -1
    jj = -1
    px = 0
    for i in range(2):
        for j in range(i+1,3):
            px = mulpliplyVectors(Apom[i],Apom[j])
            if not np.array_equal(np.array([0,0,0]),px):
                ii = i
                jj = j
                break

        if(ii != -1):
            break

    if ii == -1 and jj == -1:
        print("All rows are dependent")
        sys.exit(1)

    koef = math.sqrt(px[0] ** 2 + px[1] ** 2 + px[2] ** 2)

    for i in range(3):
        px[i] /= koef


    uV = Apom[ii]
    uVPrim = A.dot(uV)
    phi = math.acos((uV[0]*uVPrim[0]+uV[1]*uVPrim[1]+uV[2]*uVPrim[2])/(math.sqrt(uV[0]**2+uV[1]**2+uV[2]**2)*(math.sqrt(uVPrim[0]**2+uVPrim[1]**2+uVPrim[2]**2))))

    print("(p, alpha) = ", end=" ")
    print("({}, {})".format(px, phi))
    return px, phi

def Rodrigez(px,phi):

    koef = math.sqrt(px[0] ** 2 + px[1] ** 2 + px[2] ** 2)
    for i in range(3):
        px[i] /= koef

    px = px.reshape(3,1)
    pt = px.transpose()
    ppt = px.dot(pt)
    cos = math.cos(phi) * (np.identity(3) - ppt)
    pcross =math.sin(phi) * (np.array( [[0,-px[2][0],px[1][0]],[px[2][0],0,-px[0][0]],[-px[1][0],px[0][0],0]]))
    Rp = ppt + cos + pcross

    print("Matrix Aprim: ")
    print(Rp)

    return Rp

def A2Euler(A):

    if round(np.linalg.det(A), 0) != 1:
        print(np.linalg.det(A))
        print("Determinant not equal to 1!")
        sys.exit(1)

    Atr = np.transpose(A)
    if np.array_equal(Atr.dot(A), np.identity(3)):
        print("Matrix is not orthogonal!")
        sys.exit(1)

    phi = 0
    teta = 0
    psi = 0

    if A[2][0] < 1:
        # jedinstveno resenje
        if A[2][0] > -1:
            psi = math.atan2(A[1][0], A[0][0])
            teta = math.asin(-A[2][0])
            phi = math.atan2(A[2][1], A[2][2])
        else:
            psi = math.atan2(A[2][1], A[1][1])
            teta = math.pi / 2
            phi = 0
    else:
        psi = math.atan2(A[2][1], A[1][1])
        teta = - math.pi / 2
        phi = 0

    print ("(phi, teta, psi) = ", end=" ")
    print(phi, teta, psi)
    return phi, teta, psi

def AxisAngle2Q(p,phi):
    w = math.cos(phi/2)

    koef = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    for i in range(3):
        p[i] /= koef
        p[i] *= math.sin(phi/2)

    q = np.array([p[0], p[1], p[2], w])

    print("q = {}".format(q))
    return q

def Q2AngleAxis(q):
    koef = math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3]**2)
    for i in range(4):
        q[i] /= koef

    if(q[3] < 0):
        q = [-q[i] for i in range(4)]
    phi = 2 * math.acos(q[3])
    if abs(q[3]) == 1:
        p = np.array([1,0,0])

    else:
        p = np.array([q[i] for i in range(3)])
        koef = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
        for i in range(3):
            p[i] /= koef

    print("(pprim, alphaprim) = ({} {})".format(p, phi))
    return p, phi


def main():
    print("Euler2A(-atan(1/4), -asin(8/9), atan(4)): ")
    A = Euler2A(-math.atan(1/4), -math.asin(8/9), math.atan(4))
    print("---------------------------------------------------------------------------------------------")
    print("AxisAngle(A): ")
    (p, alpha) = AxisAngle(A)
    print("---------------------------------------------------------------------------------------------")

    print("Rodrigez(p, alpha): ")
    Aprim = Rodrigez(p, alpha)
    print("---------------------------------------------------------------------------------------------")
    print("A2Euler(Aprim): ")
    (phi, teta, psi) = A2Euler(Aprim)
    print("---------------------------------------------------------------------------------------------")
    print("AxisAngle2Q(p, alpha)")
    q = AxisAngle2Q(p, alpha)
    print("---------------------------------------------------------------------------------------------")
    print("Q2AngleAxis(q): ")
    (pprim, alphaprim) = Q2AngleAxis(q)






if __name__ == "__main__":
    main()