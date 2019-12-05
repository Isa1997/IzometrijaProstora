import sys

import numpy as np
import math
import copy


def Euler2A(xa, ya, za):
    rx = [[1, 0, 0], [0, math.cos(xa), -math.sin(xa)], [0, math.sin(xa), math.cos(xa)]]
    Rx = np.array(rx)

    ry = [[math.cos(ya), 0, math.sin(ya)], [0, 1, 0], [-math.sin(ya), 0, math.cos(ya)]]
    Ry = np.array(ry)

    rz = [[math.cos(za), -math.sin(za), 0], [math.sin(za), math.cos(za), 0], [0, 0, 1]]
    Rz = np.array(rz)

    matrixA = Rz.dot(Ry.dot(Rx))
    print("Matrix A:")
    matrixA = np.round(matrixA, 6)
    print(matrixA)
    return matrixA


def mulpliplyVectors(A, B):
    x = A[1] * B[2] - A[2] * B[1]
    y = -A[0] * B[2] + A[2] * B[0]
    z = A[0] * B[1] - A[1] * B[0]
    return np.array([x, y, z])


def AxisAngle(A):
    A = A.astype(float)
    if round(np.linalg.det(A), 3) != 1.0:
        print("Determinant not equal to 1!")
        print(np.linalg.det(A))
        sys.exit(1)

    Atr = np.transpose(A)
    product = Atr.dot(A)
    product = product.astype(float)
    product = np.round(product, 2)
    if not np.array_equal(product, np.identity(3)):
        print(product)
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
        for j in range(i + 1, 3):
            px = mulpliplyVectors(Apom[i], Apom[j])
            if not np.array_equal(np.array([0, 0, 0]), px):
                ii = i
                jj = j
                break

        if (ii != -1):
            break

    if ii == -1 and jj == -1:
        print("All rows are dependent")
        sys.exit(1)

    koef = math.sqrt(px[0] ** 2 + px[1] ** 2 + px[2] ** 2)

    for i in range(3):
        px[i] /= koef

    uV = Apom[ii]
    uVPrim = A.dot(uV)
    phi = math.acos((uV[0] * uVPrim[0] + uV[1] * uVPrim[1] + uV[2] * uVPrim[2]) / (
            math.sqrt(uV[0] ** 2 + uV[1] ** 2 + uV[2] ** 2) * (
        math.sqrt(uVPrim[0] ** 2 + uVPrim[1] ** 2 + uVPrim[2] ** 2))))

    # Provera mesovitog proivoda [u, up, p] < 0 => p = -p
    mesoviti = uV[0] * uVPrim[1] * px[2] + uVPrim[0] * px[1] * uV[2] + px[0] * uV[1] * uVPrim[2] - uVPrim[0] * uV[1] * \
               px[2] - uV[0] * px[1] * uVPrim[2] - px[0] * uVPrim[1] * uV[2]

    if uV[0]**2 + uV[1]**2 + uV[2]**2 == uVPrim[0]**2 + uVPrim[1]**2 + uVPrim[2]**2: # da ga ne okrece za dzabe jer nije bitno
        px = 1*px
    elif mesoviti < 0:
        px = (-1) * px

    print("(p, alpha) = ", end=" ")
    print("({}, {})".format(px, phi))
    return px, phi


def Rodrigez(px, phi):
    koef = math.sqrt(px[0] ** 2 + px[1] ** 2 + px[2] ** 2)
    for i in range(3):
        px[i] /= koef

    px = px.reshape(3, 1)
    pt = px.transpose()
    ppt = px.dot(pt)
    cos = math.cos(phi) * (np.identity(3) - ppt)
    pcross = math.sin(phi) * (np.array([[0, -px[2][0], px[1][0]], [px[2][0], 0, -px[0][0]], [-px[1][0], px[0][0], 0]]))
    Rp = ppt + cos + pcross

    print("Matrix Aprim: ")
    Rp = np.round(Rp, 6)
    print(Rp)

    return Rp


def A2Euler(A):
    A = A.astype(float)
    if round(np.linalg.det(A), 3) != 1.0:
        print(np.linalg.det(A))
        print("Determinant not equal to 1!")
        sys.exit(1)

    Atr = np.transpose(A)
    product = Atr.dot(A)
    product = product.astype(float)
    product = np.round(product, 2)
    if not np.array_equal(product, np.identity(3)):
        print(product)
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

    print("(phi, teta, psi) = ", end=" ")
    print(phi, teta, psi)
    return phi, teta, psi


def AxisAngle2Q(p, phi):
    w = math.cos(phi / 2)

    pp = [1.0, 1.0, 1.0, 1.0]

    koef = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    for i in range(3):
        p[i] /= koef

    for i in range(3):
        pp[i] = p[i] * math.sin(phi / 2) * 1.0

    q = np.array([pp[0], pp[1], pp[2], w])

    print("q = {}".format(q))
    return q


def Q2AngleAxis(q):
    koef = math.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
    for i in range(4):
        q[i] /= koef

    if (q[3] < 0):
        q = [-q[i] for i in range(4)]
    phi = 2 * math.acos(q[3])
    if abs(q[3]) == 1:
        p = np.array([1, 0, 0])

    else:
        p = np.array([q[i] for i in range(3)])
        koef = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
        for i in range(3):
            p[i] /= koef

    print("(pprim, alphaprim) = ({} {})".format(p, phi))
    return p, phi


def main():
    print("Euler2A(-atan(1/4), -asin(8/9), atan(4)): ")
    A = Euler2A(-math.atan(1 / 4), -math.asin(8 / 9), math.atan(4))
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
    A1 = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    print("A1:")
    # print(A1)
    # (phi1, teta1, psi1) = A2Euler(A1)
    #
    # A2 = np.array([[3/4,1/4,math.sqrt(6)/4],[1/4,3/4,-math.sqrt(6)/4],[-math.sqrt(6)/4,math.sqrt(6)/4,2/4]])
    # print("A2")
    # print(A2)
    # (phi2, teta2, psi2) = A2Euler(A2)

    # q = AxisAngle2Q(np.array([1,0,0]), math.pi/2)
    # (p, phi) = Q2AngleAxis(q)

    # D = Euler2A(phi1, teta1, psi1)

    p = np.array([0, 3, 0])
    alpha = math.pi / 4
    print(p, alpha)
    #
    Rp = Rodrigez(p, alpha)
    #
    (p, alpha) = AxisAngle(Rp)

    # 1) Matrica rotacije oko x,pa y, pa z ose za ugao 10*pi/6 => RADI
    ugao = 10*math.pi/6
    print("Oko x:")
    okoX = Euler2A(ugao,0,0)
    print("Oko y:")
    okoY = Euler2A(0,ugao,0)
    print("Oko z:")
    okoZ = Euler2A(0,0,ugao)

    # 2) Rotacija oko z za 3pi/2 i x za pi/2 i rotacija oko y za pi/2 i x za pi- sopstvena rotacija => RADI
    okoXZ = Euler2A(math.pi/2,0,3*math.pi/2)
    okoXY = Euler2A(math.pi, math.pi/2, 0)

    # 4) Rodrigez za dato p i ugap => RADI
    p = np.array([math.sqrt(3)/3, math.sqrt(3)/3, math.sqrt(3)/3])
    Rp = Rodrigez(p, math.pi*2/3)

    # 5) Rodrigez opet
    p = np.array([1/3, 2/3, 2/3])
    Rp = Rodrigez(p, math.pi * 3 / 2)
    p = np.array([math.sqrt(2) / 2, math.sqrt(2) / 2, 0])
    Rp = Rodrigez(p, math.pi / 3)

    # 6) Izvlacenje vektora i ugla iz date matrice => RADI
    A = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    (p, alpha) = AxisAngle(A)
    A = np.array([[1 / 9, 8 / 9, -4 / 9], [-4 / 9, 4 / 9, 7 / 9], [8 / 9, 1 / 9, 4 / 9]])
    (p, alpha) = AxisAngle(A)
    A = np.array([[3 / 4, 1 / 4, math.sqrt(6) / 4], [1 / 4, 3 / 4, -math.sqrt(6) / 4], [-math.sqrt(6) / 4, math.sqrt(6) / 4, 2 / 4]])
    (p, alpha) = AxisAngle(A)

    # 7) Opet Rodgrigez => RADI
    p = np.array([-math.sqrt(2)/2, math.sqrt(2)/2, 0])
    ugao = math.pi/3
    Rp = Rodrigez(p, ugao)

    # 8) Izvlaci uglove iz matrice => RADI
    A = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    (phi, teta, psi) = A2Euler(A)
    A = np.array([[1 / 9, 8 / 9, -4 / 9], [-4 / 9, 4 / 9, 7 / 9], [8 / 9, 1 / 9, 4 / 9]])
    (phi, teta, psi) = A2Euler(A)



if __name__ == "__main__":
    main()
