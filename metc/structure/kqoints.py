# -*- coding: utf-8 -*-
import sys
import numpy as np

def get_kqlist(args):
    numargs = len(args)
    kqlist = []
    if numargs < 3 or numargs > 4:
        print("usage: n1 n2 n3 [wan]")
        print("       n1  - divisions along 1st recip vector")
        print("       n2  - divisions along 2nd recip vector")
        print("       n3  - divisions along 3rd recip vector")
        print("       wan - omit the kpoint weight (optional)")
        sys.exit()

    # 转换为整数
    try:
        n1 = int(args[0])
        n2 = int(args[1])
        n3 = int(args[2])
    except ValueError:
        print("n1, n2, and n3 must be integers")
        sys.exit()

    if n1 <= 0:
        print("n1 must be >0")
        sys.exit()
    if n2 <= 0:
        print("n2 must be >0")
        sys.exit()
    if n3 <= 0:
        print("n3 must be >0")
        sys.exit()

    total_points = n1 * n2 * n3

    if numargs == 3:
        print("K_POINTS crystal")
        print(total_points)
        for x in range(n1):
            for y in range(n2):
                for z in range(n3):
                    kx = x / n1
                    ky = y / n2
                    kz = z / n3
                    weight = 1 / total_points
                    # print(f"{kx:12.8f}{ky:12.8f}{kz:12.8f}{weight:14.6e}")
                    kqlist.append(np.array([kx, ky, kz, weight]))
    else:
        for x in range(n1):
            for y in range(n2):
                for z in range(n3):
                    kx = x / n1
                    ky = y / n2
                    kz = z / n3
                    kqlist.append(np.array([kx, ky, kz, 1.0]))
                    # print(f"{kx:12.8f}{ky:12.8f}{kz:12.8f}")
    return np.array(kqlist)

class kqmesh:
    def __init__(self):
        pass

if __name__== '__main__':
    args = sys.argv[1:]
    kqlist = get_kqlist(args)
    print(kqlist)


