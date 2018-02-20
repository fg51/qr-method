# -*- coding: utf-8 -*-
from typing import Iterator, List
from householder import householder
from qr_decomp import qr_decomp

EPS = 1.0E-8
N = 4


def main():
    xss = [
        [0.20, .32, .12, .30],
        [0.10, .15, .24, .32],
        [0.20, .24, .46, .36],
        [0.60, .40, .32, .20],
    ]
    print(" matrix ")
    for i in show_matrix(xss):
        print(i)

    householder(xss, N)
    print(" after house holder")
    for i in show_matrix(xss):
        print(i)

    qr_decomp(xss, N)

    print(" after QR method（対角項が固有値）")
    for i in show_matrix(xss):
        print(i)



def show_matrix(xss: List[List[float]]) -> Iterator:
    for xs in xss:
        yield " ".join("{0:10.6f}".format(j) for j in xs)


if __name__ == '__main__':
    main()
