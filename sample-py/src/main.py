# -*- coding: utf-8 -*-
from typing import Iterator, List
from householder import householder
from qr_decomp import qr_decomp

EPS = 1.0E-8
N = 4


def main():
    xs = [
        0.2, .32, .12, .3,
        0.1, .15, .24, .32,
        0.2, .24, .46, .36,
        0.6, .4,  .32, .2,
    ]
    print(" matrix ")
    for i in show_matrix(xs, N, N):
        print(i)

    is_error = householder(xs, N)
    if is_error:
        return

    print(" after house holder")
    for i in show_matrix(xs, N, N):
        print(i)

    qr_decomp(xs, N)
    print(" after QR method（対角項が固有値）")
    for i in show_matrix(xs, N, N):
        print(i)



def show_matrix(xs: List[float], n: int, m: int) -> Iterator:
    for i in range(n):
        yield " ".join("{0:10.6f}".format(xs[i * m + j]) for j in range(m))


if __name__ == '__main__':
    main()
