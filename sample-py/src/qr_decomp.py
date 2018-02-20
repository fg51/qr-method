# -*- coding: utf-8 -*-
from typing import List

from math import sqrt

__all__ = "qr_decomp"

EPS = 1.0E-8


def qr_decomp(a: List[float], n: int) -> None:
    q = [0.0 for i in range(n * n)]
    w = [0.0 for i in range(n)]

    m = n

    while m != 1:
        if abs(a[n * (m - 1) + m - 2]) < EPS:
            m -= 1
            continue

        # 原点移動 mu
        mu = cal_mu(
            a[n * (m - 2) + m - 2], a[n * (m - 2) + m - 1],
            a[n * (m - 1) + m - 2], a[n * (m - 1) + m - 1])
        for i in range(m):
            a[n * i + i] -= mu

        # QR method
        init_q(q, m)

        for i in range(m - 1):
            sum1 = sqrt(a[n * i + i] * a[n * i + i] + a[n * i + n + i] * a[n * i + n + i])
            sinx = 0.0 if abs(sum1) < EPS else a[n * i + n + i] / sum1
            cosx = 0.0 if abs(sum1) < EPS else a[n * i + i] / sum1

            update_a(a, sinx, cosx, i, m, n)
            a[n * i + i] = sum1
            update_q(q, m, i, sinx, cosx)

        for i in range(m):
            for j in range(m):
                w[j] = a[n * i + j]
            for j in range(m):
                a[n * i + j] = cal_sum_w_q(w, q, i, m, j)

        for i in range(m):
            a[n * i + i] += mu


def cal_wa(x00: float, x01: float, x10: float, x11: float) -> float:
    sum1 = x00 + x11
    sum2 = x00 * x11 - x01 * x10
    wa = sum1 * sum1 - 4.0 * sum2
    return 0.0 if wa < 0.0 else sqrt(wa)


def cal_mu(a00: float, a01: float, a10: float, a11: float) -> float:
    wa = cal_wa(a00, a01, a10, a11)

    sum1 = a00 + a11
    sum2 = a00 * a11 - a01 * a10
    lam1 = 0.5 * (sum1 + wa)
    lam2 = sum2 / lam1
    if abs(a11 - lam1) < abs(a11 - lam2):
        return a11 - lam1
    return a11 - lam2


def init_q(q: List[float], m: int) -> None:
    for i in range(m * m):
        q[i] = 0.0
    for i in range(m):
        q[m * i + i] = 1.0


def update_a(a: List[float], sinx: float, cosx: float, begin: int, end: int, n: int) -> None:
    for j in range(begin + 1, end):
        nij = n * begin + j
        sum2 = a[nij] * cosx + a[nij + n] * sinx
        a[nij + n] = -a[nij] * sinx + a[nij + n] * cosx
        a[nij] = sum2
    a[n * begin + n + begin] = 0.0


def update_q(q: List[float], end: int, i: int, sinx: float, cosx: float) -> None:
    for j in range(end):
        sum2 = q[end * j + i] * cosx + q[end * j + i + 1] * sinx
        q[end * j + i + 1] = -q[end * j + i] * sinx + q[end * j + i + 1] * cosx
        q[end * j + i] = sum2


def cal_sum_w_q(w: List[float], q: List[float], begin: int, end: int, j: int) -> float:
    return sum(w[i] * q[end * i + j] for i in range(begin, end))
