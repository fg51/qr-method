# -*- coding: utf-8 -*-
from typing import List, Tuple

from math import sqrt

__all__ = ["qr_decomp"]

G_EPS = 1.0E-8


def set_EPS(eps):
    global G_EPS  # pylint: disable=global-statement
    G_EPS = eps


def qr_decomp(xss: List[List[float]], size: int) -> None:
    m = size
    while True:
        if m == 1:
            return
        if abs(xss[m - 1][m - 2]) < G_EPS:
            m -= 1
            continue

        mu = cal_mu(xss, m)  # mu is shift-eigen
        for i in range(m):
            xss[i][i] -= mu

        qr_method(xss, m, size)
        for i in range(m):
            xss[i][i] += mu


def cal_mu(xss: List[List[float]], m: int) -> float:
    pos1, pos2 = m - 2, m - 1
    a00, a01 = xss[pos1][pos1], xss[pos1][pos2]
    a10, a11 = xss[pos2][pos1], xss[pos2][pos2]
    yss = [
        [xss[pos1][pos1], xss[pos1][pos2]],
        [xss[pos2][pos1], xss[pos2][pos2]],
    ]
    wa = cal_wa(yss)

    x, y = a00 + a11, a00 * a11 - a01 * a10
    lam1 = 0.5 * (x + wa)
    lam2 = y / lam1
    if abs(a11 - lam1) < abs(a11 - lam2):
        return a11 - lam1
    return a11 - lam2


# def cal_wa(x00: float, x01: float, x10: float, x11: float) -> float:
#     a, b = x00 + x11, x00 * x11 - x01 * x10
#     wa = a * a - 4.0 * b
#     return 0.0 if wa < 0.0 else sqrt(wa)

def cal_wa(xss: List[List[float]]) -> float:
    a = xss[0][0] + xss[1][1]
    b = xss[0][0] * xss[1][1] - xss[0][1] * xss[1][0]
    wa = a * a - 4.0 * b
    return 0.0 if wa < 0.0 else sqrt(wa)


def qr_method(a: List[List[float]], m: int, n: int) -> None:
    q = create_matrix_I(m)

    for i in range(m - 1):
        norm = to_norm([a[i][i], a[i + 1][i]])
        cosx, sinx = cal_rotation(norm, a[i][i], a[i + 1][i])
        update_a(a, sinx, cosx, i, m, n)
        a[i][i] = norm

        update_q(q, m, i, sinx, cosx)
    for i in range(m):
        w = [a[i][j] for j in range(m)]
        for j in range(m):
            a[i][j] = cal_sum_w_q(w, q, i, m, j)


def create_matrix_I(m: int) -> List[float]:
    xs = [0.0 for i in range(m * m)]
    for i in range(m):
        xs[m * i + i] = 1.0
    return xs


def to_norm(xs: List[float]) -> float:
    return sqrt(sum(i * i for i in xs))


def cal_rotation(norm: float, x: float, y: float) -> Tuple[float, float]:
    cosx = 0.0 if abs(norm) < G_EPS else x / norm
    sinx = 0.0 if abs(norm) < G_EPS else y / norm
    return cosx, sinx


def update_a(a: List[List[float]], sinx: float, cosx: float, begin: int, end: int, _: int) -> None:
    """
    rotation ?
    """
    for i in range(begin + 1, end):
        x, y = a[begin][i], a[begin + 1][i]
        a[begin][i] = x * cosx + y * sinx
        a[begin + 1][i] = -x * sinx + y * cosx
    a[begin + 1][begin] = 0.0


def update_q(q: List[float], end: int, i: int, sinx: float, cosx: float) -> None:
    """
    rotation ?
    """
    for j in range(0, end * end, end):
        eji = j + i
        a, b = q[eji], q[eji + 1]
        q[eji] = a * cosx + b * sinx
        q[eji + 1] = -a * sinx + b * cosx


def cal_sum_w_q(w: List[float], q: List[float], begin: int, end: int, j: int) -> float:
    return sum(w[i] * q[end * i + j] for i in range(begin, end))
