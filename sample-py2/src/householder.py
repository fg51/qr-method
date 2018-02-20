# -*- coding: utf-8 -*-
r"""file: householder.py
x = H \cdot y
y = H \cdot x

H = I - \frac{2 (x - y) (x - y) ^ T}{(abs(x - y) ^ 2}
"""

from typing import Iterator, List, Tuple

from math import sqrt

__all__ = ["householder"]

EPS = 1.0E-8


def householder(xss: List[List[float]], size: int) -> None:
    for i in range(size - 2):
        #  変換行列 H の構築
        u, is_eps = create_u(init_u(xss, i, size), i, size)
        if is_eps is True:
            continue

        #  similarity transform
        d: List[float] = list(create_d(xss, u, i, size))
        ds: List[float] = list(create_ds(xss, u, i, size))
        update_hessenberg(xss, u, d, ds, size)


def init_u(xss: List[List[float]], k: int, size: int) -> Iterator[float]:
    for i in range(k + 1):
        yield 0.0
    for i in range(k + 1, size):
        yield xss[i][k]


def create_u(iter_u: Iterator[float], index: int, size: int) -> Tuple[List[float], bool]:
    u = list(iter_u)
    total = cal_sum_pow(u, index)
    next_index = index + 1
    if abs(u[next_index]) < EPS:
        return u, True
    sigma = sqrt(total) * u[next_index] / abs(u[next_index])
    u[next_index] += sigma
    norm = sqrt(2.0 * sigma * u[next_index])
    for i in range(next_index, size):
        u[i] /= norm
    return u, False


def cal_sum_pow(xs: List[float], begin: int) -> float:
    return sum(i * i for i in xs[begin + 1:])


def create_d(xss: List[List[float]], u: List[float], i: int, length: int) -> Iterator[float]:
    for v in update_d(init_d(xss, u, i), u, i + 1, length):
        yield v


def init_d(xss: List[List[float]], u: List[float], begin: int) -> Iterator[float]:
    for xs in xss:
        yield sum(xj * uj for xj, uj in zip(xs[begin + 1:], u[begin + 1:]))


def create_ds(xss: List[List[float]], u: List[float], i: int, length: int) -> Iterator[float]:
    for v in update_d(init_ds(xss, u, i, length), u, i + 1, length):
        yield v


def init_ds(xss: List[List[float]], u: List[float], begin: int, length: int) -> Iterator[float]:
    for i in range(length):
        yield sum(xj[i] * uj for xj, uj in zip(xss[begin + 1:], u[begin + 1:]))


def update_d(iter_ys: Iterator[float], xs: List[float], begin: int, end: int) -> Iterator[float]:
    ys = list(iter_ys)
    ud  = product_sum(xs, ys, begin, end)
    for i in range(end):
        yield 2.0 * (ys[i] - ud * xs[i])


def product_sum(xs: List[float], ys: List[float], begin: int, end: int) -> float:
    return sum(i * j for i, j in zip(xs[begin: end], ys[begin: end]))


def update_hessenberg(
        xss: List[List[float]],
        u:  List[float],
        d:  List[float],
        ds: List[float],
        num: int
) -> None:
    for i in range(num):
        for j in range(num):
            xss[i][j] -= u[i] * ds[j] + d[i] * u[j]
