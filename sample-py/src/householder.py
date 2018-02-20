# -*- coding: utf-8 -*-

from typing import List

from math import sqrt

__all__ = ["householder"]

EPS = 1.0E-8


def householder(xs: List[float], num: int) -> bool:
    u  = [0.0 for i in range(num)]
    d  = [0.0 for i in range(num)]
    ds = [0.0 for i in range(num)]

    for i in range(num - 2):
        #  変換行列 H の構築
        if create_u(u, xs, i, num) is True:
            continue

        #  similarity transform
        init_d(d, xs, u, i, num)
        update_d(d, u, i + 1, num)

        init_ds(ds, xs, u, i, num)
        update_d(ds, u, i + 1, num)

        update_hessenberg(xs, u, d, ds, num)
    return False


def create_u(u: List[float], xs: List[float], k: int, num: int) -> bool:
    init_u(u, xs, k, num)

    total = cal_sum(u, k, num)
    if abs(u[k + 1]) < EPS:
        return True

    sigma = sqrt(total) * u[k + 1] / abs(u[k + 1])
    u[k + 1] += sigma
    norm = sqrt(2.0 * sigma * u[k + 1])
    for i in range(k + 1, num):
        u[i] /= norm
    return False


def init_u(ys: List[float], xs: List[float], k: int, num: int) -> None:
    for i in range(k + 1):
        ys[i] = 0.0
    for i in range(k + 1, num):
        ys[i] = xs[num * i + k]


def cal_sum(xs: List[float], begin: int, num: int) -> float:
    return sum(i * i for i in xs[begin + 1: num])
    # total = 0.0
    # for (int i = begin + 1; i < num; i++) {
    #     total += xs[i] * xs[i]
    # return total


def init_d(d: List[float], xs: List[float], u: List[float], k: int, length: int) -> None:
    for i in range(length):
        d[i] = 0.0
        for j in range(k + 1, length):
            d[i] += xs[length * i + j] * u[j]


def init_ds(ds: List[float], xs: List[float], u: List[float], k: int, length: int) -> None:
    for i in range(length):
        ds[i] = 0.0
        for j in range(k + 1, length):
            ds[i] += xs[length * j + i] * u[j]


def product_sum(xs: List[float], ys: List[float], begin: int, end: int) -> float:
    return sum(i * j for i, j in zip(xs[begin: end], ys[begin: end]))


def update_d(ys: List[float], xs: List[float], begin: int, end: int) -> None:
    ud  = product_sum(xs, ys, begin, end)
    for i in range(end):
        ys[i] = 2.0 * (ys[i] - ud * xs[i])


def update_hessenberg(
        xs: List[float],
        u:  List[float],
        d:  List[float],
        ds: List[float],
        num: int
) -> None:
    for i in range(num):
        for j in range(num):
            xs[num * i + j] -= u[i] * ds[j] + d[i] * u[j]
