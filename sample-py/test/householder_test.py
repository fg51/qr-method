# -*- coding: utf-8 -*-
import sys
sys.path.append("src")


def import_module():
    from householder import householder  # pylint: disable=import-error
    return householder


def test_qr_decomp():
    EPS = 1E-2
    N = 4
    xs = list(create_target())
    householder = import_module()
    householder(xs, N)
    for i, j in zip(xs, gen_expects()):
        if abs(i - j) >= EPS:
            assert False
    assert True


def create_target():
    for i in [
            +0.200000, +0.320000, +0.120000, +0.300000,
            +0.100000, +0.150000, +0.240000, +0.320000,
            +0.200000, +0.240000, +0.460000, +0.360000,
            +0.600000, +0.400000, +0.320000, +0.200000,
    ]:
        yield i


def gen_expects():
    for i in [
            +0.200000, -0.368570, -0.184136,  -0.192484,
            -0.640312, +0.551951, +0.516187,  +0.081127,
            -0.000000, +0.495143, +0.204566,  -0.104097,
            -0.000000, +0.000000, -0.072862,  +0.053483,
    ]:
        yield i
