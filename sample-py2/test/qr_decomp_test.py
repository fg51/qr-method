# -*- coding: utf-8 -*-
import sys
sys.path.append("src")


def import_module():
    from qr_decomp import qr_decomp  # pylint: disable=import-error
    return qr_decomp


def test_qr_decomp():
    EPS = 1E-2
    N = 4
    xss = list(create_target())
    qr_decomp = import_module()
    qr_decomp(xss, N)
    for xs, ys in zip(xss, gen_expects()):
        for i, j in zip(xs, ys):
            if abs(i - j) >= EPS:
                assert False
    assert True


def create_target():
    for i in [
            [+0.200000, -0.368570, -0.184136,  -0.192484],
            [-0.640312, +0.551951, +0.516187,  +0.081127],
            [-0.000000, +0.495143, +0.204566,  -0.104097],
            [-0.000000, +0.000000, -0.072862,  +0.053483],
    ]:
        yield i


def gen_expects():
    for i in [
            [+1.138579, +0.181886, -0.160414, -0.007475],
            [-0.000000, -0.227957, +0.245937, -0.148040],
            [+0.000000, -0.000000, +0.180906, +0.107371],
            [+0.000000, +0.000000, -0.000000, -0.081527],
    ]:
        yield i
