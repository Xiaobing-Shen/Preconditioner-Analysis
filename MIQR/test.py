# -*- coding: utf-8 -*-
# @Time    : 18-4-10 下午9:23
# @Author  : Icy Shen
# @Email   : SAH1949@126.com
from pcgnr.python3 import IC
import numpy as np
a = np.random.random(size = (5, 5))
b = np.random.random(size=(5, ))
temp = IC.preconIMQR(A = a, b = b, o = 1)
print(temp)