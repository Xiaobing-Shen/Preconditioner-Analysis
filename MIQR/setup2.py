# -*- coding: utf-8 -*-
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
from Cython.Build import cythonize
import platform
import sys
version = int(platform.python_version()[0])
if version == 2:
    from imp import reload
    reload(sys)
    sys.setdefaultencoding('utf8')

ext = cythonize(Extension(name = 'imqrCG',sources = ["imqrCG.pyx","auxf.c",
	"imqr.c","matfun.c", "pcgnr.c"]
    ))
setup(
    name="imqrCG",
    url = "http://git-best.top/linear_system_solver/indirect.git",
    version="1.0.0",
    author="sah",
    author_email = 'SAH1949@126.com',
    description = 'preconditioned CG ',
    license="SUFE",
    requires=["numpy","scipy"],
    ext_modules=ext,
),
