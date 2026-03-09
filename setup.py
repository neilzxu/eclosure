from pathlib import Path

import pybind11
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class BuildExt(build_ext):
    c_opts = {
        "msvc": ["/EHsc", "/std:c++17"],
        "unix": ["-std=c++17", "-O3"],
    }
    l_opts = {"msvc": [], "unix": []}

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == "unix":
            opts.append("-fvisibility=hidden")
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = self.l_opts.get(ct, [])
        super().build_extensions()


ROOT = Path(__file__).resolve().parent
ext_modules = [
    Extension(
        "eclosure._cebh",
        sources=[
            "python/src/eclosure/_cebh_module.cpp",
            "r/eclosure/src/ceBH.cpp",
            "r/eclosure/src/cBY.cpp",
        ],
        include_dirs=[
            pybind11.get_include(),
            pybind11.get_include(user=True),
            str(ROOT / "r" / "eclosure" / "src"),
        ],
        define_macros=[("E_CLOSURE_NO_RCPP", None)],
        language="c++",
    )
]

setup(
    packages=find_packages(where="python/src"),
    package_dir={"": "python/src"},
    include_package_data=True,
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
)
