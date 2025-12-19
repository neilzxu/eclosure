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
            "r/eclosure/src/closedebh.cpp",
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
    name="eclosure",
    version="0.1.0",
    description="Python bindings and utilities for closed e-BH and closed BY procedures",
    author="BNP multiple testing team",
    license="GPL-3.0-only",
    python_requires=">=3.9",
    packages=find_packages(where="python/src"),
    package_dir={"": "python/src"},
    include_package_data=True,
    install_requires=[
        "numpy>=1.22",
        "pandas>=1.5",
        "scipy>=1.9",
        "matplotlib>=3.6",
        "seaborn>=0.12",
        "tqdm>=4.65",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pybind11>=2.11",
        ],
    },
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
)
