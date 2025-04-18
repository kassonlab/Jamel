from setuptools import setup, Extension
import pybind11
import test

ext_modules = [
    Extension(
        "chi_cpp",  # Module name
        ["/mnt/c/Users/jamel/PycharmProjects/Jamel/Chimeragenesis/Scripts/test.cpp"],  # Source files
        include_dirs=[pybind11.get_include()],  # Include pybind11 headers
        language="c++",  # Specify C++ as the language
    )
]

setup(
    name="chi_cpp",
    version="0.1",
    ext_modules=ext_modules,
)