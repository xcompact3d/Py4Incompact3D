import setuptools

REQUIREMENTS=[
    "numpy",
    "mpi4py",
    # "decomp-2d4py"
]

setuptools.setup(
    name = "py4incompact3d",
    author = "Xcompact3d project",
    packages = ["py4incompact3d",
                "py4incompact3d.parallel",
                "py4incompact3d.postprocess",
                "py4incompact3d.tools",
                "py4incompact3d.deriv"],
    install_requires = REQUIREMENTS
)
