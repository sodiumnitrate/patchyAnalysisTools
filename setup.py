from skbuild import setup

setup(
    name="patchyAnalysisTools",
    version="2.0.0",
    packages=["patchyAnalysisTools"],
    package_dir={"": "src"},
    cmake_install_dir="src/patchyAnalysisTools",
)