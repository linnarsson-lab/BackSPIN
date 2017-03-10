from setuptools import setup, find_packages

__version__ = "0.2.0"

setup(
    name="backspinpy",
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'future',
        'numpy',
        'scipy',
        'scikit-learn'
    ],
    # command
    scripts=['backspinpy/backspin'],
    # metadata
    author="Linnarsson Lab",
    author_email="gioelelamanno@gmail.com",
    description="backSPIN clustering algorythm",
    license="MIT",
    url="https://github.com/linnarsson-lab/BackSPIN",
)
