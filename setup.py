import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="twirl",
    version="0.0.3",
    author="Lionel J. Garcia",
    description="python-only astrometric plate solving",
    long_description=README,
    long_description_content_type="text/markdown",
    py_modules=["twirl"],
    license="MIT",
    url="https://github.com/lgrcia/twirl",
    # entry_points="""
    #     [console_scripts]
    #     prose=main:cli
    # """,
    
    install_requires=[
        "numpy",
        "scipy",
        "astropy==4.0",
        "matplotlib",
        "scikit-image",
    ],
    zip_safe=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
