import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="twirl",
    version="0.2.0",
    author="Lionel J. Garcia",
    description="twirl is an astrometric plate solving package for Python.",
    long_description=README,
    long_description_content_type="text/markdown",
    packages=["twirl"],
    license="MIT",
    url="https://github.com/lgrcia/twirl",
    install_requires=["numpy", "astropy>=4.3", "astroquery"],
    extras_require={
        "docs": [
            "sphinx",
            "docutils",
            "jupyterlab",
            "myst-parser",
            "twine",
            "sphinx-book-theme",
            "black",
            "myst-nb",
            "sphinx-copybutton",
        ]
    },
    zip_safe=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
