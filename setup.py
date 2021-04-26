from setuptools import setup

setup(
    name="twirl",
    version="0.0.1",
    author="Lionel J. Garcia",
    description="python-only astrometric plate solving",
    py_modules=["twirl"],
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
