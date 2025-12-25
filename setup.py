import setuptools
import usvcam

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="usvcam",
    version=usvcam.__version__,
    author="Jumpei Matsumoto",
    author_email="jm@u-toyama",
    description="A system for localization and assignment of ultrasonic vocalizations in rodents",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/MatsumotoJ/usvcam",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    entry_points={
        'console_scripts': [
            'rec=usvcam.recorder_fadc:main', 
            'rec_legacy=usvcam.recorder_legacy:main', 
            'config=usvcam.open_config:main', 
            'config_legacy=usvcam.open_config:main_legacy', 
            ]
    },
    package_data={"usvcam": ["*.yaml", "*.dll",]}
)
