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
    url="https://github.com/xxxxxxxxx",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
    entry_points={
        'console_scripts': [
            'rec=usvcam.recorder:main', 
            'config=usvcam.open_config:main', 
            ]
    },
    data_files=[('etc/usvcam', ['usvcam/config.yaml', 'usvcam/D.mat'])] 
)
