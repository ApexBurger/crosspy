import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='crosspy',
    version="0.1.0",
    author="Bergsmo & McAuliffe",
    author_email="",
    description="A package for Fourier-space cross correlation of images",
    long_description=long_description,
    long_description_content_type="",
    url="",
    packages=['crosspy'],
    install_requires=['matplotlib','scipy','numpy','pillow','h5py','pathlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache 2.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)