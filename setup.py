import setuptools

# with open("README.md", "r") as fh:
#     long_description = fh.read()

setuptools.setup(
    name="qcweather_Florent-H",
    version="0.0.1",
    author="Florent Herbinger",
    author_email="florent.herbinger@outlook.com",
    description="A python package to automatically control the quality of yearly weather files",
    long_description="long_description",
    long_description_content_type="text/markdown",
    url="https://github.com/Florent-H/qcweather",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)