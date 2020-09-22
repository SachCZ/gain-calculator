import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gain-calculator",
    version="0.0.5",
    author="Martin Sach",
    author_email="martin.sachin@gmail.com",
    description="A package providing a wrapper around Flexible Atomic Code specifically for gain predictions.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SachCZ/gain-calculator",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
    install_requires=[
        "numpy",
        "ray",
        "typing"
    ]
)
