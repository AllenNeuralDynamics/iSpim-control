import setuptools

with open("README.md", "r", newline="", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mesospim",
    version="0.0.1",
    author="Joshua Vasquez",
    author_email="joshua.vasquez@alleninstitute.org",
    description="mesospim prime abstraction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    license="MIT",
    keywords=['microscope', 'spim'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.6',
    # TODO: figure out how to make pip skip over already-installed packages
    #   that come from github rather than reinstalling from the repo every time.
    install_requires=['TigerASI @ git+ssh://git@github.com/AllenNeuralDynamics/TigerASI@main',
                      'obis_laser @ git+ssh://git@github.com/AllenNeuralDynamics/obis-laser@main',
                      'opto @ git+ssh://git@github.com/AllenNeuralDynamics/opto@master',
                      'scipy', 'numpy', 'matplotlib', 'tomlkit', 'mock',
                      'coloredlogs', 'nidaqmx' ]
)
