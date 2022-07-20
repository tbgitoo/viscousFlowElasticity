from setuptools import setup
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='viscousFlowElasticity',
      version='1.0',
      description='Simulation of visoelastic flow in a Bingham platic',
      url='https://github.com/tbgitoo/viscousFlowElasticity',
      author='Thomas Braschler',
      author_email='thomas.braschler@gmail.com',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='MIT',
      packages=["simulation_1mL_minimal","viscousFlowElasticity","linearElasticity","smoothing", "volumeChange", "plasticDeformation","geometry"],
      install_requires=[],
      zip_safe=False,
      project_urls={
          "Bug Tracker": "https://github.com/tbgitoo/viscousFlowElasticity/issues"
      },
      classifiers=[
          "Development Status :: 4 - Beta",
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          "Operating System :: OS Independent",
      ],
    python_requires='>=3.5'
      )