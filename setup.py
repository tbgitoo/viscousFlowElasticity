from setuptools import setup


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='viscousFlowElasticity',
      version='1.0',
      description='Simulation of visoelastic flow, specifically regarding Bingham plastics with finite yield strains',
      url='https://github.com/tbgitoo/viscousFlowElasticity',
      author='Thomas Braschler',
      author_email='thomas.braschler@gmail.com',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='MIT',
      package_dir = {
            'viscousFlowElasticity': 'viscousFlowElasticity',
            'viscousFlowElasticity.simulation_1mL_minimal': 'viscousFlowElasticity/simulation_1mL_minimal',
            'viscousFlowElasticity.linearElasticity': 'viscousFlowElasticity/linearElasticity',
            'viscousFlowElasticity.smoothing': 'viscousFlowElasticity/smoothing',
            'viscousFlowElasticity.volumeChange': 'viscousFlowElasticity/volumeChange',
            'viscousFlowElasticity.geometry': 'viscousFlowElasticity/geometry'
            },
      packages=["viscousFlowElasticity","viscousFlowElasticity.simulation_1mL_minimal",
                "viscousFlowElasticity.linearElasticity","viscousFlowElasticity.smoothing",
                "viscousFlowElasticity.volumeChange", "viscousFlowElasticity.plasticDeformation",
                "viscousFlowElasticity.geometry"],
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