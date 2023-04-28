from setuptools import find_packages, setup

VERSION = "1.0.0"

setup(name='freewillson', version=VERSION,
      description='Free Wilson analysis using the RDKit and Scikit Learn', author='Brian Kelley',
      author_email='fustigator@gmail.com', install_requires=['scikit-learn', 'tqdm'],
      packages=find_packages(), py_modules=['freewilson'])
