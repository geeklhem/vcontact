from setuptools import setup

setup(name='vcontact',
      version=open('VERSION').read(),
      description='Viral Contig Automatic Clutering and Taxonomy',
      url='http://www.eleves.ens.fr/home/doulcier/projects/virus/',
      author='Guilhem Doulcier',
      long_description=open('README').read(),
      author_email='guilhem.doulcier@ens.fr',
      license='GPLv3',
      packages=['vcontact'],
      scripts=['bin/vcontact','bin/vcontact-pcs'],
      install_requires=[
        'networkx',
        'numpy',
        'scipy',
        'pandas',
        'scikit-learn'
      ],
)
