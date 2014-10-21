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
        'networkx==1.8.1',
        'numpy==1.8.0',
        'scipy==0.13.3',
        'pandas==0.13.1',
        'scikit-learn==0.14.1',
        'biopython==1.63'
      ],
)
