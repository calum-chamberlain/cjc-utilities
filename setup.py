from setuptools import setup

setup(
    name='cjc_utilities',
    version='0.1',
    description='Utility functions written by Calum',
    url='https://bitbucket.org/calum-chamberlain/utilities',
    author='Calum Chamberlain',
    author_email='calum.chamberlain@vuw.ac.nz',
    license='GPL',
    packages=[
        'cjc_utilities', 'cjc_utilities.animator', 'cjc_utilities.get_data',
        'cjc_utilities.sac2nordic', 'cjc_utilities.coordinates',
        'cjc_utilities.plot_event', 'cjc_utilities.io'],
    zip_safe=False
)
