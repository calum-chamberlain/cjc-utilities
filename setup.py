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
        'cjc_utilities.coordinates', 'cjc_utilities.picker',
        'cjc_utilities.plot_event', 'cjc_utilities.io', 
        'cjc_utilities.magnitude_inversion'],
    zip_safe=False,
    scripts=["cjc_utilities/plot_event/plot_event.py", "cjc_utilities/picker/adjust_picks.py"]
)
