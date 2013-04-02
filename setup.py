###############################################################################
#					Setup File for Mixed Membership Block Model
#					Currently Generic Template
#					To Do: Get proper pip, virtualenv, distribute in conjunctio
#					n with git
###############################################################################
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'My Project',
    'author': 'My Name',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'My email.',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['NAME'],
    'scripts': [],
    'name': 'projectname'
}

setup(**config)
