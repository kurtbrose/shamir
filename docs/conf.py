import os
import sys
sys.path.insert(0, os.path.abspath('..'))

project = 'shamir'
author = 'Kurt Rose'

# The short X.Y version
version = '17.12.0'
# The full version, including alpha/beta/rc tags
release = '17.12.0'

extensions = ['sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'alabaster'
html_static_path = ['_static']
