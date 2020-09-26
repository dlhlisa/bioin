"""
bioin
Study notes for beginners bioinformatics
"""

# Add imports here
from .bioin import *
from . import motif
from . import replication


# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
