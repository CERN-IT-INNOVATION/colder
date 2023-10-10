from importlib.metadata import metadata

PACKAGE = "colder"

__version__ = metadata(PACKAGE)["version"]