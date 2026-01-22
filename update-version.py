import tomllib

with open('pyproject.toml', 'rb') as fid:
	v = tomllib.load(fid)['project']['version']

with open('src/carbonatix/__version__.py', 'w') as fid:
	fid.write(f'__version__ = "{v}"')
