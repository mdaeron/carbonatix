all: lock version

lock:
	uv lock
	cp pyproject.toml src/carbonatix/assets/pyproject.toml
	cp uv.lock src/carbonatix/assets/uv.lock

version:
	uv run python update-version.py
