lock:
	cp pyproject.toml src/carbonatix/assets/pyproject.toml
	cp uv.lock src/carbonatix/assets/uv.lock

update-version:
	uv run python update-version.py

