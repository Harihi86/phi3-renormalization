"""Keep this repository's test scratch space out of shared user temp histories."""

from uuid import uuid4


def pytest_configure(config):
    if config.option.basetemp is None:
        parent = config.rootpath / "runs" / "pytest-tmp"
        parent.mkdir(parents=True, exist_ok=True)
        config.option.basetemp = str(parent / uuid4().hex[:12])
