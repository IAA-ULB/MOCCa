import pytest


def pytest_addoption(parser):
    """Add options to pytest.

    '--no-plot' option allows inhibiting plotting using the `no_plot` fixture (see below).
    """
    parser.addoption(
        "--no-plot", action="store_true",
        help="Do not produce and save plots (to speed up)"
    )

    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def no_plot(request):
    return request.config.getoption("--no-plot")