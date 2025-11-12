import pytest


def pytest_addoption(parser):
    """Add options to pytest.

    '--no-plot' option allows inhibiting plotting using the `no_plot` fixture (see below).
    """
    parser.addoption(
        "--no-plot", action="store_true",
        help="Do not produce and save plots (to speed up)"
    )


@pytest.fixture
def no_plot(request):
    return request.config.getoption("--no-plot")