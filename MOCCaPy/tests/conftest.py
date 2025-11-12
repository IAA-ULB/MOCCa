import pytest

print(f" <<< Executing {__file__} >>>")


def pytest_addoption(parser):
    parser.addoption(
        "--no-plot", action="store_true",
        help="Do not produce and save plots (to speed up)"
    )


@pytest.fixture
def no_plot(request):
    return request.config.getoption("--no-plot")