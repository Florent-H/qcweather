import pytest


class TestsDifferenceQA:
    """A class that encapsulates the different tests for the functions found in
    the sensitivity_analysis.py module
    """

    """
    ~~~ PYTEST FIXTURES ~~~
    """

    @pytest.fixture()
    def a(self):
        yield a

    """
    ~~~ FUNCTION TESTS ~~~
    """

    def test_get_params_from_reg(self, capsys, a):
        with capsys.disabled():
            from qcweather import x

    def test_do_multi_lin_reg(self, capsys, a):
        with capsys.disabled():
            from qcweather import y

    def test_write_csv_table(self, a):
        with capsys.disabled():
            from qcweather import z