import pytest
import os
from pathlib import Path
from qcweather import Weather
from qcweather import get_weather

"""
~~~ PYTEST FIXTURES ~~~
"""


@pytest.fixture()
def weather_path():
    from pathlib import Path

    # CAN-QC - Montreal YUL 716270 - ISD 2015
    yield "testfiles\\epw_files\\tampered.epw"
    # yield "datafiles\\epw_files\\CAN_QC_Montreal-McTavish.716120_CWEC2016.epw"
    # yield "datafiles\\epw_files\\CAN-QC - Montreal YUL 716270 - ISD 2015.epw"

    # yield Path("datafiles/epw_files/tampered.epw")
    # yield Path("datafiles/epw_files/CAN_QC_Montreal-McTavish.716120_CWEC2016.epw")
    # yield Path("datafiles/epw_files/CAN-QC - Montreal YUL 716270 - ISD 2015.epw")


@pytest.fixture()
def drive_letter():
    yield "E"


@pytest.fixture()
def ashrae_path():
    yield "testfiles\\ashrae\\2017DesignConditions_short.xlsx"


@pytest.fixture()
def output_dir():
    yield "testfiles\\figures"


@pytest.fixture()
def weather(weather_path, drive_letter, ashrae_path, output_dir):
    weather = get_weather(weather_path, drive_letter, ashrae_path, output_dir)
    yield weather


@pytest.fixture()
def config(output_dir):
    import shutil

    if os.path.isdir(Path(output_dir) / "raw_tbl"):
        shutil.rmtree(Path(output_dir) / "raw_tbl")
    (Path(output_dir) / "raw_tbl").mkdir()


class TestsQualityControl:
    """A class that encapsulates the different tests for the functions found in
    the quality_control.py module
    """

    """
    ~~~ FUNCTION TESTS ~~~
    """

    def test_get_weather(
        self, weather_path, drive_letter, ashrae_path, output_dir, capsys
    ):

        from qcweather import get_weather

        with capsys.disabled():
            weather = get_weather(weather_path, drive_letter, ashrae_path, output_dir)

    def test_get_meteo_params(
        self, weather_path, drive_letter, ashrae_path, output_dir, capsys
    ):
        # create barebones Weather object
        weather = Weather(weather_path, drive_letter, ashrae_path, output_dir)

        # call the get_meteo_params() method to get the meteorological variables in a pandas Dataframe
        weather.get_meteo_params()

        # determine if the meteorological parameters (longitude, wmo station, etc.) are correctly set
        assert weather.location == (0.7941248096574199, 1.2842132636174277)
        assert weather.time_zone == 5.0
        assert weather.wmo == 716120

        # assert that the pandas Dataframe meteo_vars has X columns
        assert weather.meteo_vars.columns.tolist() == [
            "dry_bulb",
            "dew_point",
            "rel_hum",
            "pressure",
            "ext_hor_rad",
            "ext_dir_norm",
            "hor_inf_rad_int",
            "glob_hor_rad",
            "dir_norm_rad",
            "diff_hor_rad",
            "glob_hor_ill",
            "dir_norm_ill",
            "diff_hor_ill",
            "zen_lum",
            "wind_dir",
            "wind_speed",
            "tot_sky_cov",
            "op_sky_cov",
            "visibility",
            "ceil_height",
            "pres_wth_obs",
            "pres_wth_cod",
            "prec_water",
            "aero_opt_dep",
            "snow_depth",
            "days_lst_snow",
            "albedo",
            "liq_prec_dep",
            "liq_prec_qty",
            "hourly_flags",
            "step_flags",
            "daily_flags",
            "monthly_flags",
        ]
        assert len(weather.meteo_vars.index) == 8760

    def test_get_solar_angles(self, weather_path, drive_letter, ashrae_path, output_dir, capsys):
        # create barebones Weather object
        weather = Weather(weather_path, drive_letter, ashrae_path, output_dir)

        # call get_meteo_params() for the meteo_vars index and location tuple required for get_solar_angles()
        weather.get_meteo_params()

        # call the get_solar_angles() method to get the solar angles as pandas Dataframe
        weather.get_solar_angles()

        # determine if the solar angles are correctly calculated
        assert weather.solar_angles.columns.tolist() == ["declination", "hour_angle", "zenith", "day"]
        assert len(weather.solar_angles.index) == 8760
        assert weather.solar_angles.iloc[0, 0] == -0.4030647383382704
        assert weather.solar_angles.iloc[4380, 0] == 0.40335665337559756
        assert weather.solar_angles.iloc[0, 1] == -2.997640086680232
        assert weather.solar_angles.iloc[4380, 1] == 0.13968758543712329
        assert weather.solar_angles.iloc[0, 2] == 2.733392502374424
        assert weather.solar_angles.iloc[4380, 2] == 0.40693700654613624
        assert weather.solar_angles.iloc[0, 3] is False
        assert weather.solar_angles.iloc[4380, 3] is True

    def test_get_get_historical(self, weather_path, drive_letter, ashrae_path, output_dir, capsys):
        # create barebones Weather object
        weather = Weather(weather_path, drive_letter, ashrae_path, output_dir)

        # call get_historical to get historical probability distribution functions along with monthly averages and
        # standard deviations
        weather.get_historical()





    def test_run_quality_control(
        self, config, capsys, weather_path, drive_letter, ashrae_path, output_dir
    ):

        with capsys.disabled():
            weather = Weather.get_weather(
                weather_path, drive_letter, ashrae_path, output_dir
            )
            weather.run_quality_control()

    def test_get_month_rad_graphs(self, config, capsys):

        with capsys.disabled():
            from qcweather import get_month_rad_graphs

            get_month_rad_graphs(drive_letter, ashrae_path, output_dir)

    def test_get_month_rad_extreme_graphs(
        self, config, capsys, drive_letter, ashrae_path, output_dir
    ):
        with capsys.disabled():
            from qcweather import get_month_rad_extreme_graphs

            get_month_rad_extreme_graphs(drive_letter, ashrae_path, output_dir)

    def test_get_graph(self, config, capsys, weather_path):
        with capsys.disabled():
            weather = Weather.get_weather(
                weather_path, drive_letter, ashrae_path, output_dir
            )
            weather.get_graph("hourly")
