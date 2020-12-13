from qcweather.weather import Weather
from pathlib import Path

def run_quality_control(file_path):
    # create pandas dataframe with the meteorological variables from weather file
    weather = Weather.get_weather(file_path)

    # # run quality assurance on database
    # weather.database_qc()

    # run quality assurance on the hourly magnitudes of the meteorological variables
    weather.hour_qc()

    # run quality assurance on the hourly steps of the meteorological variables
    weather.step_qc()

    # run quality assurance on the hourly steps of the meteorological variables
    weather.day_qc()

    # run quality assurance on the monthly magnitudes of the meteorological variables
    weather.month_qc()

    # weather.get_graph("hourly", "dew_point", month=2)
    weather.get_graph("daily", "dry_bulb", month=5)
    weather.get_graph("monthly", "glob_hor_rad")

    csv_path = Path(f"datafiles/results/{weather.file_path.stem}_results.csv")
    weather.meteo_vars.to_csv(csv_path)


# # run quality assurance on the difference in magnitude between adjacent data points for each meteorological
    # # variables
    # weather.step_qc()
    #
    # # run quality assurance on the time series for each meteorological variable
    # weather.time_qc()









def get_wet_bulb(ep_path):
    pass


def get_dew_point(ep_path):
    pass


def get_rel_hum(ep_path):
    pass


def get_dir_irrad(ep_path):
    pass


def get_diff_irrad(ep_path):
    pass


def get_glob_irrad(ep_path):
    pass


def get_wind_speed(ep_path):
    pass


def get_wind_dir(ep_path):
    pass
