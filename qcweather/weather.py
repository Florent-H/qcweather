import csv
import os
import re
import subprocess

import pandas as pd
import numpy as np
import matplotlib as plt
import calendar
import datetime
from pathlib import Path
import epw
from collections import OrderedDict
import xlrd

class Weather(object):
    def __init__(self, file_path):
        self.file_path = file_path
        self.meteo_vars = None
        self.wmo = None
        self.folder = None
        self.historical = None

    @classmethod
    def get_weather(cls, file_path):

        # instantiate Weather object
        weather = cls(file_path)

        # if the weather file is an energy plus weather file
        if file_path.suffix == ".epw":
            # open energy plus weather file
            with open(file_path, "r") as fp:
                reader = csv.reader(fp)
                # iterate through each line of the energyplus file
                for i, line in enumerate(reader):

                    # for the first line of code, record the World Meteorological Organization station number and create
                    # the folder for the raw historical weather data
                    if re.match(r"(?i)location", line[0]):
                        weather.wmo = line[5]
                        weather.folder = Path(f"datafiles/raw_tbl/{weather.wmo}")
                        os.mkdir(weather.folder)

                    # if the line starts with 4 digits (the year), then we have passed through the header and can start
                    # recording the values of the meteorological variables to the class's meteo_vars attribute
                    elif re.match(r"\d{4}", line[0]):
                        if not weather.meteo_vars:
                            # get number of hours in the year as it could be a leap year (i.e. 8760 vs 8784 hours)
                            hours = len(list(reader)) + 1
                            first_date = datetime.datetime(
                                int(line[0]), int(line[1]), int(line[2]), int(line[3])
                            )
                            # create dataframe template for meteorological variables for the energyplus weather file
                            weather.meteo_vars = pd.DataFrame(
                                OrderedDict(
                                    (
                                        ("uncertainty", np.empty(hours)),
                                        ("flags", np.empty(hours)),
                                        ("dry_bulb", np.empty(hours)),
                                        ("dew_point", np.empty(hours)),
                                        ("rel_hum", np.empty(hours)),
                                        ("pressure", np.empty(hours)),
                                        ("ext_hor_rad", np.empty(hours)),
                                        ("ext_dir_norm", np.empty(hours)),
                                        ("hor_inf_rad_int", np.empty(hours)),
                                        ("glob_hor_rad", np.empty(hours)),
                                        ("dir_norm_rad", np.empty(hours)),
                                        ("diff_hor_rad", np.empty(hours)),
                                        ("glob_hor_ill", np.empty(hours)),
                                        ("dir_norm_ill", np.empty(hours)),
                                        ("diff_hor_ill", np.empty(hours)),
                                        ("zen_lum", np.empty(hours)),
                                        ("wind_dir", np.empty(hours)),
                                        ("wind_speed", np.empty(hours)),
                                        ("tot_sky_cov", np.empty(hours)),
                                        ("op_sky_cov", np.empty(hours)),
                                        ("visibility", np.empty(hours)),
                                        ("ceil_height", np.empty(hours)),
                                        ("pres_wth_obs", np.empty(hours)),
                                        ("pres_wth_cod", np.empty(hours)),
                                        ("prec_water", np.empty(hours)),
                                        ("aero_opt_dep", np.empty(hours)),
                                        ("snow_depth", np.empty(hours)),
                                        ("days_lst_snow", np.empty(hours)),
                                        ("albedo", np.empty(hours)),
                                        ("liq_prec_dep", np.empty(hours)),
                                        ("liq_prec_qty", np.empty(hours)),
                                        ("hourly_flags", np.empty(hours)),
                                        ("daily_flags", np.empty(hours)),
                                        ("monthly_flags", np.empty(hours)),
                                    )
                                ),
                                index=pd.date_range(
                                    first_date, periods=hours, freq="H"
                                ),
                            )
                            # keep track of the hour number
                            hour = 0

                        # if the meteorological variable dataframe has been created
                        if isinstance(weather.meteo_vars, pd.DataFrame):
                            # write the hour's data to the corresponding columns
                            for j, (col, arr) in enumerate(weather.meteo_vars.iteritems()):
                                if col not in ("hourly_flags", "daily_flags", "monthly_flags"):
                                    try:
                                        # missing value
                                        if line[j + 4] == "9999":
                                            weather.meteo_vars["hourly_flags"][hour] = "missing value"
                                            weather.meteo_vars[j][hour] = None

                                        weather.meteo_vars[j][hour] = float(line[j + 4])
                                    except:  # todo: check error
                                        weather.meteo_vars[j][hour] = str(line[j + 4])
                            hour += 1

        # gather historical data for the meteorological variables
        # run tblexpand to gather raw historical weather data
        run_string = (
            str(Path("E:/tblxpand.exe"))
            + f" {weather.folder} "
            + str(Path(f"E:/data/{weather.wmo}.wdv"))
            + " ALL 99 SI"
        )
        subprocess.run(run_string, cwd=Path.cwd(), capture_output=True)

        # get monthly CDFs
        db_cdf = get_cdf(weather, "dry_bulb")
        dp_cdf = get_cdf(weather, "dew_point")
        wd_cdf = get_cdf(weather, "wind_dir")
        ws_cdf = get_cdf(weather, "wind_speed")

        # get irradiation data


        weather.historical = pd.DataFrame(
            OrderedDict(
                (
                    ("dry_bulb_cdf", db_cdf),
                    ("dew_point_cdf", dp_cdf),
                    ("dew_point_cdf", wd_cdf),
                    ("rel_hum_cdf", ws_cdf),
                    ("pressure", np.empty(hours)),
                    ("ext_hor_rad", np.empty(hours)),
                    ("ext_dir_norm", np.empty(hours)),
                    ("hor_inf_rad_int", np.empty(hours)),
                    ("glob_hor_rad", np.empty(hours)),
                    ("dir_norm_rad", np.empty(hours)),
                    ("diff_hor_rad", np.empty(hours)),
                    ("glob_hor_ill", np.empty(hours)),
                    ("dir_norm_ill", np.empty(hours)),
                    ("diff_hor_ill", np.empty(hours)),
                    ("zen_lum", np.empty(hours)),
                    ("wind_dir", np.empty(hours)),
                    ("wind_speed", np.empty(hours)),
                    ("tot_sky_cov", np.empty(hours)),
                    ("op_sky_cov", np.empty(hours)),
                    ("visibility", np.empty(hours)),
                    ("ceil_height", np.empty(hours)),
                    ("pres_wth_obs", np.empty(hours)),
                    ("pres_wth_cod", np.empty(hours)),
                    ("prec_water", np.empty(hours)),
                    ("aero_opt_dep", np.empty(hours)),
                    ("snow_depth", np.empty(hours)),
                    ("days_lst_snow", np.empty(hours)),
                    ("albedo", np.empty(hours)),
                    ("liq_prec_dep", np.empty(hours)),
                    ("liq_prec_qty", np.empty(hours)),
                    ("hourly_flags", np.empty(hours)),
                    ("daily_flags", np.empty(hours)),
                    ("monthly_flags", np.empty(hours)),
                )
            ),
            index=pd.date_range(
                first_date, periods=hours, freq="H"
            ),
        )



        return weather

    def database_qc(self):
        pass

    def magnitude_qc(self):
        # Hourly checks

        # Checks that Global Horizontal Radiation falls in the prescribed range
        for time, rad in self.meteo_vars["glob_hor_rad"]:
            if


    def step_qc(self):
        pass

    def time_qc(self):
        pass

def get_cdf(weather, meteo_var):

    if meteo_var == "dry_bulb":
        meteo_id = "DB"

    elif meteo_var == "dew_point":
        meteo_id = "DP"

    elif meteo_var == "wind_dir":
        meteo_id = "WD"

    elif meteo_var == "wind_speed":
        meteo_id = "WS"

    else:
        raise ValueError("meteo_var must be either 'dry_bulb', 'dew_point', 'wind_dir', or 'wind_speed'")

    cdf = [{"bin": [], "freq": [], "cdf": []} for i in range(12)]

    for i in range(12):
        if i < 9:
            month = f"0{i + 1}"
        else:
            month = f"{i + 1}"
        month_path = weather.folder / f"{weather.wmo}_{meteo_id}_{month}.txt"
        with open(month_path, "r") as fp:
            reader = csv.reader(fp, delimeter=' ')
            # skip the five-line header
            for n in range(5):
                next(reader)

            for line in reader:
                cdf[i]["bin"].append(float(line[0]))
                cdf[i]["freq"].append(float(line[1]))
                cdf[i]["cdf"].append(float(line[2]))

    return cdf




