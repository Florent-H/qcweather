import csv
import itertools
import math
import os
import re
import subprocess

import pandas as pd
import numpy as np
import matplotlib as plt
import calendar
import datetime
from pathlib import Path
# import epw
from collections import OrderedDict
import xlrd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
from scipy.stats import norm


class Weather(object):
    def __init__(self, file_path):
        self.file_path = file_path
        self.wmo = None
        self.folder = None
        self.location = None
        self.time_zone = None
        self.meteo_vars = None
        self.solar_angles = None
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
                data = list(reader)

        # for the first line of code, record the World Meteorological Organization station number and create
        # the folder for the raw historical weather data
        for line in data:
            if re.match(r"(?i)location", line[0]):
                weather.wmo = int(line[5])
                weather.folder = Path(f"datafiles/raw_tbl/{weather.wmo}")
                os.mkdir(weather.folder)
                # (latitude, longitude) in degrees
                # using convention that to the west is positive
                weather.location = (
                    math.radians(float(line[6])),
                    math.radians(-float(line[7])),
                )
                weather.time_zone = -float(line[8])

            # if the line starts with 4 digits (the year), then we have passed through the header and can start
            # recording the values of the meteorological variables to the class's meteo_vars attribute
            elif re.match(r"\d{4}", line[0]):

                if weather.meteo_vars is None:
                    # get number of hours in the year as it could be a leap year (i.e. 8760 vs 8784 hours)
                    hours = len(data) - 8
                    first_date = datetime.datetime(
                        int(line[0]), int(line[1]), int(line[2]), int(line[3])
                    )

                    weather.meteo_vars = pd.DataFrame(
                        OrderedDict(
                            (
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
                                ("hourly_flags", [[] for i in range(hours)]),
                                ("step_flags", [[] for i in range(hours)]),
                                ("daily_flags", [[] for i in range(hours)]),
                                ("monthly_flags", [[] for i in range(hours)]),
                            )
                        ),
                        index=pd.date_range(first_date, periods=hours, freq="H"),
                    )
                    # keep track of the hour number
                    hour = 0

                # if the meteorological variable dataframe has been created
                if isinstance(weather.meteo_vars, pd.DataFrame):
                    # write the hour's data to the corresponding columns
                    for j, (col, arr) in enumerate(weather.meteo_vars.iteritems()):
                        if col not in (
                            "hourly_flags",
                            "step_flags",
                            "daily_flags",
                            "monthly_flags",
                        ):
                            # missing value
                            if line[j + 6] == "9999":
                                weather.meteo_vars["hourly_flags"][hour].append(
                                    "missing value"
                                )
                                weather.meteo_vars[col][hour] = None

                            weather.meteo_vars[col][hour] = float(line[j + 6])
                    hour += 1

        # get solar angles
        weather.solar_angles = get_solar_angles(weather)

        # gather historical data for the meteorological variables

        # run tblexpand to gather raw historical weather data
        run_string = (
            str(Path("E:/tblxpand/tblxpand.exe"))
            + f" {weather.folder} "
            + str(Path(f"E:/data/{weather.wmo}.wdv"))
            + " ALL 99 SI"
        )
        subprocess.run(run_string, cwd=Path.cwd(), capture_output=True)

        # get monthly CDFs
        db_pdf = get_pdf(weather, "dry_bulb")
        dp_pdf = get_pdf(weather, "dew_point")
        rel_pdf = get_pdf(weather, "rel_hum")
        wd_pdf = get_pdf(weather, "wind_dir")
        ws_pdf = get_pdf(weather, "wind_speed")
        dbtd_pdf = get_pdf(weather, "dry_bulb_time")

        # get irradiation data
        # todo: generalize workbook location
        w_b = xlrd.open_workbook(Path("datafiles/ashrae/2017DesignConditions_s.xlsx"))

        # read first and only sheet of the excel document
        # todo: generalize selection of sheet index
        sheet_index = 0
        excel_sheet = w_b.sheet_by_index(sheet_index)

        clr_sky_dir_norm = get_month_data(excel_sheet, weather.wmo, "clr_sky_dir_norm")
        clr_sky_diff_hor = get_month_data(excel_sheet, weather.wmo, "clr_sky_diff_hor")
        dir_norm_taub = get_month_data(excel_sheet, weather.wmo, "dir_norm_taub")
        diff_hor_taud = get_month_data(excel_sheet, weather.wmo, "diff_hor_taud")
        all_sky_glob_avg = get_month_data(excel_sheet, weather.wmo, "all_sky_glob_avg")
        all_sky_glob_std = get_month_data(excel_sheet, weather.wmo, "all_sky_glob_std")
        prec_dep_avg = get_month_data(excel_sheet, weather.wmo, "prec_dep_avg")
        prec_dep_max = get_month_data(excel_sheet, weather.wmo, "prec_dep_max")
        prec_dep_min = get_month_data(excel_sheet, weather.wmo, "prec_dep_min")
        prec_dep_std = get_month_data(excel_sheet, weather.wmo, "prec_dep_std")

        weather.historical = OrderedDict(
            (
                ("dry_bulb_pdf", db_pdf),
                ("dew_point_pdf", dp_pdf),
                ("rel_hum_pdf", rel_pdf),
                ("wind_dir_pdf", wd_pdf),
                ("wind_speed_pdf", ws_pdf),
                ("dry_bulb_time_pdf", dbtd_pdf),
                ("clr_sky_dir_norm", clr_sky_dir_norm),
                ("clr_sky_diff_hor", clr_sky_diff_hor),
                ("dir_norm_taub", dir_norm_taub),
                ("diff_hor_taud", diff_hor_taud),
                ("all_sky_glob_avg", all_sky_glob_avg),
                ("all_sky_glob_std", all_sky_glob_std),
                ("prec_dep_avg", prec_dep_avg),
                ("prec_dep_max", prec_dep_max),
                ("prec_dep_min", prec_dep_min),
                ("prec_dep_std", prec_dep_std),
            )
        )

        return weather

    def database_qc(self):
        pass

    def hour_qc(self):
        # check that meteorological variable values do not exceed thresholds

        for time, col in self.meteo_vars.iterrows():

            self.rad_hour_checks(time)
            self.db_hour_checks(time)
            self.dp_hour_check(time)
            self.rh_hour_check(time)
            self.ws_hour_check(time)

    def rad_hour_checks(self, time):

        # -----------------------
        # RADIATION / IRRADIANCE
        # -----------------------
        # (since radiation is in units of Wh/m^2 in a time series with hour time steps, both
        # are equivalent)

        # Solar constant
        G_sc = 1361

        # todo: add smaller than 0.03 * GHItoa (G_o)

        # if it is day and the solar zenith angle is less than 90 degrees
        if (
            self.solar_angles["day"][time]
            and self.solar_angles["zenith"][time] < math.pi / 2
        ):

            # calculate extraterrestrial solar radiation
            n = time.dayofyear - 0.5 + (time.hour - 0.5) / 24
            G_on = G_sc * (1 + 0.0334 * math.cos(math.radians((360 * (n - 3)) / 365)))
            G_o = G_on * math.cos(self.solar_angles["zenith"][time])

            # get month number
            month = time.month - 1

            # ---------------------------------------
            # global horizontal radiation/illuminance
            # ---------------------------------------

            # impossible physics check
            if self.meteo_vars["glob_hor_rad"][time] > min(
                1.2 * G_on,
                1.5 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 100,
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is physically impossible [3 4 7]"
                )

            # rare check
            if (
                self.meteo_vars["glob_hor_rad"][time]
                > 1.2 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 50
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is rare [4 5 7]"
                )

            # if it's not the first hour before sunrise or the final hour before sunset
            if (
                self.meteo_vars["glob_hor_rad"][time] < 0.03 * G_o
                and self.solar_angles["day"][time + datetime.timedelta(hours=-1)]
                and self.solar_angles["day"][time + datetime.timedelta(hours=+1)]
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation is rare [4 5 7]"
                )

            # -----------------------------------
            # direct normal radiation/illuminance
            # -----------------------------------

            # impossible physics check
            if self.meteo_vars["dir_norm_rad"][time] > G_on:
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is physically impossible [4]"
                )

            # rare check
            if (
                self.meteo_vars["dir_norm_rad"][time]
                > 0.95 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 0.2 + 10
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is rare [4]"
                )

            # if it's not the first hour before sunrise or the final hour before sunset
            if (
                self.meteo_vars["dir_norm_rad"][time] == 0
                and self.solar_angles["day"][time + datetime.timedelta(hours=-1)]
                and self.solar_angles["day"][time + datetime.timedelta(hours=+1)]
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation is rare [4]"
                )

            # ----------------------------------------
            # diffuse horizontal radiation/illuminance
            # ----------------------------------------

            # impossible physics check
            if self.meteo_vars["diff_hor_rad"][time] > min(
                0.8 * G_on,
                0.95 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 50,
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is physically impossible [3 4 5 7]"
                )

            # rare check
            if (
                self.meteo_vars["diff_hor_rad"][time]
                > 0.75 * G_on * math.cos(self.solar_angles["zenith"][time]) ** 1.2 + 30
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is rare [3 4]"
                )

            # if it's not the first hour before sunrise
            if (
                self.meteo_vars["diff_hor_rad"][time] < 0.03 * G_o
                and self.solar_angles["day"][time + datetime.timedelta(hours=-1)]
            ):
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation is rare [3 4]"
                )

        # if it is night
        if not self.solar_angles["day"][time]:
            # non-zero radiation during night check
            if self.meteo_vars["glob_hor_rad"][time] > 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "global horizontal radiation during night"
                )

            # non-zero radiation during night check
            if self.meteo_vars["dir_norm_rad"][time] > 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "direct normal radiation during night"
                )

            # non-zero radiation during night check
            if self.meteo_vars["diff_hor_rad"][time] > 0:
                self.meteo_vars["hourly_flags"][time].append(
                    "diffuse horizontal radiation during night"
                )

    def db_hour_checks(self, time):

        # --------------------
        # DRY BULB TEMPERATURE
        # --------------------

        # impossible physics check
        if (
            self.meteo_vars["dry_bulb"][time] < -90
            or self.meteo_vars["dry_bulb"][time] > 60
        ):

            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is physically impossible or exceedingly rare, breaking world records "
                "[10 14]"
            )

        # rare check
        if (
            self.meteo_vars["dry_bulb"][time] < -80
            or self.meteo_vars["dry_bulb"][time] > 50
        ):
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is rare [10 14]"
            )

        # 1st and 99th percentile check
        month = int(time.month) - 1

        perc_1st_idx = min(
            range(len(self.historical["dry_bulb_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dry_bulb_pdf"]["cdf"][month][i] - 0.01),
        )
        perc_1st = self.historical["dry_bulb_pdf"]["bin"][month][perc_1st_idx]

        perc_99th_idx = min(
            range(len(self.historical["dry_bulb_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dry_bulb_pdf"]["cdf"][month][i] - 0.99),
        )
        perc_99th = self.historical["dry_bulb_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["dry_bulb"][time] < perc_1st:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the lower 1st percentile"
            )

        if self.meteo_vars["dry_bulb"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 99th percentile"
            )

    def dp_hour_check(self, time):

        # --------------------
        # DEWPOINT TEMPERATURE
        # --------------------

        # impossible physics check
        # todo: check that it is inferior to dry bulb

        # todo: check that it is superior to the dew_point as calculated by the dry bulb at 0% relative humidity

        month = int(time.month) - 1

        # 1st and 99th percentile check
        perc_1st_idx = min(
            range(len(self.historical["dew_point_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dew_point_pdf"]["cdf"][month][i] - 0.01),
        )
        perc_1st = self.historical["dew_point_pdf"]["bin"][month][perc_1st_idx]

        perc_99th_idx = min(
            range(len(self.historical["dew_point_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["dew_point_pdf"]["cdf"][month][i] - 0.99),
        )
        perc_99th = self.historical["dew_point_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["dew_point"][time] < perc_1st:
            self.meteo_vars["hourly_flags"][time].append(
                "dew point is in the lower 1st percentile"
            )

        if self.meteo_vars["dew_point"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "dew point is in the upper 99th percentile"
            )

    def rh_hour_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        # impossible physics check
        if (
            self.meteo_vars["rel_hum"][time] < 0
            or self.meteo_vars["rel_hum"][time] > 100
        ):

            self.meteo_vars["hourly_flags"][time].append(
                "relative humidity is physically impossible"
            )

        month = int(time.month) - 1

        # 1st and 99th percentile check
        perc_1st_idx = min(
            range(len(self.historical["rel_hum_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["rel_hum_pdf"]["cdf"][month][i] - 0.01),
        )
        perc_1st = self.historical["rel_hum_pdf"]["bin"][month][perc_1st_idx]

        perc_99th_idx = min(
            range(len(self.historical["rel_hum_pdf"]["cdf"][month])),
            key=lambda i: abs(self.historical["rel_hum_pdf"]["cdf"][month][i] - 0.99),
        )
        perc_99th = self.historical["rel_hum_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["rel_hum"][time] < perc_1st:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the lower 1st percentile"
            )

        if self.meteo_vars["rel_hum"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "dry bulb temperature is in the upper 99th percentile"
            )

    def ws_hour_check(self, time):
        # 99th percentile check
        month = int(time.month) - 1

        perc_99th_idx = min(
            range(len(self.historical["wind_speed_pdf"]["cdf"][month])),
            key=lambda i: abs(
                self.historical["wind_speed_pdf"]["cdf"][month][i] - 0.99
            ),
        )
        perc_99th = self.historical["wind_speed_pdf"]["bin"][month][perc_99th_idx]

        if self.meteo_vars["wind_speed"][time] > perc_99th:
            self.meteo_vars["hourly_flags"][time].append(
                "wind speed is in the upper 99th percentile"
            )

    def step_qc(self):

        for time, _ in list(self.meteo_vars.iterrows())[:-1]:
            self.db_step_check(time)
            self.rh_step_check(time)
            self.ws_step_check(time)

    def db_step_check(self, time):

        # --------------------
        # DRY BULB TEMPERATURE
        # --------------------

        idx = self.meteo_vars.index.get_loc(time)
        db_step = (
            self.meteo_vars["dry_bulb"][idx + 1] - self.meteo_vars["dry_bulb"][idx]
        )

        # Exceedingly rare check
        if db_step > 8:
            self.meteo_vars["step_flags"][idx].append(
                "dry bulb temperature step is rare [8 9]"
            )
            self.meteo_vars["step_flags"][idx + 1].append(
                "dry bulb temperature step is rare [8 9]"
            )

    def rh_step_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        idx = self.meteo_vars.index.get_loc(time)
        rh_step = self.meteo_vars["rel_hum"][idx + 1] - self.meteo_vars["rel_hum"][idx]

        # Exceedingly rare check
        if rh_step > 30:
            self.meteo_vars["step_flags"][idx].append(
                "relative humidity temperature step is rare [9]"
            )
            self.meteo_vars["step_flags"][idx + 1].append(
                "relative humidity temperature step is rare [9]"
            )

    def ws_step_check(self, time):

        # -----------------
        # RELATIVE HUMIDITY
        # -----------------

        idx = self.meteo_vars.index.get_loc(time)
        ws_step = (
            self.meteo_vars["wind_speed"][idx + 1] - self.meteo_vars["wind_speed"][idx]
        )

        # Exceedingly rare check
        if ws_step > 15:
            self.meteo_vars["step_flags"][idx].append("wind speed step is rare [8 9]")
            self.meteo_vars["step_flags"][idx + 1].append(
                "wind speed step is rare [8 9]"
            )

    def day_qc(self):
        for time in pd.date_range(
            start=self.meteo_vars.index[0], end=self.meteo_vars.index[-1], freq="D"
        ):
            self.db_day_check(time)

    def db_day_check(self, time):
        # get all indices
        idx = self.meteo_vars.index.get_loc(time)
        day_slice = slice(idx, idx + 24, 1)
        month = int(time.month) - 1

        # start a counter with the number of hours that are in the 1st and 99th percentile
        suspect_flag_count = []

        for i, (hour, _) in enumerate(self.meteo_vars[day_slice].iterrows()):
            # 1st and 99th percentile check
            perc_1st_idx = min(
                range(len(self.historical["dry_bulb_time_pdf"]["cdf"][month][i])),
                key=lambda j: abs(
                    self.historical["dry_bulb_time_pdf"]["cdf"][month][i][j] - 0.01
                ),
            )
            perc_1st = self.historical["dry_bulb_time_pdf"]["bin"][month][i][
                perc_1st_idx
            ]

            perc_99th_idx = min(
                range(len(self.historical["dry_bulb_time_pdf"]["cdf"][month][i])),
                key=lambda j: abs(
                    self.historical["dry_bulb_time_pdf"]["cdf"][month][i][j] - 0.99
                ),
            )
            perc_99th = self.historical["dry_bulb_time_pdf"]["bin"][month][i][
                perc_99th_idx
            ]

            if self.meteo_vars["dry_bulb"][hour] < perc_1st:
                suspect_flag_count.append(hour)

            if self.meteo_vars["dry_bulb"][hour] > perc_99th:
                suspect_flag_count.append(hour)

        # if six or more hours during the day are suspect because their corresponding temperatures fall in either the
        # 1st and 99th percentile, flag as a suspect daily profile
        if len(suspect_flag_count) >= 6:
            for t in suspect_flag_count:
                self.meteo_vars["daily_flags"][t].append(
                    "The daily temperature profile of this day is suspect"
                )

    def month_qc(self):
        for time in pd.date_range(
            start=self.meteo_vars.index[0], periods=12, freq="MS"
        ):
            # self.prec_month_check(time)
            self.rad_month_check(time)

    def prec_month_check(self, time):
        # -------------------
        # WATER PRECIPITATION
        # -------------------

        # get all indices
        idx = self.meteo_vars.index.get_loc(time)
        days_in_month = calendar.monthrange(time.year, time.month)[1]
        month_slice = slice(idx, idx + days_in_month * 24, 1)
        month = int(time.month) - 1

        # max and min check
        if (
            sum(self.meteo_vars["liq_prec_dep"][month_slice])
            < self.historical["prec_dep_min"][month]
        ):
            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is less than the minimum observed"
                )

        if (
            sum(self.meteo_vars["liq_prec_dep"][month_slice])
            > self.historical["prec_dep_max"][month]
        ):

            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is more than the maximum observed"
                )

        # 1st and 99th percentile check
        perc_1st = (
            self.historical["prec_dep_avg"][month]
            - 2.32635 * self.historical["prec_dep_std"][month]
        )
        perc_99th = (
            self.historical["prec_dep_avg"][month]
            + 2.32635 * self.historical["prec_dep_std"][month]
        )

        if sum(self.meteo_vars["liq_prec_dep"][month_slice]) < perc_1st:
            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is in the lower 1st percentile"
                )

        if sum(self.meteo_vars["liq_prec_dep"][month_slice]) > perc_99th:
            for hour in range(month_slice.start, month_slice.stop, month_slice.step):
                self.meteo_vars["monthly_flags"][hour].append(
                    "water precipitation is in the upper 99th percentile"
                )

    def rad_month_check(self, time):
        # -----------------------
        # RADIATION / IRRADIANCE
        # -----------------------
        # (since radiation is in units of Wh/m^2 in a time series with hour time steps, both
        # are equivalent)

        # todo: add smaller than 0.03 * GHItoa (G_o)

        # get all indices
        idx = self.meteo_vars.index.get_loc(time)
        days_in_month = calendar.monthrange(time.year, time.month)[1]
        days_sum = []
        for day in pd.date_range(start=time, periods=days_in_month, freq="D"):
            day_idx = self.meteo_vars.index.get_loc(day)
            day_slice = slice(day_idx, day_idx + 24, 1)
            days_sum.append(sum(self.meteo_vars["glob_hor_rad"][day_slice]))
        month = int(time.month) - 1

        # # calculate extraterrestrial solar radiation
        # G_sc = 1361
        # n = time.dayofyear + 20 + (11.5) / 24
        # G_on = G_sc * (1 + 0.0334 * math.cos(math.radians((360 * (n - 3)) / 365)))
        # G_o = G_on * math.cos(self.solar_angles["zenith"][time])

        # ---------------------------------------
        # global horizontal radiation/illuminance
        # ---------------------------------------

        # 1st and 99th percentile check
        perc_all_1st = (
            self.historical["all_sky_glob_avg"][month]
            - 2.32635 * self.historical["all_sky_glob_std"][month]
        )
        perc_all_99th = (
            self.historical["all_sky_glob_avg"][month]
            + 2.32635 * self.historical["all_sky_glob_std"][month]
        )

        if np.mean(days_sum) / 1000 < perc_all_1st:
            for hour in range(idx, idx + days_in_month * 24):
                self.meteo_vars["monthly_flags"][hour].append(
                    "daily mean global horizontal radiation in the month is in the lower 1st percentile"
                )

        if np.mean(days_sum) / 1000 > perc_all_99th:

            for hour in range(idx, idx + days_in_month * 24):
                self.meteo_vars["monthly_flags"][hour].append(
                    "daily mean global horizontal radiation in the month is in the upper 99th percentile"
                )

        # # -----------------------------------
        # # direct normal radiation/illuminance
        # # -----------------------------------
        #
        # daylight_hours = sum(self.solar_angles["day"][month_slice])
        # day_adj_perc_1st = perc_all_1st * (days_in_month * 24) / daylight_hours
        # day_adj_perc_99th = perc_all_99th * (days_in_month * 24) / daylight_hours
        #
        # # Use Erbs model to split global horizontal average in diffuse and direct
        # k_t_1st = day_adj_perc_1st / G_o
        # Id_over_I_1st = self.get_Id_over_I(k_t_1st)
        #
        # k_t_99th = day_adj_perc_1st / G_o
        # Id_over_I_99th = self.get_Id_over_I(k_t_99th)
        #
        # perc_dir_1st =
        #
        #
        # ab = (
        #     1.454
        #     - 0.406 * self.historical["dir_nor_taub"][month_idx]
        #     - 0.268 * self.historical["diff_hor_taud"][month_idx]
        #     + 0.021
        #     * self.historical["dir_nor_taub"][month_idx]
        #     * self.historical["diff_hor_taud"][month_idx]
        # )
        #
        # ad = (
        #     0.507
        #     + 0.205 * self.historical["dir_nor_taub"][month_idx]
        #     - 0.08 * self.historical["diff_hor_taud"][month_idx]
        #     - 0.19
        #     * self.historical["dir_nor_taub"][month_idx]
        #     * self.historical["diff_hor_taud"][month_idx]
        # )
        #
        # # 1st and 99th percentile check
        # perc_dir_1st = perc_all_1st * x
        # perc_dir_99th = perc_all_99th * x
        #
        # if sum(self.meteo_vars["glob_hor_rad"][month_slice]) < perc_dir_1st:
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )
        #
        # if sum(self.meteo_vars["glob_hor_rad"][month_slice]) > perc_dir_99th:
        #
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )
        #
        # # ----------------------------------------
        # # diffuse horizontal radiation/illuminance
        # # ----------------------------------------
        #
        # # 1st and 99th percentile check
        # # perc_diff_1st = perc_all_1st * x
        # # perc_diff_99th = perc_all_99th * x
        #
        # if sum(self.meteo_vars["glob_hor_rad"][month_slice]) < perc_diff_1st:
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )
        #
        # elif sum(self.meteo_vars["glob_hor_rad"][month_slice]) > perc_diff_99th:
        #
        #     for hour in range(month_slice.start, month_slice.stop, month_slice.step):
        #         self.meteo_vars["monthly_flags"][hour].append(
        #             "monthly aggegrated global horizontal radiation is in the lower 1st percentile"
        #         )

    # def get_Id_over_I(self, k_t):
    #     if k_t <= 0.22:
    #         Id_over_I = 1.0 - 0.09 * k_t
    #     elif k_t > 0.22 and k_t <= 0.80:
    #         Id_over_I = 0.9511 - 0.1604 * k_t + 4.388 * k_t ** 2 - 16.638 * k_t ** 3 + 12.336 * k_t ** 4
    #     elif k_t > 0.8:
    #         Id_over_I = 0.165

    def get_graph(self, check, meteo_var, month=None):

        if check == "hourly":

            string_check = meteo_var.split("_")[0]
            value_list = []
            time_list = []

            # get flagged values for the meteorological variable
            for time, value in self.meteo_vars[meteo_var].iteritems():
                for flag in self.meteo_vars["hourly_flags"][time]:
                    if re.match(fr".*{string_check}", flag) and time.month == month + 1:
                        value_list.append(value)
                        time_list.append(time)

            # create a day slice of the indices to get the daily profile of the meteorological variable
            midnight = datetime.datetime(
                time_list[0].year,
                time_list[0].month,
                time_list[0].day,
            )
            midnight_idx = self.meteo_vars.index.get_loc(midnight)
            day_slice = slice(midnight_idx, midnight_idx + 25, 1)
            daily_profile = self.meteo_vars[meteo_var][day_slice]

            pdf_name = meteo_var + "_pdf"

            # get 1st percentile and 99th percentile values
            perc_1st_idx = min(
                range(len(self.historical[pdf_name]["cdf"][month])),
                key=lambda i: abs(self.historical[pdf_name]["cdf"][month][i] - 0.01),
            )
            perc_1st = self.historical[pdf_name]["bin"][month][perc_1st_idx]

            perc_99th_idx = min(
                range(len(self.historical[pdf_name]["cdf"][month])),
                key=lambda i: abs(self.historical[pdf_name]["cdf"][month][i] - 0.99),
            )
            perc_99th = self.historical[pdf_name]["bin"][month][perc_99th_idx]

            bins = self.historical[pdf_name]["bin"][month]
            pdf = self.historical[pdf_name]["pdf"][month]
            fig = plt.figure(figsize=(3.3, 2.1))
            ax1 = fig.add_subplot(111)
            ax2 = ax1.twinx()
            ax1.grid(alpha=0.2)
            plt.rcParams["font.family"] = "Times New Roman"
            csfont = {"fontname": "Times New Roman"}
            fontsz = {"fontsize": 10}

            ax1.hist(bins, bins=len(bins), weights=pdf, color="gray")
            ax2.plot(daily_profile, range(25), color="blue")
            ax1.set_ylabel("Probability", **csfont, **fontsz)
            ax2.set_ylabel("Hour of the day", color="blue", **csfont, **fontsz)
            ax1.set_xlabel("Dewpoint temperature [°C]", **csfont, **fontsz)

            ax1.set_xticks(np.arange(-35, 20, 5).astype(int))
            ax1.set_xlim([-35, 15])
            ax1.set_yticks(np.arange(0, 0.036, 0.004))
            ax2.set_yticks(np.arange(0, 27, 3).astype(int))
            ax1.set_ylim([0, 0.032])
            ax2.set_ylim([0, 24])

            ax1.set_xticklabels(ax1.get_xticks(), **csfont, **fontsz)
            ax1.set_yticklabels(ax1.get_yticks(), **csfont, **fontsz)
            ax2.set_yticklabels(ax2.get_yticks(), color="blue", **csfont, **fontsz)
            plt.vlines(perc_1st, 0, 24, colors="red")
            plt.vlines(perc_99th, 0, 24, colors="red")

            if self.file_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                plt.title("Typical Meteorological Year", fontweight="bold", **csfont, **fontsz)

            if self.file_path.stem == "tampered":
                plt.title("Tampered Year", fontweight="bold", **csfont, **fontsz)

            plt.tight_layout(pad=0.05)
            plt.savefig(Path(f"datafiles/figures/dew_point_{self.file_path.stem}.svg"))

        elif check == "daily":

            value_list = []
            time_list = []

            # get flagged values for the meteorological variable
            for time, value in self.meteo_vars[meteo_var].iteritems():
                for flag in self.meteo_vars["daily_flags"][time]:
                    if re.match(r".*profile", flag) and time.month == month + 1:
                        value_list.append(value)
                        time_list.append(time)

            # create a day slice of the indices to get the daily profile of the meteorological variable
            midnight = datetime.datetime(
                time_list[0].year,
                time_list[0].month,
                time_list[0].day,
            )
            midnight_idx = self.meteo_vars.index.get_loc(midnight)
            day_slice = slice(midnight_idx, midnight_idx + 25, 1)
            daily_profile = self.meteo_vars[meteo_var][day_slice]

            pdf_name = meteo_var + "_time_pdf"

            plt.figure(figsize=(3.3, 5))
            gs1 = gridspec.GridSpec(24, 1)
            gs1.update(wspace=0, hspace=0)  # set the spacing between axes.
            plt.rcParams["font.family"] = "Times New Roman"
            csfont = {"fontname": "Times New Roman"}
            fontsz = {"fontsize": 10}

            # make a subplot of the histogram for each hour of the month
            for hour in range(24):
                # get 1st percentile and 99th percentile values
                perc_1st_idx = min(
                    range(len(self.historical[pdf_name]["cdf"][month][hour])),
                    key=lambda i: abs(self.historical[pdf_name]["cdf"][month][hour][i] - 0.01),
                )
                perc_1st = self.historical[pdf_name]["bin"][month][hour][perc_1st_idx]

                perc_99th_idx = min(
                    range(len(self.historical[pdf_name]["cdf"][month][hour])),
                    key=lambda i: abs(self.historical[pdf_name]["cdf"][month][hour][i] - 0.99),
                )
                perc_99th = self.historical[pdf_name]["bin"][month][hour][perc_99th_idx]

                bins = self.historical[pdf_name]["bin"][month][hour]
                pdf = self.historical[pdf_name]["pdf"][month][hour]

                plt.axis('on')
                ax1 = plt.subplot(gs1[23 - hour])
                ax2 = ax1.twinx()
                ax1.xaxis.grid(True, alpha=0.2)
                ax1.hist(bins, bins=len(bins), weights=pdf, color="gray")
                ax2.plot(daily_profile[hour:hour + 2], np.arange(hour, hour + 2, 1), color="blue")
                if hour == 12:
                    ax1.set_ylabel("Probability", **csfont, **fontsz)
                    ax2.set_ylabel("Hour of the day", **csfont, **fontsz, color="blue")
                else:
                    for tic in ax1.yaxis.get_major_ticks():
                        tic.tick1line.set_visible(False)
                        tic.tick2line.set_visible(False)
                        tic.label1.set_visible(False)
                        tic.label2.set_visible(False)

                ax1.set_xticks(np.arange(4, 40, 4).astype(int))
                ax1.set_xlim([4, 36])
                ax1.set_yticks([0, 0.04, 0.08])
                ax1.set_yticklabels(["0", "", "0.08"], **csfont, **fontsz)
                ax1.set_ylim([0, 0.08])
                ml = MultipleLocator(2)
                ax1.xaxis.set_minor_locator(ml)

                ax2.set_ylim([hour, hour + 1])
                ax2.set_yticklabels([f"{int(y)}" for y in ax2.get_yticks()], color="blue", **csfont, **fontsz)

                if hour == 0:
                    ax1.set_xlabel("Dry bulb temperature [°C]", **csfont, **fontsz)
                    ax1.set_xticklabels(ax1.get_xticks(), **csfont, **fontsz)
                    ax2.set_yticks([int(hour)])
                else:
                    for tic in ax1.xaxis.get_major_ticks():
                        tic.tick1line.set_visible(False)
                        tic.tick2line.set_visible(False)
                        tic.label1.set_visible(False)
                        tic.label2.set_visible(False)
                    for tic in ax1.xaxis.get_minor_ticks():
                        tic.tick1line.set_visible(False)
                        tic.tick2line.set_visible(False)
                        # tic.label1.set_visible(False)
                        # tic.label2.set_visible(False)

                    ax2.set_yticks([int(hour), int(hour + 1)])

                ax2.vlines(perc_1st, hour, hour + 1, colors="red")
                ax2.vlines(perc_99th, hour, hour + 1, colors="red")

            if self.file_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                plt.title("Typical Meteorological Year", fontweight="bold", **csfont, **fontsz)

            if self.file_path.stem == "tampered":
                plt.title("Tampered Year", fontweight="bold", **csfont, **fontsz)

            plt.tight_layout(pad=0.05)
            plt.savefig(Path(f"datafiles/figures/dry_bulb_profile_{self.file_path.stem}.svg"))

        elif check == "monthly":

            plt.figure(figsize=(3.3, 4))
            gs1 = gridspec.GridSpec(12, 1)
            gs1.update(wspace=0, hspace=0)  # set the spacing between axes.
            plt.rcParams["font.family"] = "Times New Roman"
            plt.rcParams['mathtext.fontset'] = 'cm'

            csfont = {"fontname": "Times New Roman"}
            fontsz = {"fontsize": 10}

            past_days = 364

            for month in range(11, -1, -1):
                mu = self.historical["all_sky_glob_avg"][month]
                sigma = self.historical["all_sky_glob_std"][month]

                perc_1st = mu - 2.32635 * sigma
                perc_99th = mu + 2.32635 * sigma

                z1 = mu - 3 * sigma
                z2 = mu + 3 * sigma

                x = np.arange(z1, z2, 0.001)
                y = norm.pdf(x, mu, sigma)

                # get all indices
                days_in_month = calendar.monthrange(self.meteo_vars.index[0].year, month + 1)[1]
                days_sum = []
                for day in range(past_days, past_days - days_in_month, -1):
                    day_slice = slice(day*24, (day + 1) * 24, 1)
                    days_sum.append(sum(self.meteo_vars["glob_hor_rad"][day_slice]))
                past_days -= days_in_month

                month_average = np.mean(days_sum) / 1000

                plt.axis('on')
                ax1 = plt.subplot(gs1[month])
                ax2 = ax1.twinx()
                ax1.xaxis.grid(True, alpha=0.2, which="both")
                ax1.barh(0, month_average, align="center", color="blue")
                ax2.plot(x, y, color="gray")

                if self.file_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                    ax1.set_xticks(range(8))
                    ax1.set_xlim([0, 7])

                if self.file_path.stem == "CAN-QC - Montreal YUL 716270 - ISD 2015":
                    ax1.set_xticks(range(9))
                    ax1.set_xlim([0, 8])
                ml = MultipleLocator(0.5)
                ax1.xaxis.set_minor_locator(ml)

                ax1.set_yticks([-1, 0, 1])
                ax1.set_yticklabels(['', calendar.month_name[month + 1], ''], color="blue", **csfont, **fontsz)
                for i, tic in enumerate(ax1.yaxis.get_major_ticks()):
                    if i != 1:
                        tic.label1.set_visible(False)
                        tic.label2.set_visible(False)
                    tic.tick1line.set_visible(False)
                    tic.tick2line.set_visible(False)

                if month == 11:
                    ax1.set_xlabel(r"Monthly average" "\n" r"global horizontal irradiation [kWh/$\rm m^2$]", **csfont, **fontsz)
                    ax1.set_xticklabels(ax1.get_xticks(), **csfont, **fontsz)

                else:
                    for tic in ax1.xaxis.get_major_ticks():
                        tic.tick1line.set_visible(False)
                        tic.tick2line.set_visible(False)
                        tic.label1.set_visible(False)
                        tic.label2.set_visible(False)
                    for tic in ax1.xaxis.get_minor_ticks():
                        tic.tick1line.set_visible(False)
                        tic.tick2line.set_visible(False)
                        tic.label1.set_visible(False)
                        tic.label2.set_visible(False)

                ax2.get_yaxis().set_visible(False)

                ax2.vlines(perc_1st, 0, max(y), colors="red")
                ax2.vlines(perc_99th, 0, max(y), colors="red")

            if self.file_path.stem == "CAN_QC_Montreal-McTavish.716120_CWEC2016":
                plt.title("Typical Meteorological Year", fontweight="bold", **csfont, **fontsz)

            if self.file_path.stem == "CAN-QC - Montreal YUL 716270 - ISD 2015":
                plt.title("2015", fontweight="bold", **csfont, **fontsz)

            plt.tight_layout(pad=0.05)
            plt.savefig(Path(f"datafiles/figures/glob_rad_months_{self.file_path.stem}.svg"))


def get_pdf(weather, meteo_var):

    if meteo_var == "dry_bulb":
        meteo_id = "DB"

    elif meteo_var == "dew_point":
        meteo_id = "DP"

    elif meteo_var == "rel_hum":
        meteo_id = "DBDP"

    elif meteo_var == "wind_dir":
        meteo_id = "WD"

    elif meteo_var == "wind_speed":
        meteo_id = "WS"

    elif meteo_var == "dry_bulb_time":
        meteo_id = "DBTD"
    else:
        raise ValueError(
            "meteo_var must be either 'dry_bulb', 'dew_point', 'rel_hum', 'wind_dir', 'wind_speed', or 'dry_bulb_time'"
        )

    if meteo_id == "DBTD":
        pdf = {
            "bin": [[[] for i in range(24)] for i in range(12)],
            "freq": [[[] for i in range(24)] for i in range(12)],
            "cdf": [[[] for i in range(24)] for i in range(12)],
            "pdf": [[[] for i in range(24)] for i in range(12)],
        }

    else:
        pdf = {
            "bin": [[] for i in range(12)],
            "freq": [[] for i in range(12)],
            "cdf": [[] for i in range(12)],
            "pdf": [[] for i in range(12)],
        }

    for m in range(12):
        if m < 9:
            month = f"0{m + 1}"
        else:
            month = f"{m + 1}"
        month_path = weather.folder / f"{weather.wmo}_{meteo_id}_{month}.txt"

        with open(month_path, "r") as fp:
            # skip the five-line or four-line header
            if meteo_id in ("DBDP", "DBTD"):
                skip = 4
            else:
                skip = 5

            for n in range(skip):
                next(fp)

            for j, line in enumerate(fp):
                if meteo_id == "DBDP":
                    if j == 0:
                        dew_point = [float(s) for s in line.split()[1:]]
                        dry_bulb = []
                        freq = []

                    else:
                        dry_bulb.append(float(line.split()[0]))
                        freq.append([float(s) for s in line.split()[1:]])

                elif meteo_id == "DBTD":
                    if j == 0:
                        time_day = [float(s) for s in line.split()[1:]]
                        dry_bulb_day = []
                        freq_day = [[] for s in line.split()[1:]]

                    else:
                        dry_bulb_day.append(float(line.split()[0]))
                        for s, _ in enumerate(freq_day):
                            freq_day[s].append(float(line.split()[s + 1]))
                else:
                    pdf["bin"][m].append(float(line.split()[0]))
                    pdf["freq"][m].append(float(line.split()[1]))
                    pdf["cdf"][m].append(float(line.split()[2]))

        # todo: remove
        if meteo_id == "DBDP":
            bin_unsorted = []
            freq_unsorted = []
            run_freq_sum = [0]

            for x, db in enumerate(dry_bulb):
                for y, dp in enumerate(dew_point):
                    bin_unsorted.append(
                        100
                        * (17.625 * dp / (243.04 + dp))
                        / (17.625 * (db + 0.00000000000000000000001) / (243.04 + db))
                    )
                    freq_unsorted.append(freq[x][y])
                    run_freq_sum.append(run_freq_sum[-1] + freq[x][y])

            # sort bins and corresponding frequencies
            pdf["bin"][m] = sorted(bin_unsorted)
            pdf["freq"][m] = [f for _, f in sorted(zip(pdf["bin"][m], freq_unsorted))]
            pdf["cdf"][m] = [
                c / run_freq_sum[-1]
                for _, c in sorted(zip(pdf["bin"][m], run_freq_sum))
            ]

        if meteo_id is not "DBTD":
            # calculate probability distribution function
            total_obs = sum(pdf["freq"][m])
            pdf["pdf"][m] = [f / total_obs for f in pdf["freq"][m]]

        else:
            for h, _ in enumerate(time_day):
                run_freq_sum = 0
                for x, _ in enumerate(dry_bulb_day):
                    pdf["bin"][m][h].append(dry_bulb_day[x])
                    pdf["freq"][m][h].append(freq_day[h][x])
                    pdf["pdf"][m][h].append(freq_day[h][x] / sum(freq_day[h]))
                    run_freq_sum += freq_day[h][x]
                    pdf["cdf"][m][h].append(run_freq_sum / sum(freq_day[h]))

    return pdf


def get_month_data(excel_sheet, wmo, var_type):

    # find the corresponding row of data for the weather station
    wmo_col = excel_sheet.col_values(4)[4:]
    row = wmo_col.index(str(wmo))

    # read the values from the corresponding columns in the row of the weather station
    if var_type == "clr_sky_dir_norm":
        cells = excel_sheet.row_slice(row, 533, 545)

    elif var_type == "clr_sky_diff_hor":
        cells = excel_sheet.row_slice(row, 545, 557)

    elif var_type == "dir_norm_taub":
        cells = excel_sheet.row_slice(row, 509, 521)

    elif var_type == "diff_hor_taud":
        cells = excel_sheet.row_slice(row, 521, 533)

    elif var_type == "all_sky_glob_avg":
        cells = excel_sheet.row_slice(row, 557, 569)

    elif var_type == "all_sky_glob_std":
        cells = excel_sheet.row_slice(row, 569, 581)

    elif var_type == "prec_dep_avg":
        cells = excel_sheet.row_slice(row, 206, 218)

    elif var_type == "prec_dep_max":
        cells = excel_sheet.row_slice(row, 219, 231)

    elif var_type == "prec_dep_min":
        cells = excel_sheet.row_slice(row, 232, 244)

    elif var_type == "prec_dep_std":
        cells = excel_sheet.row_slice(row, 245, 257)

    else:
        raise ValueError(
            "rad_type must be either 'clr_sky_dir_norm', 'clr_sky_diff_hor', 'dir_norm_taub', 'diff_hor_taud', "
            "'all_sky_glob_avg', 'all_sky_glob_std', 'prec_dep_avg', 'prec_dep_max', 'prec_dep_min', or "
            "'prec_dep_std'"
        )

    # turn N/A values into None and other numerical strings into floats
    month_data = [None if c.value == "N/A" else float(c.value) for c in cells]

    return month_data


def get_solar_angles(weather):

    # get hours in meteorological file
    hours = len(weather.meteo_vars["dry_bulb"])
    latitude = weather.location[0]
    longitude = weather.location[1]

    solar_angles = pd.DataFrame(
        OrderedDict(
            (
                ("declination", np.empty(hours)),
                ("hour_angle", np.empty(hours)),
                ("zenith", np.empty(hours)),
                ("altitude", np.empty(hours)),
                ("azimuth", np.empty(hours)),
                ("day", np.empty(hours)),
            )
        ),
        index=weather.meteo_vars.index,
    )

    for time, _ in weather.meteo_vars.iterrows():
        # fractional day of the year
        n = time.dayofyear - 0.5 + (time.hour - 0.5) / 24
        # equation of time
        B = math.radians((n - 1) * 360 / 365)
        E = 2.2918 * (
            0.0075
            + 0.1868 * math.cos(B)
            - 3.2077 * math.sin(B)
            - 1.4615 * math.cos(2 * B)
            - 4.089 * math.sin(2 * B)
        )

        # solar declination angle
        declin = (
            0.006918
            - 0.399912 * math.cos(B)
            + 0.070257 * math.sin(B)
            - 0.006758 * math.cos(2 * B)
            + 0.000907 * math.sin(2 * B)
            - 0.002697 * math.cos(3 * B)
            + 0.00148 * math.sin(3 * B)
        )
        solar_angles["declination"][time] = declin

        # solar time
        t_s = (
            time.hour
            - 0.5
            + (4 * (weather.time_zone * 15 - math.degrees(longitude)) + E) / 60
        ) % 24
        # solar hour angle
        hour_angle = math.radians((t_s - 12) * 15)
        solar_angles["hour_angle"][time] = hour_angle

        # calculate sunrise and sunset angles
        sunset = math.acos(-math.sin(latitude) * math.tan(declin))
        sunrise = -sunset
        if hour_angle < sunrise:
            solar_angles["day"][time] = False

        elif hour_angle > sunset:
            solar_angles["day"][time] = False

        else:
            solar_angles["day"][time] = True

        # solar zenith angle
        zenith = math.acos(
            math.cos(latitude) * math.cos(declin) * math.cos(hour_angle)
            + math.sin(latitude) * math.sin(declin)
        )

        solar_angles["zenith"][time] = zenith

        altitude = math.asin(math.cos(zenith))
        solar_angles["altitude"][time] = altitude

        # todo: find equivalent of sign() in python (not math.copysign)

        azimuth = math.copysign(
            abs(
                math.acos(
                    (math.cos(zenith) * math.sin(latitude) - math.sin(declin))
                    / (math.sin(zenith) * math.cos(latitude)),
                )
            ),
            hour_angle,
        )
        solar_angles["azimuth"][time] = azimuth

    return solar_angles
