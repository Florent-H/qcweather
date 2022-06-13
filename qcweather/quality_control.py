from qcweather.weather import Weather
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import gridspec
import calendar
import numpy as np
from matplotlib.ticker import MultipleLocator
import math


def run_quality_control(weather):
    # run quality assurance on the hourly magnitudes of the meteorological variables
    weather.hour_qc()

    # run quality assurance on the hourly steps of the meteorological variables
    weather.step_qc()

    # run quality assurance on the hourly steps of the meteorological variables
    weather.day_qc()

    # run quality assurance on the monthly magnitudes of the meteorological variables
    weather.month_qc()

    # get graphs for different time spans and meteorological variables
    weather.get_graph("hourly", "dew_point", month=2)
    weather.get_graph("daily", "dry_bulb", month=5)
    weather.get_graph("monthly", "glob_hor_rad")

    # write results csv file
    csv_path = Path(weather.output_dir) / f"{weather.weather_path.stem}_results.csv"
    weather.meteo_vars.to_csv(csv_path)


def get_month_rad_extreme_graphs(drive_letter, ashrae_path, output_dir):
    rad = {
        "2015": Weather.get_weather(
            Path("datafiles/epw_files/CAN-QC - Montreal YUL 716270 - ISD 2015.epw"),
            drive_letter,
            ashrae_path,
            output_dir,
        ),
        "tmy": Weather.get_weather(
            Path("datafiles/epw_files/CAN_QC_Montreal-McTavish.716120_CWEC2016.epw"),
            drive_letter,
            ashrae_path,
            output_dir,
        ),
        "extreme": Weather(
            Path("datafiles/epw_files/CAN-QC - Montreal YUL 716270 - ISD 2015.epw"),
            drive_letter,
            ashrae_path,
            output_dir,
        ),
    }

    rad["extreme"].meteo_vars = rad["2015"].meteo_vars.copy(deep=True)
    rad["extreme"].solar_angles = rad["2015"].solar_angles

    G_sc = 1361

    for time, _ in rad["extreme"].meteo_vars.iterrows():
        if (
            rad["extreme"].solar_angles["day"][time]
            and rad["extreme"].solar_angles["zenith"][time] < math.pi / 2
        ):

            n = time.dayofyear - 0.5 + (time.hour - 0.5) / 24
            G_on = G_sc * (1 + 0.0334 * math.cos(math.radians((360 * (n - 3)) / 365)))

            rad["extreme"].meteo_vars.loc[time, "glob_hor_rad"] = (
                1.2
                * G_on
                * math.cos(rad["extreme"].solar_angles["zenith"][time]) ** 1.2
                + 50
            )

    csfont = {"fontname": "Arial"}
    fontsz = {"fontsize": 10}
    color = {
        "glob_hor_rad": {"2015": "orange", "tmy": "red", "extreme": "purple"},
        "diff_hor_rad": {"2015": "green", "tmy": "blue"},
        "dir_norm_rad": {"2015": "orange", "tmy": "red"},
    }

    avg_hours = {
        year: rad[year]
        .meteo_vars.groupby(
            [rad[year].meteo_vars.index.month, rad[year].meteo_vars.index.hour]
        )
        .mean()
        for year in rad.keys()
    }
    for key in avg_hours.keys():
        avg_hours[key].index.names = ["month", "hour"]

    for met_vars in ["glob_hor_rad", ["dir_norm_rad", "diff_hor_rad"]]:
        # Monthly irradiation
        # Direct normal irradiation
        fig = plt.figure(figsize=(7, 4), tight_layout=True)
        gs1 = fig.add_gridspec(3, 4, left=0.12, bottom=0.12, wspace=0.18, hspace=0.35)
        if met_vars == "glob_hor_rad":
            ylabel = "Global horizontal irradiation"
        else:
            ylabel = "Solar irradiation"
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams["mathtext.fontset"] = "cm"
        fig.supxlabel("Hour of the day [h]", **csfont, **fontsz)
        fig.supylabel(rf"{ylabel} [kWh/$\rm m^2$]", **csfont, **fontsz)

        for month in range(12):
            ax = plt.subplot(gs1[month])
            ax.grid(alpha=0.2)

            if met_vars == "glob_hor_rad":
                for year in rad.keys():
                    if year == "2015":
                        linestyle = "solid"
                        label = "2015"
                    elif year == "tmy":
                        linestyle = "dashed"
                        label = "TMY"
                    else:
                        linestyle = "dashdot"
                        label = "Extreme"

                    ax.plot(
                        range(24),
                        avg_hours[year].loc[month + 1, met_vars],
                        color=color[met_vars][year],
                        linestyle=linestyle,
                        label=label,
                    )

            else:
                for year in ["2015", "tmy"]:
                    for met_v in met_vars:
                        # get all indices
                        if year == "2015":
                            linestyle = "-"
                        else:
                            linestyle = "--"
                        if met_v == "glob_hor_rad":
                            if year == "tmy":
                                label = "typical meteorological year"
                            else:
                                label = year
                        elif met_v == "diff_hor_rad":
                            if year == "tmy":
                                label = "TMY: Diffuse horizontal irradiation"
                            else:
                                label = "2015: Diffuse horizontal irradiation"
                        else:
                            if year == "tmy":
                                label = "TMY: Direct normal irradiation"
                            else:
                                label = "2015: Direct normal irradiation"

                        ax.plot(
                            range(24),
                            avg_hours[year].loc[month + 1, met_v],
                            color=color[met_v][year],
                            linestyle=linestyle,
                            label=label,
                        )

            ax.set_xticks(np.arange(0, 30, 6).astype(int))
            ax.set_xlim([0, 24])

            if met_vars == "glob_hor_rad":
                ax.set_yticks(np.arange(0, 1750, 250).astype(int))
                ax.set_ylim([0, 1500])
                ml_y = MultipleLocator(250)

            else:
                ax.set_yticks(np.arange(0, 1000, 200).astype(int))
                ax.set_ylim([0, 800])
                ml_y = MultipleLocator(100)

            ax.yaxis.set_minor_locator(ml_y)
            ml_x = MultipleLocator(3)
            ax.xaxis.set_minor_locator(ml_x)
            ax.set_title(calendar.month_name[month + 1], **csfont, **fontsz)
            ax.set_xticklabels(ax.get_xticks(), **csfont, **fontsz)
            ax.set_yticklabels(ax.get_yticks(), **csfont, **fontsz)

            if month not in (8, 9, 10, 11):
                for tic in ax.xaxis.get_major_ticks():
                    tic.label1.set_visible(False)
                    tic.label2.set_visible(False)

            if month % 4 != 0:
                for tic in ax.yaxis.get_major_ticks():
                    tic.label1.set_visible(False)
                    tic.label2.set_visible(False)

        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)

        handles, labels = ax.get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            bbox_to_anchor=(0.55, 1.08),
            loc="lower center",
            labelspacing=0.1,
        )
        plt.savefig(
            f"datafiles/figures/monthly_irradiation_{met_vars}.svg",
            bbox_inches="tight",
            pad_inches=0,
        )


def get_month_rad_graphs(drive_letter, ashrae_path, output_dir):
    rad = {
        "2015": Weather.get_weather(
            Path("datafiles/epw_files/CAN-QC - Montreal YUL 716270 - ISD 2015.epw"),
            drive_letter,
            ashrae_path,
            output_dir,
        ),
        "tmy": Weather.get_weather(
            Path("datafiles/epw_files/CAN_QC_Montreal-McTavish.716120_CWEC2016.epw"),
            drive_letter,
            ashrae_path,
            output_dir,
        ),
    }
    csfont = {"fontname": "Times New Roman"}
    fontsz = {"fontsize": 10}
    color = {
        "glob_hor_rad": {"2015": "orange", "tmy": "red"},
        "diff_hor_rad": {"2015": "green", "tmy": "blue"},
        "dir_norm_rad": {"2015": "orange", "tmy": "red"},
    }

    # def get_hour_avg(year, met_v, month, past_hours, rad):
    #     hours_in_month = (
    #             calendar.monthrange(rad[year].meteo_vars.index[0].year, month + 1)[1] * 24
    #     )
    #     month_slice = slice(past_hours, past_hours + hours_in_month, 1)
    #     past_hours += hours_in_month
    #     all_hours = [[] for i in range(24)]
    #     for h, r in enumerate(rad[year].meteo_vars[met_v][month_slice]):
    #         hour = (h + 1) % 24
    #         all_hours[hour - 1].append(r)
    #     hour_avg = [np.mean(all_hours[i]) for i in range(24)]
    #
    #     return hour_avg

    avg_hours = {
        year: rad[year]
        .meteo_vars.groupby(
            [rad[year].meteo_vars.index.month, rad[year].meteo_vars.index.hour]
        )
        .mean()
        for year in rad.keys()
    }
    for key in avg_hours.keys():
        avg_hours[key].index.names = ["month", "hour"]

    for met_vars in [["glob_hor_rad"], ["dir_norm_rad", "diff_hor_rad"]]:

        # Monthly irradiation
        # Global horizontal irradiation
        # fig, axes = plt.subplots(6, 2, sharex=True, sharey=True, figsize=(3.3, 6))
        fig = plt.figure(figsize=(3.3, 5.3), tight_layout=True)
        gs1 = fig.add_gridspec(6, 2, left=0.18, bottom=0.08, wspace=0.15, hspace=0.45)
        if "glob_hor_rad" in met_vars:
            ylabel = "Global horizontal irradiation"
        else:
            ylabel = "Solar irradiation"

        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["mathtext.fontset"] = "cm"
        fig.supxlabel("Hour of the day [h]", **csfont, **fontsz)
        fig.supylabel(rf"{ylabel} [kWh/$\rm m^2$]", **csfont, **fontsz)

        for month in range(12):
            # x = int(np.floor(month / 2))
            # y = month % 2
            ax = plt.subplot(gs1[month])
            ax.grid(alpha=0.2)

            for met_v in met_vars:
                for year in rad.keys():
                    # get all indices
                    if year == "2015":
                        linestyle = "-"
                    else:
                        linestyle = "--"
                    if met_v == "glob_hor_rad":
                        if year == "tmy":
                            label = "typical meteorological year"
                        else:
                            label = year
                    elif met_v == "diff_hor_rad":
                        if year == "tmy":
                            label = "TMY: Diffuse horizontal irradiation"
                        else:
                            label = "2015: Diffuse horizontal irradiation"
                    else:
                        if year == "tmy":
                            label = "TMY: Direct normal irradiation"
                        else:
                            label = "2015: Direct normal irradiation"

                    ax.plot(
                        range(24),
                        avg_hours[year].loc[month + 1, met_v],
                        color=color[met_v][year],
                        linestyle=linestyle,
                        label=label,
                    )

            ax.set_xticks(np.arange(0, 30, 6).astype(int))
            ax.set_xlim([0, 24])
            if "glob_hor_rad" in met_vars:
                ax.set_yticks(np.arange(0, 1250, 250).astype(int))
                ax.set_ylim([0, 1000])
                ml_y = MultipleLocator(125)

            else:
                ax.set_yticks(np.arange(0, 1000, 200).astype(int))
                ax.set_ylim([0, 800])
                ml_y = MultipleLocator(100)

            ax.yaxis.set_minor_locator(ml_y)
            ml_x = MultipleLocator(3)
            ax.xaxis.set_minor_locator(ml_x)
            ax.set_title(calendar.month_name[month + 1], **csfont, **fontsz)
            ax.set_xticklabels(ax.get_xticks(), **csfont, **fontsz)
            ax.set_yticklabels(ax.get_yticks(), **csfont, **fontsz)

            if month not in (10, 11):
                for tic in ax.xaxis.get_major_ticks():
                    tic.label1.set_visible(False)
                    tic.label2.set_visible(False)

            if month % 2 != 0:
                for tic in ax.yaxis.get_major_ticks():
                    tic.label1.set_visible(False)
                    tic.label2.set_visible(False)

        # if "diff_hor_rad" in met_vars:
        #     gs1.update(top=0.84)
        # else:
        #     gs1.update(top=0.88)

        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)

        handles, labels = ax.get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            bbox_to_anchor=(0.55, 1.03),
            loc="lower center",
            labelspacing=0.1,
        )
        plt.savefig(
            f"datafiles/figures/monthly_irradiation_{met_vars}.svg",
            bbox_inches="tight",
            pad_inches=0,
        )


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
