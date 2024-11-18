import os
from pathlib import Path
from sharc.results import Results
from sharc.post_processor import PostProcessor
import plotly.graph_objects as go
from sharc.parameters.parameters import Parameters

import glob
from sharc.antenna.antenna_s465 import AntennaS465
from sharc.antenna.antenna_s580 import AntennaS580
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt, PlotAntennaPattern

campaign_base_dir = str((Path(__file__) / ".." / "..").resolve())
dl_dir = os.path.join(campaign_base_dir, "output_dl_16x8")

post_processor = PostProcessor()

# Add a legend to results in folder that match the pattern
# This could easily come from a config file
post_processor\
    .add_plot_legend_pattern(
        dir_name_contains="sat_jpss_1_0m",
        legend="ES for Sat. JPSS-1 (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_1_10000m",
        legend="ES for Sat. JPSS-1 (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_1_15000m",
        legend="ES for Sat. JPSS-1 (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_1_5000m",
        legend="ES for Sat. JPSS-1 (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_1_20000m",
        legend="ES for Sat. JPSS-1 (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_2_0m",
        legend="ES for Sat. JPSS-2 (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_2_10000m",
        legend="ES for Sat. JPSS-2 (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_2_15000m",
        legend="ES for Sat. JPSS-2 (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_2_5000m",
        legend="ES for Sat. JPSS-2 (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="sat_jpss_2_20000m",
        legend="ES for Sat. JPSS-2 (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_F_0m",
        legend="ES for Sat. F (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_F_10000m",
        legend="ES for Sat. F (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_F_15000m",
        legend="ES for Sat. F (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_F_5000m",
        legend="ES for Sat. F (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_F_20000m",
        legend="ES for Sat. F (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7820_0m",
        legend="ES for Sat. E and G (7820 MHz) (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7820_10000m",
        legend="ES for Sat. E and G (7820 MHz) (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7820_15000m",
        legend="ES for Sat. E and G (7820 MHz) (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7820_5000m",
        legend="ES for Sat. E and G (7820 MHz) (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7820_20000m",
        legend="ES for Sat. E and G (7820 MHz) (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7780_0m",
        legend="ES for Sat. E and G (7780 MHz) (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7780_10000m",
        legend="ES for Sat. E and G (7780 MHz) (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7780_15000m",
        legend="ES for Sat. E and G (7780 MHz) (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7780_5000m",
        legend="ES for Sat. E and G (7780 MHz) (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_E_and_G_7780_20000m",
        legend="ES for Sat. E and G (7780 MHz) (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_C_and_S_0m",
        legend="ES for Sat. C and S (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_C_and_S_10000m",
        legend="ES for Sat. C and S (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_C_and_S_15000m",
        legend="ES for Sat. C and S (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_C_and_S_5000m",
        legend="ES for Sat. C and S (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_C_and_S_20000m",
        legend="ES for Sat. C and S (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_AN_0m",
        legend="ES for Sat. AN (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_AN_10000m",
        legend="ES for Sat. AN (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_AN_15000m",
        legend="ES for Sat. AN (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_AN_5000m",
        legend="ES for Sat. AN (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_AN_20000m",
        legend="ES for Sat. AN (20km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_BF_0m",
        legend="ES for Sat. BF (0m)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_BF_10000m",
        legend="ES for Sat. BF (10km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_BF_15000m",
        legend="ES for Sat. BF (15km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_BF_5000m",
        legend="ES for Sat. BF (5km)"
    ).add_plot_legend_pattern(
        dir_name_contains="non_gso_sat_BF_20000m",
        legend="ES for Sat. BF (20km)"
    )

attributes_to_plot = [
    "system_imt_antenna_gain",
    "imt_system_path_loss",
    "imt_system_antenna_gain",
    "system_dl_interf_power_per_mhz",
    "system_ul_interf_power_per_mhz",
]

dl_results = Results.load_many_from_dir(dl_dir, only_latest=True, only_samples=attributes_to_plot)
ul_results = Results.load_many_from_dir(os.path.join(campaign_base_dir, "output_ul"), only_latest=True, only_samples=attributes_to_plot)
# ^: typing.List[Results]

all_results = [*dl_results, *ul_results]

post_processor.add_results(all_results)

plots = post_processor.generate_cdf_plots_from_results(
    all_results
)

post_processor.add_plots(plots)

# Show a single plot:
system_dl_interf_power_plot = post_processor\
    .get_plot_by_results_attribute_name("system_dl_interf_power_per_mhz")

system_ul_interf_power_plot = post_processor\
    .get_plot_by_results_attribute_name("system_ul_interf_power_per_mhz")


# Add a protection criteria line:
# dB to dBm (+ 30)
# the following conversion makes the criteria more strict, so there may not be a problem
# / 10MHz to / 1MHz  (-10)
interf_protection_criteria = -148 + 30 - 10
if system_dl_interf_power_plot:
    system_dl_interf_power_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="20% criteria"
    )

if system_ul_interf_power_plot:
    system_ul_interf_power_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="20% criteria"
    )

interf_protection_criteria = -127 + 30 - 10

if system_dl_interf_power_plot:
    system_dl_interf_power_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="0.0016% criteria"
    )

if system_ul_interf_power_plot:
    system_ul_interf_power_plot.add_vline(
        interf_protection_criteria, line_dash="dash",
        name="0.0016% criteria"
    )

# # Aggregate results.
# # TODO: Need to update building entry loss before simulating UL

# aggregated_plot = None
# if system_ul_interf_power_plot and system_dl_interf_power_plot:
#     aggregated_plot = go.Figure()

#     for dl_r in dl_results:
#         legend1 = post_processor.get_results_possible_legends(dl_r)[0]
#         ul_r = None
#         for maybe in ul_results:
#             legend2 = post_processor.get_results_possible_legends(maybe)[0]
#             if legend1 == legend2:
#                 ul_r = maybe
#                 break
#         if ul_r is None:
#             # raise Exception(f"Cannot aggregate {legend1} and {legend2}")
#             continue
            
#         aggregated_results = PostProcessor.aggregate_results(
#             dl_samples=dl_r.system_dl_interf_power,
#             ul_samples=ul_r.system_ul_interf_power,
#             ul_tdd_factor=0.25,
#             n_bs_sim=1,
#             n_bs_actual=1
#         )
#         x, y = PostProcessor.cdf_from(aggregated_results)

#         aggregated_plot.add_trace(
#             go.Scatter(x=x, y=y, mode='lines', name=f'{legend1["legend"]}',),
#         )

#     # Add a protection criteria line:
#     # dB to dBm (+ 30)
#     # the following conversion makes the criteria more strict, so there may not be a problem
#     interf_protection_criteria = -154 + 30

#     aggregated_plot.add_vline(
#         interf_protection_criteria, line_dash="dash",
#         name="1% criteria"
#     )

# this is not optimal at all, but it works, so...
antenna_legend_start = "Radiation Pattern of "
i = 0

antenna_bs = None
antenna_ue = None
for result in all_results:
    if "_5000m_" not in result.output_directory:
        continue

    params_file = glob.glob(result.output_directory + "/*.yaml")[0]
    params = Parameters()
    params.set_file_name(params_file)
    params.read_params()

    # TODO: use antenna factory here if it ever exists
    legend = post_processor.get_results_possible_legends(result)[0]
    if params.single_earth_station.antenna.pattern == "ITU-R S.465":
        antenna = AntennaS465(params.single_earth_station.antenna.itu_r_s_465)
    if params.single_earth_station.antenna.pattern == "ITU-R S.580":
        antenna = AntennaS580(params.single_earth_station.antenna.itu_r_s_580)
    PostProcessor.generate_antenna_radiation_pattern_plot(
        antenna, legend["legend"]
    ).show()
    if i == 0:
        antenna_bs = AntennaBeamformingImt(
            params.imt.bs.antenna.get_antenna_parameters(),
            0,
            0,
            # -params.imt.bs.antenna.downtilt
        )
        antenna_ue = AntennaBeamformingImt(
            params.imt.ue.antenna.get_antenna_parameters(),
            0,
            0
        )

    i += 1


# Show a single plot:
for pl_name in attributes_to_plot:
    plot = post_processor\
        .get_plot_by_results_attribute_name(pl_name)
    if plot:
        plot.show()

# if aggregated_plot:
#     aggregated_plot.show()


plot_antenna_imt = PlotAntennaPattern("")

# Plot BS TX radiation patterns
f = plot_antenna_imt.plot_element_pattern(antenna_bs, "BS", "ELEMENT")
# f.savefig(figs_dir + "BS_element.pdf", bbox_inches='tight')
f = plot_antenna_imt.plot_element_pattern(antenna_bs, "TX", "ARRAY")
# f.savefig(figs_dir + "BS_array.pdf", bbox_inches='tight')

# Plot UE TX radiation patterns
plot_antenna_imt.plot_element_pattern(antenna_ue, "UE", "ELEMENT")
plot_antenna_imt.plot_element_pattern(antenna_ue, "UE", "ARRAY")

# Plot every plot:
# for plot in plots:
#     plot.show()

full_results = ""

for result in all_results:
    # This generates the mean, median, variance, etc
    stats = PostProcessor.generate_statistics(
        result=result
    ).write_to_results_dir()

    full_results += str(stats) + "\n"
    # # do whatever you want here:
    # if "fspl_45deg" in stats.results_output_dir:
    #     get some stat and do something

with open(dl_dir + "/stats.txt", "w") as f:
    f.write(full_results)
