# TWDTW_codes

## Time-Weighted Dynamic Time Warping for satellite image time series analysis

The TWDTW codes provides a full calculation of the Time-Weighted Dynamic Time Warping (TWDTW) method for maize mapping using fused NDVI time series (Maus et al. 2016). TWDTW codes provides full cycle of maize classification using image time series, ranging from calculating the distance between two series to classifying through the threshold (Peng, et al. 2023).

## Instruction：

(1) Use 'distance_calculate.py' first to obtain the distance file of each province, then use 'maize_classified.jl' to obtain the classified maize file through the distance file and threshold of each province.

(2) It should be noted that 'twdtw.f90' is the main code for calculating the distance, which should be called by 'distance_calculate.py'. Besides, 'recognize_core.jl' is the main code for classifying maize by the distance file, which should be called by 'maize_classified.jl'.

## References

Maus, Victor, Gilberto Camara, Ricardo Cartaxo, Alber Sanchez, Fernando M. Ramos, and Gilberto R. de Queiroz. 2016. “A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.” IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing 9 (8): 3729–39. https://doi.org/10.1109/JSTARS.2016.2517118.

Peng Q, Shen R, Dong J, Han W, Huang J, Ye T, Zhao W and Yuan W (2023), A new method for classifyingmaize by combining the phenological information of multiple satellite-based spectral bands. Front. Environ. Sci. 10:1089007. https://doi:10.3389/fenvs.2022.1089007
