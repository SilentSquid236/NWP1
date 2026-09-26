# References

Work by other people that this project uses: methods, formulas, data,
software, and designs it imitates. The list is in the format of the
American Meteorological Society (AMS) journals. The list is alphabetical by
first author. Text in the code and docs cites these works author–year,
also in AMS style: "Davies (1976)", "(Wicker and Skamarock 2002)",
"Koch et al. (1983)", and "(Miles 1961; Howard 1961)" for more than one.
Citations use "and", never "&".

**What gets a reference.** Anything the code implements from a named
source, and every data set, service, software library and interface
design the project draws on. Textbook physics (hydrostatic balance, the
hypsometric equation, Poisson's equation for theta, the Coriolis
parameter) does not. When in doubt, cite.

**How the entries were checked.** On 26 September 2026 every DOI below
was resolved against its registry record (CrossRef, or DataCite for
ETOPO1 and Herbie). Author initials, year, title, volume and pages were
taken from that record, not from memory. Two corrections to the
registry:
- The CrossRef record for Buizza et al. (1999) misspells the second
  author "Milleer". The published author is M. Miller.
- CrossRef has no DOI for the PyTorch paper. Its volume, editors and
  pages were checked against the NeurIPS proceedings page and dblp.

Web services have no fixed version, so their entries carry the access
date.

To add a reference, give it a DOI where one exists and resolve it before
it goes in. Then add its in-text citation where the work is used, and a
row in the table at the end.

---

Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus form approximation of saturation vapor pressure. *J. Appl. Meteor.*, **35**, 601–609, https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2.

Arakawa, A., and C. S. Konor, 1996: Vertical differencing of the primitive equations based on the Charney–Phillips grid in hybrid σ–*p* vertical coordinates. *Mon. Wea. Rev.*, **124**, 511–528, https://doi.org/10.1175/1520-0493(1996)124<0511:VDOTPE>2.0.CO;2.

Arakawa, A., and V. R. Lamb, 1977: Computational design of the basic dynamical processes of the UCLA general circulation model. *Methods in Computational Physics: Advances in Research and Applications*, Vol. 17, Academic Press, 173–265, https://doi.org/10.1016/B978-0-12-460817-7.50009-4.

Barnes, S. L., 1964: A technique for maximizing details in numerical weather map analysis. *J. Appl. Meteor.*, **3**, 396–409, https://doi.org/10.1175/1520-0450(1964)003<0396:ATFMDI>2.0.CO;2.

Blaylock, B. K., 2026: Herbie: Retrieve numerical weather prediction model data. Zenodo, accessed 26 September 2026, https://doi.org/10.5281/zenodo.4567540.

Blumberg, W. G., K. T. Halbert, T. A. Supinie, P. T. Marsh, R. L. Thompson, and J. A. Hart, 2017: SHARPpy: An open-source sounding analysis toolkit for the atmospheric sciences. *Bull. Amer. Meteor. Soc.*, **98**, 1625–1636, https://doi.org/10.1175/BAMS-D-15-00309.1.

Bolton, D., 1980: The computation of equivalent potential temperature. *Mon. Wea. Rev.*, **108**, 1046–1053, https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2.

Bougeault, P., 1983: A non-reflective upper boundary condition for limited-height hydrostatic models. *Mon. Wea. Rev.*, **111**, 420–429, https://doi.org/10.1175/1520-0493(1983)111<0420:ANRUBC>2.0.CO;2.

Buizza, R., M. Miller, and T. N. Palmer, 1999: Stochastic representation of model uncertainties in the ECMWF ensemble prediction system. *Quart. J. Roy. Meteor. Soc.*, **125**, 2887–2908, https://doi.org/10.1002/qj.49712556006.

Cressman, G. P., 1959: An operational objective analysis system. *Mon. Wea. Rev.*, **87**, 367–374, https://doi.org/10.1175/1520-0493(1959)087<0367:AOOAS>2.0.CO;2.

Davies, H. C., 1976: A lateral boundary formulation for multi-level prediction models. *Quart. J. Roy. Meteor. Soc.*, **102**, 405–418, https://doi.org/10.1002/qj.49710243210.

Dowell, D. C., and Coauthors, 2022: The High-Resolution Rapid Refresh (HRRR): An hourly updating convection-allowing forecast model. Part I: Motivation and system description. *Wea. Forecasting*, **37**, 1371–1395, https://doi.org/10.1175/WAF-D-21-0151.1.

ECMWF, 2026: cfgrib. GitHub repository, accessed 26 September 2026, https://github.com/ecmwf/cfgrib.

Harris, C. R., and Coauthors, 2020: Array programming with NumPy. *Nature*, **585**, 357–362, https://doi.org/10.1038/s41586-020-2649-2.

Howard, L. N., 1961: Note on a paper of John W. Miles. *J. Fluid Mech.*, **10**, 509–512, https://doi.org/10.1017/S0022112061000317.

Hoyer, S., and J. Hamman, 2017: xarray: N-D labeled arrays and datasets in Python. *J. Open Res. Software*, **5**, 10, https://doi.org/10.5334/jors.148.

Hunter, J. D., 2007: Matplotlib: A 2D graphics environment. *Comput. Sci. Eng.*, **9**, 90–95, https://doi.org/10.1109/MCSE.2007.55.

Iowa Environmental Mesonet, 2026a: ASOS-AWOS-METAR data download service. Iowa State University, accessed 26 September 2026, https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py.

Iowa Environmental Mesonet, 2026b: Upper-air sounding (RAOB) download service and RAOB station table. Iowa State University, accessed 26 September 2026, https://mesonet.agron.iastate.edu/cgi-bin/request/raob.py.

Klemp, J. B., and D. R. Durran, 1983: An upper boundary condition permitting internal gravity wave radiation in numerical mesoscale models. *Mon. Wea. Rev.*, **111**, 430–444, https://doi.org/10.1175/1520-0493(1983)111<0430:AUBCPI>2.0.CO;2.

Koch, S. E., M. desJardins, and P. J. Kocin, 1983: An interactive Barnes objective map analysis scheme for use with satellite and conventional data. *J. Climate Appl. Meteor.*, **22**, 1487–1503, https://doi.org/10.1175/1520-0450(1983)022<1487:AIBOMA>2.0.CO;2.

Louis, J.-F., 1979: A parametric model of vertical eddy fluxes in the atmosphere. *Bound.-Layer Meteor.*, **17**, 187–202, https://doi.org/10.1007/BF00117978.

Manabe, S., J. Smagorinsky, and R. F. Strickler, 1965: Simulated climatology of a general circulation model with a hydrologic cycle. *Mon. Wea. Rev.*, **93**, 769–798, https://doi.org/10.1175/1520-0493(1965)093<0769:SCOAGC>2.3.CO;2.

Miles, J. W., 1961: On the stability of heterogeneous shear flows. *J. Fluid Mech.*, **10**, 496–508, https://doi.org/10.1017/S0022112061000305.

Natural Earth, 2026: Natural Earth 1:50m cultural and physical vectors (natural-earth-vector repository). Accessed 26 September 2026, https://github.com/nvkelso/natural-earth-vector.

NOAA, 2026a: NOAA High-Resolution Rapid Refresh (HRRR) model. NOAA Open Data Dissemination, Registry of Open Data on AWS, accessed 26 September 2026, https://registry.opendata.aws/noaa-hrrr-pds/.

NOAA, 2026b: NOAA Multi-Radar/Multi-Sensor System (MRMS). NOAA Open Data Dissemination, Registry of Open Data on AWS, accessed 26 September 2026, https://registry.opendata.aws/noaa-mrms-pds/.

NOAA National Data Buoy Center, 2026: Active station list and 5-day standard meteorological data files. Accessed 26 September 2026, https://www.ndbc.noaa.gov/.

NOAA National Geophysical Data Center, 2009: ETOPO1 1 arc-minute global relief model. NOAA National Centers for Environmental Information, accessed 26 September 2026 through the NOAA CoastWatch ERDDAP server (dataset `etopo180`, https://coastwatch.pfeg.noaa.gov/erddap/), https://doi.org/10.7289/V5C8276M.

Paszke, A., and Coauthors, 2019: PyTorch: An imperative style, high-performance deep learning library. *Advances in Neural Information Processing Systems 32*, H. Wallach et al., Eds., Curran Associates, 8024–8035.

Phillips, N. A., 1957: A coordinate system having some special advantages for numerical forecasting. *J. Meteor.*, **14**, 184–185, https://doi.org/10.1175/1520-0469(1957)014<0184:ACSHSS>2.0.CO;2.

Pillow Contributors, 2026: Pillow: The friendly PIL fork. Accessed 26 September 2026, https://python-pillow.org/.

Pivotal Weather, 2026: Pivotal Weather forecast model maps. Accessed 26 September 2026, https://www.pivotalweather.com/.

Sadourny, R., 1975: The dynamics of finite-difference models of the shallow-water equations. *J. Atmos. Sci.*, **32**, 680–689, https://doi.org/10.1175/1520-0469(1975)032<0680:TDOFDM>2.0.CO;2.

Simmons, A. J., and D. M. Burridge, 1981: An energy and angular-momentum conserving vertical finite-difference scheme and hybrid vertical coordinates. *Mon. Wea. Rev.*, **109**, 758–766, https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2.

Skamarock, W. C., and J. B. Klemp, 1992: The stability of time-split numerical methods for the hydrostatic and the nonhydrostatic elastic equations. *Mon. Wea. Rev.*, **120**, 2109–2127, https://doi.org/10.1175/1520-0493(1992)120<2109:TSOTSN>2.0.CO;2.

Snyder, J. P., 1987: *Map Projections—A Working Manual*. U.S. Geological Survey Professional Paper 1395, https://doi.org/10.3133/pp1395.

Toth, Z., and E. Kalnay, 1993: Ensemble forecasting at NMC: The generation of perturbations. *Bull. Amer. Meteor. Soc.*, **74**, 2317–2330, https://doi.org/10.1175/1520-0477(1993)074<2317:EFANTG>2.0.CO;2.

Wicker, L. J., and W. C. Skamarock, 2002: Time-splitting methods for elastic models using forward time schemes. *Mon. Wea. Rev.*, **130**, 2088–2097, https://doi.org/10.1175/1520-0493(2002)130<2088:TSMFEM>2.0.CO;2.

Zhang, J., and Coauthors, 2016: Multi-Radar Multi-Sensor (MRMS) quantitative precipitation estimation: Initial operating capabilities. *Bull. Amer. Meteor. Soc.*, **97**, 621–638, https://doi.org/10.1175/BAMS-D-14-00174.1.

---

## Where each is used

| Reference | Used for | Where |
|---|---|---|
| Alduchov and Eskridge (1996) | Magnus constants 17.625 and 243.04 °C for RH from dewpoint | `src/verification/fetchers.py` |
| Arakawa and Konor (1996) | the vertical computational mode of the Lorenz grid (u, v and theta on the same levels), a P-60 candidate | `docs/RESEARCH_LOG.md`, `docs/PROBLEMS.md` |
| Arakawa and Lamb (1977) | C-grid staggering | `src/dynamics/grid.py` |
| Barnes (1964); Koch et al. (1983) | successive-correction analysis, multi-pass with convergence parameter gamma | `src/analysis/barnes.py` |
| Blaylock (2026) | HRRR download (seeding and diagnosis only) | `src/ingest_hrrr.py`, `diagnose_herbie.py` |
| Blumberg et al. (2017) | layout of the click-for-sounding panel (skew-T and hodograph) | `src/maps/viewer.py` |
| Bolton (1980) | saturation vapour pressure | `src/analysis/build.py` |
| Bougeault (1983); Klemp and Durran (1983) | radiating upper boundary (optional, not ported to torch) | `src/dynamics/radiation.py` |
| Buizza et al. (1999) | SPPT-style stochastic tendency perturbations (optional) | `src/dynamics/subgrid.py`, `src/forecast.py` |
| Cressman (1959) | background to the analysis design | `docs/DATA_ASSIMILATION.md` |
| Davies (1976) | lateral boundary relaxation zone | `src/dynamics/boundaries.py`, `src/forecast.py` |
| Dowell et al. (2022); NOAA (2026a) | HRRR, which may seed a forecast but never verify one | `src/ingest_hrrr.py` |
| ECMWF (2026); Hoyer and Hamman (2017) | reading HRRR GRIB2 | `src/ingest_hrrr.py` |
| Harris et al. (2020) | NumPy, the numerical base of the whole code | throughout |
| Hunter (2007) | Matplotlib, all maps and figures | `src/maps/render.py`, `src/make_maps.py` |
| Iowa Environmental Mesonet (2026a, b) | surface (ASOS) and upper-air (RAOB) observations | `src/analysis/sources.py`, `src/verification/fetchers.py` |
| Louis (1979) | shape of the Richardson-number stability functions | `src/dynamics/turbulence.py`, `src/dynamics/surface.py` |
| Manabe et al. (1965) | dry convective adjustment | `src/dynamics/convection.py` |
| Miles (1961); Howard (1961) | Ri < 0.25 as the shear-instability criterion; the mixing threshold | `src/dynamics/turbulence.py`, P-60 |
| Natural Earth (2026) | state, coast and lake lines on the maps | `src/maps/geography.py` |
| NOAA (2026b); Zhang et al. (2016) | MRMS radar products | `src/analysis/sources.py`, `src/verification/fetchers.py` |
| NOAA National Data Buoy Center (2026) | buoy and C-MAN observations | `src/analysis/sources.py` |
| NOAA National Geophysical Data Center (2009) | model terrain (ETOPO1), fetched through ERDDAP | `src/analysis/geo.py` |
| Paszke et al. (2019) | PyTorch backend | `src/dynamics/backend.py` |
| Phillips (1957) | terrain-following sigma coordinate | `src/dynamics/sigma.py` |
| Pillow Contributors (2026) | reading rendered PNGs in the map tests | `src/maps/test_maps.py` |
| Pivotal Weather (2026) | product names, units and layout (titles with UTC and US Eastern times) that the maps imitate | `src/maps/render.py` |
| Sadourny (1975) | vector-invariant momentum form of the shallow-water test core | `src/dynamics/shallow_water.py` |
| Simmons and Burridge (1981) | hydrostatically consistent vertical discretisation | `src/dynamics/sigma.py` |
| Skamarock and Klemp (1992) | divergence damping (optional, `--div-damp`; P-60 test V) | `src/dynamics/subgrid.py`, `src/forecast.py` |
| Snyder (1987) | Lambert conformal conic formulas | `src/maps/geography.py` |
| Toth and Kalnay (1993) | bred-vector idea behind the round-off difference mode | `tools/mode_structure.py`, `tools/mode_budget.py` |
| Wicker and Skamarock (2002) | third-order Runge–Kutta time stepping | `src/dynamics/shallow_water.py`, `src/dynamics/primitive_sigma.py` |
