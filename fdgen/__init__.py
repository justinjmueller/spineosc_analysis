"""fdgen: blind fake-data weight generator driven by a YAML scenario file.

Per MC event, a scenario is a ratio to the PROfit prediction: the two-flavour oscillation probability (evaluated
exactly the way PROfit evaluates it) times any number of dial weights read from the CVMC friend trees. Optional
Poisson counts follow PROfit's per-bin Poisson variation. See README.md.
"""
