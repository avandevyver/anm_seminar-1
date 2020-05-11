In order to compile the project, all that needs to be done is to execute cmake
at the top level directory.

Concerning the simulations, only 2 varialbles are supposed to be change:

- The variable TYPE, allowing to choose which of the 3 cases will be computed.
  It can be equal to FLAT_PLATE for the heat diffusion of a flat plate, 
  LID_DRIVEN for the lid driven cavity problem, or CONCENTRIC_CIRCLES for
  the more original problem explained in the report.

- The variable TO_DISPLAY, to decide which distribution shoud be displayed.
  It can be equal to DISPLAY_PRESSURE, DISPLAY_DENSITY, DISPLAY_VELOCITY or
  DISPLAY_TEMPERATURE, displaying of course the chosen distribution.