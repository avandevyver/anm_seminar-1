#ifndef SHARED_VARIABLES_H
#define SHARED_VARIABLES_H

#define DISPLAY_PRESSURE 1
#define DISPLAY_DENSITY 2
#define DISPLAY_VELOCITY 3
#define DISPLAY_TEMPERATURE 4

#define USE_PREDICTOR 0
#define USE_EULER 1

#define FLAT_PLATE 1
#define LID_DRIVEN 2
#define CONCENTRIC_CIRCLES 3

#define TO_DISPLAY DISPLAY_VELOCITY
#define TYPE CONCENTRIC_CIRCLES
#define INTEGRATION_METHOD USE_PREDICTOR

#define MASS 5
#define VELOCITY 30
#define INITDENSITY 0.5
#define ID 0.5


//for TEMPERATURE TYPE == 1
#if TYPE == 1
#define NPTS 2500
#define DT 2
#define CUT_OFF 2.5
#define MU 75
#define PSY 100
#define R_0 5.5
#define k 1
#define alph 60
#endif

//for VELOCITY TYPE == 2
#if TYPE == 2
#define NPTS 2500
#define DT 0.1
#define CUT_OFF 5.5
#define MU 1000 
#define PSY 60
#define R_0 5.5
#define k 1
#define alph 60
#endif

//for CIRCULAR TYPE == 3
#if TYPE == 3
#define NPTS 225
#define DT 0.1
#define  CUT_OFF 2.5
#define MU 75
#define PSY 100
#define R_0 5.5
#define k 1
#define alph 60
#endif
#endif