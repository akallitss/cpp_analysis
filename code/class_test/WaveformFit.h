//
// Created by dylan on 10/20/24.
//

#ifndef WAVEFORM_FIT_H
#define WAVEFORM_FIT_H

#include "TObject.h"

class WaveformFit : public TObject {
public:
    float amplitude;
    float time;
    float chi2;

    WaveformFit();  // Default constructor
    WaveformFit(float amp, float t, float chi);  // Parameterized constructor
    virtual ~WaveformFit();  // Destructor

    ClassDef(WaveformFit, 1);  // ROOT's dictionary macro
};

#endif
//WAVEFORMFIT_H
