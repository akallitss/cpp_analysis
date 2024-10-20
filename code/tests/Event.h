#ifndef EVENT_H
#define EVENT_H

#include <TObject.h>

class Event : public TObject {
public:
    float amplitude;
    float start_time;
    float charge;

    Event() : amplitude(0), start_time(0), charge(0) {}
    virtual ~Event() {}

    void SetValues(float amp, float start, float chg) {
        amplitude = amp;
        start_time = start;
        charge = chg;
    }

    ClassDef(Event, 1);  // Version 1 of the class
};

#endif

