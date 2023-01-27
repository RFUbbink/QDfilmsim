#pragma once
#include "Electrochemistry.h"

class Freezable :
    public Cell
{
private:
    std::vector<double> m_frozencations;
    std::vector<double> m_frozenanions;

public:
    Freezable(Cell& cell);

    virtual void calculatePotentialProfile();
    virtual void updateConcentrations();
    void setIonMobility(double dt);
    void fixIons();
};

