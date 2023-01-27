#pragma once
#include "Electrochemistry.h"
class Saltbinding :
    public Cell
{
private:

    enum salt
    {
        salt_salt,
        salt_toSalt,
        salt_toIons
    };

    double m_captureCoefficient{};
    double m_dissociationCoeffecient{};
    double m_dt{};

    array_type m_salt;


public:

    Saltbinding(settings_array& settings);
    void react();
    virtual void initializeConcentrations(double contaminantConcentration);
    virtual void calculateCurrents();
    virtual void midSave(std::ofstream& midf);
};

