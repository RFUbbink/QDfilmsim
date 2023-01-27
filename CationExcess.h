#pragma once
#include "Electrochemistry.h"

class CationExcess :
    public Cell
{

	enum cations
	{
		cations_toExcess,
		cations_excess,
		cations_current
	};
private:
    array_type m_cations{};

public:
	CationExcess(settings_array& settings);

	virtual void calculatePotentialProfile();
	void react();
	virtual void calculateCurrents();
	virtual void updateConcentrations();
	virtual void midSave(std::ofstream& midf);
};

