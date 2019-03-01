/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Elementi termici */

#include "mbconfig.h"		/* This goes first in every *.c,*.cc file */

#include "mynewmem.h"
#include "therm.h"
#include "thermalnode.h"
#include "thermalcapacitance.h"
#include "thermalresistance.h"
#include "thermalsource.h"
#include "dataman.h"

/* Thermal - begin */

Thermal::Thermal(unsigned int uL,
			const DofOwner* pDO, flag fOut)
: Elem(uL, fOut),
ElemWithDofs(uL, pDO, fOut)
{
	NO_OP;
}


Thermal::~Thermal(void)
{
	NO_OP;
}


/* Contributo al file di restart
 * (Nota: e' incompleta, deve essere chiamata dalla funzione corrispndente
 * relativa alla classe derivata */
std::ostream& Thermal::Restart(std::ostream& out) const {
	return out << "  thermal: " << GetLabel();
}


/* Tipo dell'elemento (usato solo per debug ecc.) */
Elem::Type Thermal::GetElemType(void) const
{
	return Elem::THERMAL;
}

/* Electric - end */





/* Legge un elemento elettrico */

Elem* ReadThermal(DataManager* pDM,
			MBDynParser& HP,
			const DofOwner* pDO,
			unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadEThermal()");

	const char* sKeyWords[] = {
		"resistance",
		"capacitance",
		"source"
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		THERMALRESISTANCE = 0,
		THERMALCAPACITANCE,
		THERMALSOURCE,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di elemento elettrico */
	KeyWords CurrKeyWord = KeyWords(HP.GetWord());

#ifdef DEBUG
	if (CurrKeyWord >= 0) {
		std::cout << "thermal element type: "
	<< sKeyWords[CurrKeyWord] << std::endl;
	}
#endif

	Elem* pEl = 0;

	switch (CurrKeyWord) {
		/*  */

		case THERMALRESISTANCE: {
			const ThermalNode* pThNode1 = pDM->ReadNode<const ThermalNode, Node::THERMAL>(HP);
			const ThermalNode* pThNode2 = pDM->ReadNode<const ThermalNode, Node::THERMAL>(HP);
			doublereal r = HP.GetReal();
			flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
			SAFENEWWITHCONSTRUCTOR(pEl,
				ThermalResistance,
				ThermalResistance(uLabel, pDO,
						pThNode1, pThNode2, 
						r, fOut));

			std::ostream& out = pDM->GetLogFile();
			out << "thermal resistance: " << uLabel
				<< " " << pThNode1->GetLabel()
				<< " " << pThNode2->GetLabel()
				<< std::endl;
			break;
		}

		case THERMALCAPACITANCE: {
			const ThermalNode* pThNode1 = pDM->ReadNode<const ThermalNode, Node::THERMAL>(HP);
			doublereal c = HP.GetReal();
			flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
			SAFENEWWITHCONSTRUCTOR(pEl,
				ThermalCapacitance,
				ThermalCapacitance(uLabel, pDO, pThNode1, c, fOut));

			std::ostream& out = pDM->GetLogFile();
			out << "thermal capacitance: " << uLabel
				<< " " << pThNode1->GetLabel()
				<< std::endl;
			break;
		}

		case THERMALSOURCE: {
			const ThermalNode* pThNode1 = pDM->ReadNode<const ThermalNode, Node::THERMAL>(HP);
			DriveCaller* pDC = HP.GetDriveCaller();
			flag fOut = pDM->fReadOutput(HP, Elem::ELECTRIC);
			SAFENEWWITHCONSTRUCTOR(pEl,
				ThermalSource,
				ThermalSource(uLabel, pDO, pThNode1, pDC, fOut));

			std::ostream& out = pDM->GetLogFile();
			out << "thermal source: " << uLabel
				<< " " << pThNode1->GetLabel()
				<< std::endl;
			break;
		}

		/* Aggiungere altri elementi elettrici */

		default: {
			silent_cerr("unknown thermal element type in thermal element " << uLabel
				<< " at line " << HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
} /* ReadThermal() */

