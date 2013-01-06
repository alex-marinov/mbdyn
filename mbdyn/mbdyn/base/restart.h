/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Macro per inserire le funzioni di restart nelle classi che le richiedono */

/* Uso: in ogni classe di base che partecipa alla scrittura del file 
 * di restart mettere la dichiarazione:
 * 
 *     DECLARE_RESTART_VIRTUAL;
 * 
 * (attenzione al punto e virgola finale)
 * In ogni classe derivata che implementa una sua funzione di restart
 * mettere la dichiarazione:
 * 
 *     DECLARE_RESTART;
 * 
 * quindi nell'implementazione della classe aggiungere la riga:
 * 
 *     IMPLEMENT_RESTART( nome_classe )
 *     {
 *        return out [ << dati in uscita] ;
 *     }
 * 
 * Si fara' in modo di mantenere l'oggetto si cui si scrive, 'out',
 * compatibile con l'operatore di scorrimento a sinistra 
 */
       

#ifndef RESTART_H
#define RESTART_H


#define DECLARE_RESTART_VIRTUAL \
  virtual ostream& Restart(ostream& out) const = 0

#define DECLARE_RESTART \
  virtual ostream& Restart(ostream& out) const

#define IMPLEMENT_RESTART( ClassName ) \
  ostream& ClassName::Restart(ostream& out) const
  
#endif
