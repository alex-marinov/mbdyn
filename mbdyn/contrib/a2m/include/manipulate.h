#ifndef MANIPULATE_H
#define MANIPULATE_H

#include <cards.h>
#include <storage.h>

void Ref_Interp (MBDyn_reference*, Id&, Id&);
void Ref_Interp (MBDyn_reference*, MBDyn_reference*, MBDyn_reference*);
void NodeSet (MBDyn_node_structural*, Id);
void NodeSet (MBDyn_node_structural*, MBDyn_reference*);

#endif
