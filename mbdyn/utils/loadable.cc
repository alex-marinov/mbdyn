#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <ctype.h>
#include <string.h>

struct list {
   	list* next;
   	list* prev;
   	char* p;
    list(void) : next(NULL), prev(NULL), p(NULL) {};
};

int eat_spaces(istream& in) 
{
    int c;
 	while (isspace(c = in.get()));
   	in.putback(c);
   	if (!in || c == EOF) {
	   	return EOF;
	}
   	return c;
}

int main(int argn, const char* const argv[]) 
{
    char *demangled_fname = "demangled.h";
   	char *mangled_fname = "mangled.tmp";
   
   	ostream* pout = &(ostream&)cout;
   
   	if (argn > 1) {
	   	demangled_fname = (char*)argv[1];
	}
   	if (argn > 2) {
	   	mangled_fname = (char*)argv[2];
	}
   
   	if (argn > 3) {
	   	pout = new ofstream((char*)argv[3]);
	}
   
   	ifstream demangled(demangled_fname);
   	ifstream mangled(mangled_fname);
   	
   	if (!demangled || !mangled) {
	   	cerr << "unable to open demangled/mangled file(s)" << endl;
	   	exit(EXIT_FAILURE);
	}
   
    const int buflen = 1024;
   	char buf[buflen];   	

   	list* d_start = new list;
   	list* d_end = new list;
   
   	d_start->next = d_end;
   	d_end->prev = d_start;
   
   	do {
 	  	int c = eat_spaces(demangled);
	   	if (c == EOF) {
		   	break;
		}
 	   	demangled.getline(buf, buflen);	   	   
	    if (buf[0] == '\0') {
		   	cerr << "error: non-null string expected" << endl;
		   	break;
		}
	   	char *sep = strrchr(buf, ',');
	   	if (sep == NULL) {
		    sep = strrchr(buf, '"');
		   	if (sep == NULL) {
			   	sep = buf+strlen(buf);
			}
		} else if (sep[-1] == '"') {		   
		    sep = sep-1;
		}
	   	sep[0] = '\0';
	   	list *p = new list;
	   	p->p = new char[strlen(buf+(buf[0] == '"' ? 1 : 0))];
	   	strcpy(p->p, buf+(buf[0] == '"' ? 1 : 0));
	   	p->prev = d_end->prev;
	   	d_end->prev->next = p;
	   	d_end->prev = p;
	   	p->next = d_end;	 
	} while (demangled);
   	if (strcmp(d_end->prev->p, "NULL") != 0) {
	   	cerr << "warning: NULL expected as last demangled function" << endl;
	} else {
	   	list* d_tmp = d_end->prev->prev;
	   	d_tmp->next = d_end;
	   	delete[] d_end->prev->p;
	   	delete d_end->prev;
	   	d_end->prev = d_tmp;
	}

   	list* m_start = new list;
   	list* m_end = new list;
   
   	m_start->next = m_end;
   	m_end->prev = m_start;
   
   	do {
 	  	int c = eat_spaces(mangled);
	   	if (c == EOF) {
		   	break;
		}
 	   	mangled.getline(buf, buflen);	 
	    if (buf[0] == '\0') {
		   	cerr << "error: non-null string expected" << endl;
		   	break;
		}
	   	list *p = new list;
	    p->p = new char[strlen(buf)+1];
	   	strcpy(p->p, buf);
	   	p->prev = m_end->prev;
	   	m_end->prev->next = p;
	   	m_end->prev = p;
	   	p->next = m_end;	 
	} while (mangled);
      	   
   	list* d_l = d_start->next;
   	list* m_curr = m_start;
   	while (d_l->next != NULL) {
	   	list* m_l = m_curr->next;	   	
	   	while (m_l->next != NULL) {
		   	if (strncmp(d_l->p, m_l->p, strlen(d_l->p)) == 0) {
		   		m_l->prev->next = m_l->next;
		   		m_l->next->prev = m_l->prev;
		   
		   		m_l->prev = m_curr;
		   		m_l->next = m_curr->next;
		   		m_curr->next->prev = m_l;			 
		   		m_curr->next = m_l;
			    m_curr = m_l;
			   	break;			   
			}
		   	m_l = m_l->next;
		}
	   	if (m_l->next == NULL) {
		   	cerr << "error: \"" << d_l->p << "\" not found in mangled!" << endl;
		   	exit(EXIT_FAILURE);
		} 
	   	d_l = d_l->next;
	}
	if (m_curr->next != m_end) {
  		cerr << "warning: extra functions in mangled file!" << endl;
	}

   	list* l = d_start->next;
    /*
   	while (l->next != NULL) {
	   	cout << l->p << endl;
	   	l = l->next;
	}
   	l = m_start->next;
   	while (l->next != NULL) {
	   	cout << l->p << endl;
	   	l = l->next;
	}
	*/
   
    l = m_start->next;
    while (l->next != NULL) {
	   	*pout << "    \"" << l->p << "\"," << endl;
	   	l = l->next;
	}
   	*pout << "    NULL" << endl;
   
 	return EXIT_SUCCESS;  
}
   
   
