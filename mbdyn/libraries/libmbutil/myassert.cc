/******************************************************************************

Macro di assert personalizzate

Uso: ASSERT(<expr>);
	- se <expr> e' vera ( != 0 ) non fa nulla;
	- se <expr> e' falsa, scrive sul flusso di errore cerr il file e la riga
		solo se DEBUG e' definita.

Uso: ASSERTMSG(<expr>, <msg>);
	- se <expr> e' vera ( != 0 ) non fa nulla;
	- se <expr> e' falsa, scrive sul flusso di errore cerr il file e la riga,
		seguiti dal messaggio <msg>, solo se DEBUG e' definita.

Entrambe chiamano la funzione _Assert(file, line, msg = NULL);
se msg e' definito, viene aggiunto in coda al messaggio di default

******************************************************************************/

#include <mbconfig.h>

#include <myassert.h>
#include <string.h>

#ifdef DEBUG

long int debug_level = MYDEBUG_ANY;
long int DEFAULT_DEBUG_LEVEL = MYDEBUG_ANY;

void _Assert(const char* file, const int line, const char* msg)
{
   cout.flush();
   
   cerr << endl << "ASSERT fault in file " << file 
     << " at line " << line;
   if (msg) { 
      cerr << ':' << endl << msg; 
   }
   cerr << endl;
   
#ifdef DEBUG_STOP
   THROW(MyAssert::ErrGeneric());
#endif   
   
   return;
}

ostream& _Out(ostream& out, const char* file, const int line)
{
   cout.flush();
   
   // out << "File <" << file << ">, line [" << line << "]: ";
   out << "[" << file << "," << line << "]: ";
   return out;
}

int get_debug_options(const char *const s, const debug_array da[])
{
   if (s == NULL || s[0] == '\0') {
      ::debug_level = DEFAULT_DEBUG_LEVEL;
      return 0;
   }
   
   char* p = (char*)s;
   while (1) {
      char* sep = strchr(p, ':');
      unsigned int l;
      if (sep != NULL) {
	 l = int(sep-p);
      } else {
	 l = strlen(p);
      }
      debug_array* w = (debug_array*)da;
      while (w->s != NULL) {
	 if (l == strlen(w->s) && strncmp(w->s, p, l) == 0) {
	    break;
	 }
	 w++;
      }
      if (w->s == NULL) {
	 if (l == 4 && strncmp("none", p, 4) == 0) {
	    ::debug_level = MYDEBUG_NONE;
	 } else if (l == 3 && strncmp("any", p, 3) == 0) {
	    ::debug_level = MYDEBUG_ANY;
	 } else {
	    cerr << "Unknown debug level \"";
	    for (unsigned int i = 0; i < l; i++) {
	       cerr << p[i];
	    }
	    cerr << "\"" << endl;
	 }
      } else {
	 ::debug_level |= w->l;
	 cerr << "debug level: " << w->s << endl;
      }
      if (sep == NULL) {
	 break;
      }
      p = sep+1;
   }
   
   return 0;
}

#endif /* DEBUG */

