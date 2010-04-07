%module mbc_py
%include "typemaps.i"
%header
%{
        #include <numpy/arrayobject.h>
        #include "mbc_py_global.h"

        extern int
        mbc_py_nodal_initialize(const char *const path,
                const char *const host, unsigned port,
                unsigned data_and_next, unsigned verbose,
                unsigned rigid, unsigned nodes,
                unsigned labels, unsigned rot, unsigned accels);

        extern int
        mbc_py_nodal_send(int last);

        extern int
        mbc_py_nodal_recv(void);

        extern void
        mbc_py_nodal_destroy(void);
%}

%include "mbc_py_global.h"
int
mbc_py_nodal_initialize(const char *path,
        const char *host, unsigned port,
        unsigned data_and_next, unsigned verbose,
        unsigned rigid, unsigned nodes,
        unsigned labels, unsigned rot, unsigned accels);

int
mbc_py_nodal_send(int last);

int
mbc_py_nodal_recv(void);

void
mbc_py_nodal_destroy(void);

%init
%{
        import_array();
%}

