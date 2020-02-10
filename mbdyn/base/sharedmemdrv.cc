/* $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/mbdyn/base/sockdrv.cc,v 1.48 2017/01/12 14:46:10 masarati Exp $ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati    <masarati@aero.polimi.it>
 * Paolo Mantegazza    <mantegazza@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "dataman.h"
#include "sharedmemdrv.h"

#ifdef USE_BOOST

#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>

#include "sharedmem.h"

using namespace boost::interprocess;

const size_t USERLEN = 32;
const size_t CREDLEN = 128;
const size_t BUFSIZE = 1024;
const char *MBDynSharedMemName = "mbdyn_shared_mem";


SharedMemDrive::SharedMemDrive(unsigned int uL, const DriveHandler* pDH,
                               const char *shm_name,
                               integer nd, const std::vector<doublereal>& v0)
    : FileDrive(uL, pDH, "socket", nd, v0),
      auth(NULL),
      pFlags(NULL)
{
    int            save_errno;

    ASSERT(path != NULL);
    ASSERT(nd > 0);

    SAFENEW(auth, NoAuth);

    /* Create the socket and set it up to accept connections. */
    SAFESTRDUP(Name, shm_name);
//       shm = mbdyn_make_named_shm (0, Name, 1, &save_errno);
//       if (shm == -1) {
//        const char    *err_msg = strerror(save_errno);
//
//              silent_cerr("SharedMemDrive(" << GetLabel()
//            << "): shared memory creation failed "
//            "(" << save_errno << ": "<< err_msg << ")"
//            << std::endl);
//              throw ErrGeneric(MBDYN_EXCEPT_ARGS);
//
//       }

    Init();
}

void
SharedMemDrive::Init(void)
{
    //use old shared memory if exists else create a new one
    shm = shared_memory_object (open_or_create, Name, read_write);

    //set the size of the memory object
    shm.truncate(sizeof(shared_memory_buffer));

    //map the whole shared memory in this process
    region = mapped_region (shm,read_write);

    //get the region address
    addr = region.get_address();

    //create a shared memory buffer in memory
    data = shared_memory_buffer (addr);

    SAFENEWARR(pFlags, int, iNumDrives + 1);
    for (int iCnt = 0; iCnt <= iNumDrives; iCnt++)
    {
        pFlags[iCnt] = SharedMemDrive::DEFAULT;
    }
}


SharedMemDrive::~SharedMemDrive(void)
{
    /* some shutdown stuff ... */
    if (auth != NULL)
    {
        SAFEDELETE(auth);
    }
}

static char *
get_line(char *buf, size_t bufsize,  FILE *fd)
{
    int len;

    if (fgets(buf, bufsize, fd) == NULL)
    {
        return NULL;
    }

    len = strlen(buf);
    if (len > 0 && buf[len-1] == '\n')
    {
        buf[len-1] = '\0';
        if (len > 1 && buf[len-2] == '\r')
        {
            buf[len-2] = '\0';
        }
    }
    else
    {
        fprintf(stderr, "buffer overflow\n");
        return NULL;
    }

    return buf;
}

/*
 * in/out:
 *     user : puntatore a buffer di 9 bytes (8+'\0')
 *     cred : puntatore a buffer di 129 bytes (128+'\0')
 *
 * out:
 *     buf : puntatore a buffer statico che contiene una copia della nuova
 *           linea nel caso sia stata accidentalmente letta
 */
int
get_auth_token(FILE *fd, char *user, char *cred, char **nextline)
{
    char buf[BUFSIZE];

    if (get_line(buf, BUFSIZE, fd) == NULL)
    {
        return -1;
    }

    user[0] = '\0';
    cred[0] = '\0';

    if (strncasecmp(buf, "user:", 5) != 0)
    {
        *nextline = (char *)malloc(sizeof(char)*(strlen(buf)+1));
        if ((*nextline) == NULL)
        {
            return -1;
        }
        strcpy(*nextline, buf);
        return 0;
    }
    else
    {
        char *p;
        unsigned int i;

        p = buf+5;
        while (isspace(*p))
        {
            p++;
        }
        for (i = 0; i < USERLEN; i++)
        {
            if ((user[i] = p[i]) == '\0' || isspace(user[i]))
            {
                break;
            }
        }
        user[i] = '\0';
    }

    if (get_line(buf, BUFSIZE, fd) == NULL)
    {
        return -1;
    }

    if (strncasecmp(buf, "password:", 9) != 0)
    {
        *nextline = (char *)malloc(sizeof(char)*(strlen(buf)+1));
        if ((*nextline) == NULL)
        {
            return -1;
        }
        strcpy(*nextline, buf);
        return 0;
    }
    else
    {
        char *p;
        unsigned int i;

        p = buf+9;
        while (isspace(*p))
        {
            p++;
        }
        for (i = 0; i < CREDLEN; i++)
        {
            if ((cred[i] = p[i]) == '\0' || isspace(cred[i]))
            {
                break;
            }
        }
        cred[i] = '\0';
    }

    return 1;
}

void
SharedMemDrive::ServePending(const doublereal& /* t */ )
{
    /* prova */
    for (integer iCnt = 1; iCnt <= iNumDrives; iCnt++)
    {
        if (pFlags[iCnt] & SharedMemDrive::IMPULSIVE)
        {
            pdVal[iCnt] = 0.;
        }
    }

    while (true)
    {
        int got_value = 0;
        char *nextline = NULL;
        const size_t bufsize = BUFSIZE;
        char buf[bufsize];

        int label;
        doublereal value;

        bool bAuthc = false;
#ifdef HAVE_SASL2
        if (dynamic_cast<SASL2_Auth*>(auth))
        {
            if (auth->Auth(cur_sock) != AuthMethod::AUTH_OK)
            {
                silent_cerr(
                    "SharedMemDrive(" << GetLabel() << "): "
                    "authentication failed" << std::endl);
                continue;
            }
            bAuthc = true;
        }
#endif /* HAVE_SASL2 */

        FILE* fd = fdopen(cur_sock, "r");

        if (!bAuthc)
        {
            char user[USERLEN + 1];
            char cred[CREDLEN + 1];

            if (get_auth_token(fd, user, cred, &nextline) == -1)
            {
                silent_cerr(
                    "SharedMemDrive(" << GetLabel() << "): "
                    "corrupted stream" << std::endl);
                fclose(fd);
                continue;
            }

            DEBUGCOUT("got auth token: user=\"" << user
                      << "\", cred=\"" << cred << "\"" << std::endl);

            if (auth->Auth(user, cred) != AuthMethod::AUTH_OK)
            {
                silent_cerr(
                    "SharedMemDrive(" << GetLabel() << "): "
                    "authentication failed" << std::endl);
                fclose(fd);
                continue;
            }
        }

        DEBUGCOUT("authenticated" << std::endl);

        /*
         * la nuova linea puo' essere gia' stata letta
         * da get_auth_token
         */
        if (nextline == NULL)
        {
            if (get_line(buf, bufsize, fd) == NULL)
            {
                silent_cerr(
                    "SharedMemDrive(" << GetLabel() << "): "
                    "corrupted stream" << std::endl);
                fclose(fd);
                continue;
            }
        }
        else
        {
            strncpy(buf, nextline, bufsize);
            free(nextline);
        }
        nextline = buf;

        /* read the label */
        if (strncasecmp(nextline, "label:", 6) != 0)
        {
            silent_cerr("SharedMemDrive(" << GetLabel() << "): "
                        "missing label" << std::endl);
            fclose(fd);
            continue;
        }

        char *p = nextline + 6;
        while (isspace(p[0]))
        {
            p++;
        }

        if (sscanf(p, "%d", &label) != 1)
        {
            silent_cerr("SharedMemDrive(" << GetLabel() << "): "
                        "unable to read label" << std::endl);
            fclose(fd);
            continue;
        }

        if (label <= 0 || label > iNumDrives)
        {
            silent_cerr("SharedMemDrive(" << GetLabel() << "): "
                        "illegal label " << label << std::endl);
            fclose(fd);
            continue;
        }

        while (true)
        {
            if (get_line(buf, bufsize, fd) == NULL)
            {
                silent_cerr(
                    "SharedMemDrive(" << GetLabel() << "): "
                    "corrupted stream" << std::endl);
                fclose(fd);
                break;
            }

            nextline = buf;

            if (nextline[0] == '.')
            {
                fclose(fd);
                break;
            }

            if (strncasecmp(nextline, "value:", 6) == 0)
            {
                char *p = nextline+6;
                while (isspace(p[0]))
                {
                    p++;
                }

                if (sscanf(p, "%lf", &value) != 1)
                {
                    silent_cerr("SharedMemDrive(" << GetLabel() << "): "
                                "unable to read value"
                                << std::endl);
                    fclose(fd);
                    break;
                }
                got_value = 1;

            }
            else if (strncasecmp(nextline, "inc:", 4) == 0)
            {
                char *p = nextline+4;
                while (isspace(p[0]))
                {
                    p++;
                }

                if (strncasecmp(p, "yes", 3) == 0)
                {
                    pFlags[label] |= SharedMemDrive::INCREMENTAL;
                }
                else if (strncasecmp(p, "no", 2) == 0)
                {
                    pFlags[label] &= !SharedMemDrive::INCREMENTAL;
                }
                else
                {
                    silent_cerr("SharedMemDrive(" << GetLabel() << "): "
                                "\"inc\" line in "
                                "\"" << nextline << "\" "
                                "looks corrupted"
                                << std::endl);
                    fclose(fd);
                    break;
                }
                nextline = NULL;

            }
            else if (strncasecmp(nextline, "imp:", 4) == 0)
            {
                char *p = nextline+4;
                while (isspace(p[0]))
                {
                    p++;
                }

                if (strncasecmp(p, "yes", 3) == 0)
                {
                    pFlags[label] |= SharedMemDrive::IMPULSIVE;
                }
                else if (strncasecmp(p, "no", 2) == 0)
                {
                    pFlags[label] &= !SharedMemDrive::IMPULSIVE;
                }
                else
                {
                    silent_cerr("SharedMemDrive(" << GetLabel() << "): "
                                "\"imp\" line" " in "
                                "\"" << nextline << "\""
                                " looks corrupted"
                                << std::endl);
                    fclose(fd);
                    break;
                }
                nextline = NULL;
            }

            /* usa i valori */
            if (got_value)
            {
                if (pFlags[label] & SharedMemDrive::INCREMENTAL)
                {
                    silent_cout("SharedMemDrive(" << GetLabel() << "): "
                                "adding " << value
                                << " to label " << label
                                << std::endl);
                    pdVal[label] += value;

                }
                else
                {
                    silent_cout("SharedMemDrive(" << GetLabel() << "): "
                                "setting label " << label
                                << " to value " << value
                                << std::endl);
                    pdVal[label] = value;
                }
            }
        }
    }
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
SharedMemDrive::Restart(std::ostream& out) const
{
    return out << "SharedMemDrive not implemented yet" << std::endl;
}

/* legge i drivers tipo socket */

Drive *
SharedMemDR::Read(unsigned uLabel, const DataManager *pDM, MBDynParser& HP)
{
    Drive* pDr = NULL;

    integer idrives = HP.GetInt();
    unsigned short int port = MBDynSharedMemDrivePort;
    const char *name = NULL;

    std::vector<doublereal> v0;
    if (HP.IsKeyWord("initial" "values"))
    {
        v0.resize(idrives);
        for (integer i = 0; i < idrives; i++)
        {
            v0[i] = HP.GetReal();
        }
    }

    name = HP.GetString();
    ASSERT(name != NULL);

    SAFENEWWITHCONSTRUCTOR(pDr,
                           SharedMemDrive,
                           SharedMemDrive(uLabel, pDM->pGetDrvHdl(),
                                          path, idrives, v0));

    return pDr;
}

#endif /* USE_BOOST */
