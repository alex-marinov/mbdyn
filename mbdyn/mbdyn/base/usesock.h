/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef USESOCK_H
#define USESOCK_H
#ifdef USE_SOCKET
#include <sys/socket.h>
#include <netinet/in.h>
#include <sys/un.h>

class UseSocket {
protected:
	int sock;
	int socket_type;
	bool create;
	bool connected;
	bool abandoned;

	mutable socklen_t socklen;

	void PostConnect(void);

public:
	UseSocket(bool c);
	UseSocket(int t, bool c);
	virtual ~UseSocket(void);
	
	virtual std::ostream& Restart(std::ostream& out) const;

	int GetSock(void) const;
	void SetSock(int s);
	virtual void Connect(void);
	virtual void ConnectSock(int s);
	bool Create(void) const;
	bool Connected(void) const;
	void Abandon(void);
	bool Abandoned(void) const;

	socklen_t& GetSocklen(void) const;
	virtual struct sockaddr *GetSockaddr(void) const = 0;

	ssize_t send(const void *buf, size_t len, int flags);
	ssize_t recv(void *buf, size_t len, int flags);
};

class UseInetSocket : public UseSocket {
protected:
	std::string host;
 	unsigned short int port;
	struct sockaddr_in addr;

	void UseInetSocket_int(void);

public:
	UseInetSocket(const std::string& h, unsigned short p, bool c);
	UseInetSocket(const std::string& h, unsigned short p, int t, bool c);
	virtual ~UseInetSocket(void);
	
	std::ostream& Restart(std::ostream& out) const;

	void Connect(void);
	void ConnectSock(int s);
	struct sockaddr *GetSockaddr(void) const;
};

class UseLocalSocket : public UseSocket {
protected:
	std::string path;
	struct sockaddr_un addr;

	void UseLocalSocket_int(void);

public:
	UseLocalSocket(const std::string& p, bool c);
	UseLocalSocket(const std::string& p, int t, bool c);
	virtual ~UseLocalSocket(void);
	
	std::ostream& Restart(std::ostream& out) const;

	void Connect(void);
	void ConnectSock(int s);
	struct sockaddr *GetSockaddr(void) const;
};
#endif
#endif /* USESOCK_H */

