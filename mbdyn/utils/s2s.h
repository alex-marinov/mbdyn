/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef S2S_H
#define S2S_H

const int	S2S_BUFSIZE = 1024;

struct s2s_sockaddr_t {
	int	ms_type;
	int	ms_domain;
	union {
		struct sockaddr		ms_addr_generic;
		struct sockaddr_in	ms_addr_inet;
		struct sockaddr_un	ms_addr_local;
	}	ms_addr;
	socklen_t	ms_len;
};

enum {
	BLOCK_UNKNOWN = 0,
	BLOCK_NO,
	BLOCK_YES
};

struct s2s_t {
	// depend on configuration
	int	nChannels;
	const char	*path;
	// path and [host:][port] are mutually exclusive
	const char	*host;
	int	port;
	s2s_sockaddr_t	addr;
	int	create;
	int	block;

	// internal data
	int	sock;

	std::string buf;
	std::vector<double>	dbuf;
	bool	stream2socket;
	const char	*progname;

	s2s_t(void);
	~s2s_t(void);

	void usage(int rc) const;
	void parse(int argc, char *argv[]);
	void prepare(void);
	void shutdown(void);

	bool is_blocking(void) const;

	ssize_t send(int flags) const;
	ssize_t recv(int flags);

	struct s2s_t	*next;
};

#endif /* S2S_H */
