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

#include <cstdlib>
#include "unistd.h"
#include <cstdio>
#include <cerrno>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <iomanip>
#include "Leap.h"

int
main(int argc, char *argv[])
{

	long int dt = -1;
	long max_retries = -1;
	long sleep_time = 1;

	struct timespec t0, t;

	while (1) {
		int opt = getopt(argc, argv, "R:S:T:");
		if (opt == -1) {
			break;
		}

		switch (opt) {
		case 'R': {
			char *next = 0;
			max_retries = strtol(optarg, &next, 10);
			if (next == optarg || max_retries <= 0) {
				/* error */
				exit(EXIT_FAILURE);
			}
			} break;

		case 'S': {
			char *next = 0;
			sleep_time = strtol(optarg, &next, 10);
			if (next == optarg || sleep_time <= 0) {
				/* error */
				exit(EXIT_FAILURE);
			}
			} break;

		case 'T': {
			char *next = 0;
			dt = strtol(optarg, &next, 10);
			if (next == optarg || dt <= 0) {
				/* error */
				exit(EXIT_FAILURE);
			}
			} break;

		default:
			exit(EXIT_FAILURE);
		}
	}

	if (dt <= 0) {
		fprintf(stderr, "dt must be defined and positive\n");
		exit(EXIT_FAILURE);
	}

	if (optind != argc) {
		// TODO: error
	}

	// connect to LM
	Leap::Controller controller;
	for (;;) {
		if (controller.isConnected()) {
			break;
		}
		sleep(sleep_time);
		if (max_retries >= 0) {
			if (max_retries == 0) {
				std::cerr << "controller not available" << std::endl;
				exit(EXIT_FAILURE);
			}
			max_retries--;
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &t0);
	t = t0;
	for (;;) {
		int c;

		t.tv_nsec += dt;
		if (t.tv_nsec >= 1000000000L) {
			t.tv_nsec -= 1000000000L;
			t.tv_sec++;
		}

		c = clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &t, NULL);

		Leap::Frame frame = controller.frame();

		Leap::HandList hands = frame.hands();
		Leap::Hand hand = hands.rightmost();
		if (hand.isValid()) {
			const char *what_hand = hand.isRight() ? "right" : "left";
			Leap::Vector pos(hand.palmPosition());
			Leap::Vector nor(hand.palmNormal());
			Leap::Vector dir(hand.direction());
			std::cout
				<< t.tv_sec << "." << std::setfill('0') << std::setw(9) << t.tv_nsec
				<< " " << what_hand << " hand,"
				" position={" << pos.x << ", " << pos.y << ", " << pos.z << "} mm,"
				" direction={" << dir.x << ", " << dir.y << ", " << dir.z << "} mm,"
				" normal={" << nor.x << ", " << nor.y << ", " << nor.z << "} mm,"
				" normal.roll=" << nor.roll() << " rad,"
				" direction.pitch=" << dir.pitch() << " rad,"
				" direction.yaw=" << dir.yaw() << " rad"
				<< std::endl;

		} else {
			std::cout
				<< t.tv_sec << "." << std::setfill('0') << std::setw(9) << t.tv_nsec
				<< " invalid hand" << std::endl;
		}
	}

	return 0;
}
