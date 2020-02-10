
#include <boost/interprocess/managed_shared_memory.hpp>
//#include <boost/interprocess/containers/vector.hpp>
//#include <boost/interprocess/allocators/allocator.hpp>
#include <cstdlib> //std::system
#include <iostream>
#include <string>

#ifndef HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP
#define HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP 1
#endif // HAVE_BOOST_INTERPROCESS_MANAGED_SHARED_MEMORY_HPP

#include "sharedmem.h"

using namespace mbdyn;
using namespace boost::interprocess;

//Main function. For parent process argc == 1, for child process argc == 2
int main(int argc, char *argv[])
{

    std::string shmName ("test_sharedmemory_h");

    int nnodes = 3;
    uint8_t maxi = 10;
    unsigned i = 0;

    if(argc == 1){ // Server process
        std::cout << "[SERVER] Beginning test program for mbdyn shared_memory_buffer" << std::endl;

        // Erase previous shared memory if it exists
        std::cout << "[SERVER] Erase previous shared memory in '" << shmName << "' if it exists" << std::endl;
        shared_memory_object::remove(shmName.c_str ());

        unsigned totbytes = 2 * mapped_region::get_page_size();

        std::cout << "[SERVER] Creating shm with memory size " << totbytes << " bytes" << std::endl;

        managed_shared_memory shm(create_only, shmName.c_str (), totbytes); // 8kb should be sufficient

        // there should only be one buffer in this shared memory region
        std::cout << "[SERVER] Create buffer" << std::endl;
        mbdyn::shared_memory_buffer* s_buffer = shm.construct<shared_memory_buffer>(unique_instance)(nnodes, MBC_ROT_NONE, shmName);

        std::cout << "[SERVER] Getting pointer to X_Data vector" << std::endl;
        shmVector_double* X_data = shm.find<shmVector_double>("X_data_vector").first;

        {
            scoped_lock<interprocess_mutex> lock(s_buffer->mutex);

            s_buffer->peer_ready_for_cmd = false;
            s_buffer->mbdyn_ready_for_cmd = false;
            s_buffer->mbdyn_cmd_available = false;
            s_buffer->peer_cmd_available = false;

        }


        // start the client process
        std::cout << "[SERVER] start the client process" << std::endl;
        std::string s(argv[0]); s += " child &";

        std::cout << "[SERVER] starting client process with command: " << s << std::endl;
        if(0 != std::system(s.c_str()))
            return 1;

        int secs = 2;

        std::cout << "[SERVER] sleeping " << secs << "seconds" << std::endl;
#ifdef _WIN32
        Sleep(secs);
#else
        sleep(secs);
#endif // _WIN32

        uint8_t c_cmd = 0;
        std::cout << "[SERVER] Starting communication using mbdyn shared_memory_buffer. " << std::endl;
        for (uint8_t i = 1; i < maxi+1; i++)
        {
            std::cout << "[SERVER] i is " << unsigned(i) << std::endl;

            {
                scoped_lock<interprocess_mutex> lock(s_buffer->mutex);

                std::cout << "[SERVER] checking if peer is ready to receive" << std::endl;
                if (!s_buffer->peer_ready_for_cmd)
                {
                    std::cout << "[SERVER] waiting for peer be ready to receive" << std::endl;
                    s_buffer->cond_peer_ready_for_cmd.wait (lock);
                }

                // send command to client
                s_buffer->cmd = i;

                std::cout << "[SERVER] put cmd " << unsigned(s_buffer->cmd) << " in buffer" << std::endl;

                std::cout << "[SERVER] notifying new command is available" << std::endl;
                s_buffer->cond_mbdyn_cmd_available.notify_all();

                s_buffer->mbdyn_cmd_available = true;

            }


            {
                scoped_lock<interprocess_mutex> lock(s_buffer->mutex);

                s_buffer->mbdyn_ready_for_cmd = true;

                std::cout << "[SERVER] notifying that server is ready to receive" << std::endl;
                s_buffer->cond_mbdyn_ready_for_cmd.notify_all();

                std::cout << "[SERVER] checking if peer cmd is available" << std::endl;
                if (!s_buffer->peer_cmd_available)
                {
                    std::cout << "[SERVER] waiting for command from client" << std::endl;
                    s_buffer->cond_peer_cmd_available.wait (lock);
                }

                std::cout << "[SERVER] getting command" << std::endl;
                // first value should be 0
                c_cmd = s_buffer->cmd;

                std::cout << "[SERVER] received cmd " << unsigned(c_cmd) << " from peer" << std::endl;

                s_buffer->mbdyn_ready_for_cmd = false;
                // no peer command is available (we've read it)
                s_buffer->peer_cmd_available = false;
            }

        }

        scoped_lock<interprocess_mutex> lock(s_buffer->mutex);

        std::cout << "[SERVER] waiting for client to release semaphore" << std::endl;
        // now test the vectors
        if (!s_buffer->peer_ready_for_cmd)
        {
            s_buffer->cond_peer_ready_for_cmd.wait (lock);
        }

        std::cout << "[SERVER] got semaphore" << std::endl;
        std::cout << "[SERVER] node x data vector contains: ";

        for (i = 0; i < X_data->size(); i++)
            std::cout << X_data->at(i) << ' ';

        std::cout << std::endl;

        std::cout << "[SERVER] putting node x data in shared memory vector: "  << std::endl;

        X_data->at(0) = 1.0;
        X_data->at(1) = 2.0;
        X_data->at(2) = 3.0;

        std::cout << "[SERVER] node x data now contains: ";

        for (i = 0; i < X_data->size(); i++)
            std::cout << X_data->at(i) << ' ';

        std::cout << std::endl;

        s_buffer->peer_cmd_available = true;

        s_buffer->cond_peer_cmd_available.notify_all ();
        std::cout << "[SERVER] released mutex" << std::endl;

    }
    else{ // Client process

        std::cout << "[CLIENT] process started" << std::endl;
        // Open the managed segment
        managed_shared_memory shm (open_only, shmName.c_str ());

        // there should only be one buffer in this shared memory region
        std::pair<shared_memory_buffer*, int> temp = shm.find<shared_memory_buffer>(unique_instance);

        mbdyn::shared_memory_buffer* c_buffer = temp.first;

        std::cout << "[CLIENT] Getting pointer to X_Data vector" << std::endl;
        shmVector_double* X_data = shm.find<shmVector_double>("X_data_vector").first;

        uint8_t s_cmd = 0;

        while (s_cmd < maxi)
        {

            {
                scoped_lock<interprocess_mutex> lock(c_buffer->mutex);

                c_buffer->peer_ready_for_cmd = true;

                std::cout << "[CLIENT] notifying that peer is ready to receive" << std::endl;
                c_buffer->cond_peer_ready_for_cmd.notify_all();


                std::cout << "[CLIENT] checking if mbdyn command is available" << std::endl;
                if (!c_buffer->mbdyn_cmd_available)
                {
                    std::cout << "[CLIENT] waiting for mbdyn command to become available" << std::endl;
                    c_buffer->cond_mbdyn_cmd_available.wait (lock);
                }

                // get the command from the server
                s_cmd = c_buffer->cmd;
                std::cout << "[CLIENT] received cmd " << unsigned(s_cmd) << " from server" << std::endl;

                c_buffer->peer_ready_for_cmd = false;
                c_buffer->mbdyn_cmd_available = false;
            }

            {
                scoped_lock<interprocess_mutex> lock(c_buffer->mutex);

                if (!c_buffer->mbdyn_ready_for_cmd)
                {
                    std::cout << "[CLIENT] waiting for mbdyn to be ready to receive command" << std::endl;
                    c_buffer->cond_mbdyn_ready_for_cmd.wait (lock);
                }

                std::cout << "[CLIENT] sending 10 * " << unsigned(s_cmd) << " to server" << std::endl;
                // multiply by 10 and return
                c_buffer->cmd = 10 * s_cmd;

                c_buffer->peer_cmd_available = true;

                std::cout << "[CLIENT] notifying peer command is available" << std::endl;
                c_buffer->cond_peer_cmd_available.notify_all();
            }

        }

        std::cout << "[CLIENT] waiting to receive node x data from server" << std::endl;
        // now test the vectors
        scoped_lock<interprocess_mutex> lock(c_buffer->mutex);

        c_buffer->peer_ready_for_cmd = true;

        c_buffer->cond_peer_ready_for_cmd.notify_all();

        if (!c_buffer->mbdyn_cmd_available)
        {
            c_buffer->cond_mbdyn_cmd_available.wait (lock);
        }

        std::cout << "[CLIENT] received node x data " << std::endl;

        for (i = 0; i < X_data->size(); i++)
            std::cout << X_data->at(i) << ' ';

        std::cout << std::endl;

    }

   return 0;
}
