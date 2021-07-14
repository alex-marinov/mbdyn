/* 
 * File:   MDUtility.cc
 * Author: tingnan
 * 
 * Created on February 8, 2011, 10:24 PM
 */

#include "md_utility.h"

std::string ftoa (double num)   //double to string
{
    std::stringstream converter;
    converter << num;
    return converter.str ();

}

std::string int2str (int num)   //integer to string
{
    std::stringstream converter;
    converter << num;
    return converter.str ();
}

int
sign (double x)
{
    if (x > 0)
        return 1;
    else
        return -1;
}
