/* 
 * File:   MDsand.h
 * Author: tingnan
 *
 * Created on January 14, 2011, 1:42 PM
 */
#include "md_header.h"
#include <vector>
#include <iostream>
#include <fstream>

#ifndef _MDSAND_H
#define	_MDSAND_H


#if MD_DIM == 2
#define MD_OFFSET_VALS                                          \
   {{0,0}, {1,0}, {1,1}, {0,1}, {-1,1}}
#define MD_N_NBR  5
#endif

#define MD_OFFSET_VALS                                          \
   { {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},              \
     {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,1,1}, {-1,0,1},    \
     {-1,-1,1}, {0,-1,1}, {1,-1,1}                              \
   }
#define MD_N_NBR 14

struct MDDisk;
struct MDSlip;

struct MDShapeParser
{
    class MDShape **pshape;
    int n;
};

//the cell, a chain structure, ref_num stores which ball is in this cell

struct MDCell
{
    int ref_num;
    MDCell *next;
};

class MDRegion
{
private:
    int n_proc_;
    int n_subregion_;
    int n_par_;
    int regionx_;
    int regiony_;
    int regionz_;
    MDDisk *sand_pool_;
    std::vector <class MDSand> subportion_;
public:
    MDRegion(int, int, int, int, int);
    void RegionSplit();
    void DataExchange();
};

class MDSand
{
private:
    int n_boulder_;
    MDDisk *sand_array_;
    MDSlip **slip_array_;
    MDCell **box_cell_;
    double sys_energy_;
    double sys_kinematic_energy_;
    double sys_potential_Energy_;
    double sys_dissip_energy_;
    //--------------------
    double box_length_x_;
    double box_length_z_;
    double box_length_y_;
    double cell_length_;
    //--------------------
    double radius_min_;
    double radius_max_;
    double mass_min_;
    double mass_max_;
    double rho_;
    double bip_ratio_;
    double kn_;
    double kt_;
    double gn_;
    double gs_;
    double mu_pp_;
    double mu_bp_;
    //---------------------
    double step_tau_;
    double sys_time_;
    double sav_intval_;
    int sys_step_;
    int sav_intval_step_;
    //---------------------
    bool if_output_;
    int boundary_para_;
    int max_contact_n_;
    int num_particles_;
    int cell_num_x_;
    int cell_num_y_;
    int cell_num_z_;
    std::ofstream sand_file_;
    std::ofstream para_file_;
    std::ofstream energy_file_;
    class IntVec * md_cell_nbr_;
public:
    MDSand(const MDSand&);
    MDSand(double, double, double, int);
    ~MDSand();
    MDSand& operator=(const MDSand&);
    void Init();
    void Interact(int, int);
    void Run(int, const MDShapeParser &, bool);
    int TransCell(int, const class IntVec &, const class IntVec &);
    MDCell *RemovefromCell(int, const class IntVec &);
    void InitCell(int, const class IntVec &);
    void LoadSand();
    void LoadSand(const char *);
    void InteractVirtualBall(int, class RealVec &, class RealVec &, class RealVec &, double,
                             double);
    void InteractExt(const MDShapeParser &);
    void PositionToCell(class RealVec &, class IntVec &);
    void GenerateSand();
    void ExForce();
    void RunOneStep();
    void SwitchOutput(bool);
    void OutputFile();
    void OutputState();
    void BoundaryCondition(unsigned int);
    void Integration();
    bool RegionOut(const class IntVec &);
    bool RegionEnd(const class IntVec &);
    int  CellIndexToDimOne(const class IntVec &);
    void CellIndexToDim23(int, class IntVec &);
    void InterateCellIndex(class IntVec &);
};

#endif
