// The units in this program are second, centimeter, gram.
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <string>
#include <sstream>
#include <cassert>
#include "md_sand.h"
#include "md_shape.h"
#include "md_utility.h"
#include "md_vec.h"

using namespace std;
const double math_pi = 3.1415926535897932384626433832795;
const int nnbr[][MD_DIM] = MD_OFFSET_VALS;
const RealVec gravity_vector(0, 0, -double (GRAVITY));

struct MDDisk
{
    IntVec index;
    RealVec pos;
    RealVec v;
    RealVec omega;
    RealVec moment;
    RealVec F;
    double r;
    double m;
    double abF;
    bool remove;

    MDDisk() : remove(0)
    {
    };
};

struct MDSlip
{
    int ref_num;
    RealVec disp;
    RealVec ndir;
    MDSlip *next;

    MDSlip() : next(NULL)
    {
    };
};

MDRegion::MDRegion(int particles, int processors, int nx, int ny, int nz) :
n_proc_(processors),
n_subregion_(processors),
n_par_(particles),
regionx_(nx),
regiony_(ny),
regionz_(nz)
{
    sand_pool_ = new MDDisk [particles];
    MDSand temp(regionx_, regiony_, regionz_, n_par_);
    subportion_.resize(n_proc_, temp);
}

void
MDRegion::RegionSplit()
{
    int npar_per_region = n_par_ / n_subregion_;
}

MDSand::MDSand(const MDSand &T)
{
    radius_min_ = T.radius_min_;
    radius_max_ = T.radius_max_;
    mass_min_ = T.mass_min_;
    mass_max_ = T.mass_max_;
    bip_ratio_ = T.bip_ratio_;
    kn_ = T.kn_;
    kt_ = T.kt_;
    gn_ = T.gn_;
    gs_ = T.gs_;
    mu_pp_ = T.mu_pp_;
    mu_bp_ = T.mu_bp_;
    cell_length_ = T.cell_length_;
    cell_num_x_ = T.cell_num_x_;
    cell_num_y_ = T.cell_num_y_;
    cell_num_z_ = T.cell_num_z_;
    sys_step_ = T.sys_step_;
    sys_time_ = T.sys_time_;
    step_tau_ = T.step_tau_;
    sav_intval_ = T.sav_intval_;
    sav_intval_step_ = T.sav_intval_step_;

    if_output_ = T.if_output_;
    boundary_para_ = T.boundary_para_;
    sys_energy_ = T.sys_energy_;
    sys_potential_Energy_ = T.sys_potential_Energy_;
    sys_kinematic_energy_ = T.sys_kinematic_energy_;
    sys_dissip_energy_ = T.sys_dissip_energy_;
    sand_array_ = new MDDisk[num_particles_];
    box_cell_ = new MDCell *[cell_num_x_ * cell_num_y_ * cell_num_z_];
    slip_array_ = new MDSlip *[num_particles_];
    md_cell_nbr_ = new IntVec [MD_N_NBR];
    max_contact_n_ = T.max_contact_n_;
    int i, j;
    for (i = 0; i < MD_N_NBR; i++)
    {
        md_cell_nbr_[i].x = nnbr[i][0];
        md_cell_nbr_[i].y = nnbr[i][1];
        md_cell_nbr_[i].z = nnbr[i][2];
    }
    for (i = 0; i < num_particles_; i++)
    {
        slip_array_[i] = new MDSlip[max_contact_n_];
        for (j = 0; j < max_contact_n_; j++)
        {
            slip_array_[i][j].ref_num = 0;
            slip_array_[i][j].next = NULL;
        }
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        box_cell_[i] = NULL;
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        if (box_cell_[i] != NULL)
        {
            cout << i << "\t" << box_cell_[i] << endl;
            assert(0);
        }
    }
}

MDSand& MDSand::operator=(const MDSand &T)
{
    radius_min_ = T.radius_min_;
    radius_max_ = T.radius_max_;
    mass_min_ = T.mass_min_;
    mass_max_ = T.mass_max_;
    bip_ratio_ = T.bip_ratio_;
    kn_ = T.kn_;
    kt_ = T.kt_;
    gn_ = T.gn_;
    gs_ = T.gs_;
    mu_pp_ = T.mu_pp_;
    mu_bp_ = T.mu_bp_;
    cell_length_ = T.cell_length_;
    cell_num_x_ = T.cell_num_x_;
    cell_num_y_ = T.cell_num_y_;
    cell_num_z_ = T.cell_num_z_;
    sys_step_ = T.sys_step_;
    sys_time_ = T.sys_time_;
    step_tau_ = T.step_tau_;
    sav_intval_ = T.sav_intval_;
    sav_intval_step_ = T.sav_intval_step_;

    if_output_ = T.if_output_;
    boundary_para_ = T.boundary_para_;
    sys_energy_ = T.sys_energy_;
    sys_potential_Energy_ = T.sys_potential_Energy_;
    sys_kinematic_energy_ = T.sys_kinematic_energy_;
    sys_dissip_energy_ = T.sys_dissip_energy_;
    sand_array_ = new MDDisk[num_particles_];
    box_cell_ = new MDCell *[cell_num_x_ * cell_num_y_ * cell_num_z_];
    slip_array_ = new MDSlip *[num_particles_];
    md_cell_nbr_ = new IntVec [MD_N_NBR];
    max_contact_n_ = T.max_contact_n_;
    int i, j;

    for (i = 0; i < MD_N_NBR; i++)
    {
        md_cell_nbr_[i].x = nnbr[i][0];
        md_cell_nbr_[i].y = nnbr[i][1];
        md_cell_nbr_[i].z = nnbr[i][2];
    }
    for (i = 0; i < num_particles_; i++)
    {
        slip_array_[i] = new MDSlip[max_contact_n_];
        for (j = 0; j < max_contact_n_; j++)
        {
            slip_array_[i][j].ref_num = 0;
            slip_array_[i][j].next = NULL;
        }
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        box_cell_[i] = NULL;
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        if (box_cell_[i] != NULL)
        {
            cout << i << "\t" << box_cell_[i] << endl;
            assert(0);
        }
    }
    return *this;
}

MDSand::MDSand(double regionx, double regiony, double regionz, int numberparticles) :
box_length_x_(regionx),
box_length_y_(regiony),
box_length_z_(regionz),
num_particles_(numberparticles)
{
    radius_min_ = (0.321 - 0.02) / 2;
    radius_max_ = (0.321 + 0.02) / 2;
    mass_min_ = 0.03555;
    mass_max_ = 0.05170;
    bip_ratio_ = 0.5;
    kn_ = 2e8;
    kt_ = 2e8;
    gn_ = 1500;
    gs_ = 0;
    mu_pp_ = 0.1;
    mu_bp_ = 0.27;
    // region geometry
    cell_length_ = 2 * radius_max_;
    cell_num_x_ = box_length_x_ / cell_length_ + 2;
    cell_num_y_ = box_length_y_ / cell_length_ + 2;
    cell_num_z_ = box_length_z_ / cell_length_ + 2;
    // initial system time and time step, time step length
    sys_step_ = 0;
    sys_time_ = 0;
    step_tau_ = 2.5e-6;
    sav_intval_ = 1e-3;
    sav_intval_step_ = int(sav_intval_ / step_tau_);

    if_output_ = 1;
    boundary_para_ = 1;
    sys_energy_ = 0;
    sys_potential_Energy_ = 0;
    sys_kinematic_energy_ = 0;
    sys_dissip_energy_ = 0;
    sand_array_ = new MDDisk[num_particles_];
    box_cell_ = new MDCell *[cell_num_x_ * cell_num_y_ * cell_num_z_];
    slip_array_ = new MDSlip *[num_particles_];
    md_cell_nbr_ = new IntVec [MD_N_NBR];
    max_contact_n_ = 15 * 0;
    int i, j;
    for (i = 0; i < MD_N_NBR; i++)
    {
        md_cell_nbr_[i].x = nnbr[i][0];
        md_cell_nbr_[i].y = nnbr[i][1];
        md_cell_nbr_[i].z = nnbr[i][2];
    }
    for (i = 0; i < num_particles_; i++)
    {
        slip_array_[i] = new MDSlip[max_contact_n_];
        for (j = 0; j < max_contact_n_; j++)
        {
            slip_array_[i][j].ref_num = 0;
            slip_array_[i][j].next = NULL;
        }
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        box_cell_[i] = NULL;
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        if (box_cell_[i] != NULL)
        {
            cout << i << "\t" << box_cell_[i] << endl;
        }
    }
    sand_file_.open("san.dat", ios::out | ios::binary);
    if (sand_file_.fail())
    {
        cout << "can't open .san.dat" << endl;
    }
    energy_file_.open("energy.txt", ios::out | ios::binary);
    if (energy_file_.fail())
    {
        cout << "can't open energy.txt" << endl;
    }
}

/*
 * void MDSand::Init(int ibatch) { int i; for (i = ibatch * 100; i < ibatch * 100 + 100; i++) { IntVec vecint;
 * PositionToCell(sand_array_[i].pos, vecint); if (vecint.out()) cout << "wrong at Init() add" << endl; InitCell(i, vecint); } }
 */

bool
MDSand::RegionOut(const class IntVec & vec)
{
    if (vec.x < 0 || vec.y < 0 || vec.z < 0 || vec.x >= cell_num_x_
        || vec.y >= cell_num_y_ || vec.z >= cell_num_z_)
        return 1;
    else
        return 0;
}

bool
MDSand::RegionEnd(const class IntVec & vec)
{
    if (vec.x == cell_num_x_ - 1 && vec.z == cell_num_z_ - 1 && vec.y == cell_num_y_ - 1)
        return 1;
    else
        return 0;
}

int
MDSand::CellIndexToDimOne(const class IntVec & vec)
{
    return vec.x + vec.y * cell_num_x_ + vec.z * cell_num_x_ * cell_num_y_;
}

void
MDSand::CellIndexToDim23(int i, IntVec & vec)
{
    vec.z = i / (cell_num_x_ * cell_num_y_);
    i = i - vec.z * cell_num_x_ * cell_num_y_;
    vec.y = i / cell_num_x_;
    vec.x = i - vec.y * cell_num_x_;
}

void
MDSand::Init()
{

    int i;
    int tempcounter = 0;
    MDCell *temppair;

    for (i = 0; i < num_particles_; i++)
    {
        IntVec vecint;
        PositionToCell(sand_array_[i].pos, vecint);
        if (RegionOut(vecint))
            cout << "wrong at Init() add" << endl;
        InitCell(i, vecint);
    }
    for (i = 0; i < cell_num_x_ * cell_num_z_ * cell_num_y_; i++)
    {
        temppair = box_cell_[i];
        while (temppair != NULL)
        {
            temppair = (*temppair).next;
            tempcounter++;
        }
    }
    assert(tempcounter == num_particles_ && tempcounter);
    cout << "total of " << tempcounter << " sand particles loaded" << endl;
    cout << "current box dimension is " << box_length_x_ << "x" << box_length_y_
        << "x" << box_length_z_ << " cm" << endl;

    for (i = 0; i < num_particles_; i++)
    {
        sys_kinematic_energy_ =
            sys_kinematic_energy_ + 0.5 * sand_array_[i].m * (sand_array_[i].v * sand_array_[i].v);
        sys_potential_Energy_ =
            sys_potential_Energy_ + sand_array_[i].m * GRAVITY * sand_array_[i].pos.z;
    }
}

void
MDSand::PositionToCell(class RealVec & sand_position, class IntVec & cell_position)
{

    cell_position.x = (int) (sand_position.x / cell_length_);
    cell_position.y = (int) (sand_position.y / cell_length_);
    cell_position.z = (int) (sand_position.z / cell_length_);
    if (RegionOut(cell_position))
    {
        cout << "error! out of cell boundary" << endl;
    }
}

int
MDSand::TransCell(int i, const class IntVec & vec1, const class IntVec & vec2)
{
    MDCell *pointa, *pointb;

    // let us keep the error check part.
    if (RegionOut(vec1) || RegionOut(vec2))
    {
        cout << "" << endl;
        cout << vec1.x << "\t" << vec2.x << endl;
        cout << vec1.y << "\t" << vec2.y << endl;
        cout << vec1.z << "\t" << vec2.z << endl;
        exit(0);
    }
    // ------------------------------------------

    pointa = RemovefromCell(i, vec1);
    if (pointa == NULL)
    {
        cout << "error! Cannot find the particle in cell." << endl;
        return -1;
    }
    pointb = box_cell_[CellIndexToDimOne(vec2)];
    if (pointb == NULL)
    {
        box_cell_[CellIndexToDimOne(vec2)] = pointa;
        return 1;
    }
    else
    {
        while ((*pointb).next != NULL)
            pointb = (*pointb).next;
        (*pointb).next = pointa;
    }
    return 0;
}

// remove an atom from a specific cell and return the pointer to the freed cell

struct MDCell *
MDSand::RemovefromCell(int i, const class IntVec & vector)
{
    struct MDCell *p, *pnext;

    // need to make sure this atom i belong to the vector specified by vector
    p = box_cell_[CellIndexToDimOne(vector)];
    assert(p && "cell is empty");
    if ((*p).ref_num == i)
    {
        box_cell_[CellIndexToDimOne(vector)] = (*p).next;
        (*p).next = NULL;
        return p;
    }
    else
    {
        pnext = (*p).next;
        while ((*pnext).ref_num != i && (*pnext).next != NULL)
        {
            p = pnext;
            pnext = (*pnext).next;
        }
        assert((*pnext).ref_num == i && "not found in the cell");
        (*p).next = (*pnext).next;
        (*pnext).next = NULL;
        return pnext;
    }
}

// Add ith ball to the cell it should belong to.
// IntVec should be calculated by PositionToCell to avoid error.

void
MDSand::InitCell(int i, const class IntVec & vec)
{

    struct MDCell *p, *n;

    // I am not sure what is this for atm
    if (RegionOut(vec))
    {
        cout << i << endl;
        cout << sand_array_[i].pos.x << '\t' << sand_array_[i].pos.y << '\t' << sand_array_[i].
            pos.z << endl;
        cout << vec.x << '\t' << vec.y << '\t' << vec.z << endl;
        cout << "error: Cell number is  is two small" << endl;
    }
    // after allocating memory remember to keep all of them under track
    // to avoid memory leak
    n = (struct MDCell *) malloc(sizeof (struct MDCell));
    (*n).ref_num = i;
    (*n).next = NULL;
    p = box_cell_[CellIndexToDimOne(vec)];

    if (p != NULL)
    {
        while ((*p).next != NULL)
            p = (*p).next;
        (*p).next = n;
    }
    else
    {
        box_cell_[CellIndexToDimOne(vec)] = n;
    }

}

// load previously saved sand into the system. The counterpart is generate new
// sand

void
MDSand::LoadSand(const char *filename)
{
    ifstream sand;
    sand.open(filename, ios::in | ios::binary);
    if (sand.fail())
    {
        cout << "Can't open " << filename << endl;
        exit(-1);
    }
    sand.seekg(sizeof (float), ios::beg);
    double tmp_z = 0;
    int i;
    for (i = 0; i < num_particles_; i++)
    {
        float temp;
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].r = temp;
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].pos.x = temp;
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].pos.y = temp;
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].pos.z = temp;
        if (tmp_z < sand_array_[i].pos.z + sand_array_[i].r)
        {
            tmp_z = sand_array_[i].pos.z + sand_array_[i].r;
        }
        sand_array_[i].m = pow((sand_array_[i].r / radius_min_), 3) * mass_min_;
        /*
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].v.x = temp;
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].v.y = temp;
        sand.read(reinterpret_cast<char *> (&temp), sizeof (float));
        sand_array_[i].v.z = temp;
        */
    }
    cout << "max height of sand surface is " << tmp_z << endl;
    sand.close();
}

// generate sand initial position and velocity, wait to be deposited.

void
MDSand::GenerateSand()
{
    srand(time(NULL));
    int i, m;
    double inner_x_lb = box_length_x_ / 2.0 - 1.5;
    double inner_x_ub = box_length_x_ / 2.0 + 1.5;
    double inner_y_lb = box_length_y_ / 2.0 - 1.5;
    double inner_y_ub = box_length_y_ / 2.0 + 1.5;
    double inner_z_lb = 6.0;
    double lastx = inner_x_lb;
    double z_row = inner_z_lb, y_row = inner_y_lb;

    for (i = 0; i < num_particles_; i++)
    {
        m = rand();
        sand_array_[i].r =
            radius_min_ + (radius_max_ - radius_min_) * static_cast<double> (m) / RAND_MAX;
        sand_array_[i].m = pow((sand_array_[i].r / radius_min_), 3) * mass_min_;
        m = rand();
        sand_array_[i].pos.x =
            lastx + sand_array_[i].r * (1 + 0.1 * static_cast<double> (m) / RAND_MAX);
        if (sand_array_[i].pos.x + sand_array_[i].r > inner_x_ub)
        {
            sand_array_[i].pos.x =
                inner_x_lb + sand_array_[i].r * (1 + 0.1 * static_cast<double> (m) / RAND_MAX);
            lastx = inner_x_lb;
            y_row = y_row + 2 * radius_max_;
            if (y_row + 2 * radius_max_ > inner_y_ub)
            {
                m = rand();
                y_row =
                    inner_y_lb + sand_array_[i].r * (1 + 0.1 * static_cast<
                    double> (m) / RAND_MAX);
                z_row += 2 * radius_max_;
            }
        }
        lastx = sand_array_[i].pos.x + sand_array_[i].r;
        sand_array_[i].pos.y = y_row + sand_array_[i].r;
        sand_array_[i].pos.z = z_row + sand_array_[i].r;
        sand_array_[i].v = 0;
    }
}

// @and the interaction is the same as between small balls.

void
MDSand::InteractVirtualBall(int j,
    RealVec & shapePos,
    RealVec & shapeVel,
    RealVec & shapeForce, double shapeRadius, double shapeMass)
{
    RealVec dr = sand_array_[j].pos - shapePos;

    if (rsquare(dr) < (shapeRadius + sand_array_[j].r) * (shapeRadius + sand_array_[j].r))
    {
        RealVec dv, vn, vs, Fn, Fs, vs_direction, normal_direction;
        double scl_vn, abs_vs, abs_Fn, abs_Fs;
        double mr, abs_dr, overlap;

        mr = 1 / (1 / shapeMass + 1 / sand_array_[j].m);
        abs_dr = sqrt(rsquare(dr));
        overlap = shapeRadius + sand_array_[j].r - abs_dr;
        normal_direction = (sand_array_[j].pos - shapePos) / abs_dr;
        dv = sand_array_[j].v - shapeVel;
        scl_vn = dv * normal_direction;
        vn = scl_vn * normal_direction;
        vs = dv - vn;
        abs_vs = sqrt(vs * vs);

        if (abs_vs > 1e-8)
        {
            vs_direction = vs / abs_vs;
        }
        else
        {
            vs_direction = 0;
        }
        abs_Fn = sqrt(overlap) * (kn_ * overlap - gn_ * scl_vn);
        abs_Fs = mu_bp_ * abs_Fn;
        Fn = -abs_Fn * normal_direction;
        Fs = abs_Fs * vs_direction;
        shapeForce = Fn + Fs;
        sand_array_[j].F = sand_array_[j].F - Fn - Fs;
    }
    else
    {
        shapeForce = 0;
    }
}

// let a pair of balls interact, adding the force between them to each one's
// total force
// always use i<j in order to correctly calculate the slip term

void
MDSand::Interact(int i, int j)
{

    bool if_overlap = false, if_found = false, if_voc = false;
    int ref_voc;
    RealVec dr_ji, dv_ji, vn, vs, Fn, Fs, vs_direction, normal_direction;
    double abs_vn, abs_vs, abs_Fn, abs_Fs, abs_Fss, mr, abs_dr, overlap;

    dr_ji = sand_array_[j].pos - sand_array_[i].pos;
    mr = 1 / (1 / sand_array_[i].m + 1 / sand_array_[j].m);
    if_overlap =
        bool (rsquare(dr_ji) <
        (sand_array_[i].r + sand_array_[j].r) * (sand_array_[i].r + sand_array_[j].r));
    /*
    int k = 0;
    for (; k < max_contact_n_; k++)
    {
        if (if_voc == false)
            if (slip_array_[i][k].ref_num == 0)
            {
                if_voc = true;
                ref_voc = k;
                continue;
            }
        if (slip_array_[i][k].ref_num == j)
        {
            if_found = true;
            break;
        }
    }
    if (k == max_contact_n_ && if_voc == false && if_found == false)
    {
        cout << "code error" << endl;
        assert(0);
    }
    */
    // while (if_found == false && pslip
    // != NULL)
    // {
    // if (i == pslip->ref_num1 || i == pslip->ref_num2)
    // {
    // if (j == pslip->ref_num1 || j == pslip->ref_num2)
    // {
    // if_found = true;
    // }
    // else
    // {
    // previous = pslip;
    // pslip = pslip->next;
    // }
    // }
    // else if (pslip->ref_num1 == 0 && pslip->ref_num2 == 0)
    // {
    // if_voc = true;
    // pvoc = pslip;
    // previous = pslip;
    // pslip = pslip->next;
    // }
    // else
    // {
    // previous = pslip;
    // pslip = pslip->next;
    // }
    // }

    if (if_overlap)
    {
        abs_dr = sqrt(rsquare(dr_ji));
        overlap = sand_array_[i].r + sand_array_[j].r - abs_dr;
        normal_direction = (sand_array_[j].pos - sand_array_[i].pos) / abs_dr;
        dv_ji =
            sand_array_[j].v - sand_array_[i].v;// +
            // CrossProduct(sand_array_[j].omega,
            // -1.0 * dr_ji * sand_array_[j].r / abs_dr) -
            // CrossProduct(sand_array_[i].omega, 1.0 * dr_ji * sand_array_[i].r / abs_dr);
        abs_vn = dv_ji * normal_direction;
        vn = abs_vn * normal_direction;
        vs = dv_ji - vn;
        abs_vs = sqrt(vs * vs);

        if (abs_vs > 1e-8)
        {
            vs_direction = vs / abs_vs;
        }
        else
        {
            vs_direction = 0;
        }
        /*
        if (if_found == true)
        {
            RealVec disp_rot, disp_proj;

            disp_rot =
                slip_array_[i][k].disp -
                slip_array_[i][k].ndir * (slip_array_[i][k].disp * normal_direction);
            disp_proj = disp_rot - normal_direction * (normal_direction * disp_rot);
            slip_array_[i][k].disp = disp_proj + vs * step_tau_;
            slip_array_[i][k].ndir = normal_direction;
        }
        else
        {
            slip_array_[i][ref_voc].ref_num = j;
            slip_array_[i][ref_voc].disp = 0;
            slip_array_[i][ref_voc].ndir = normal_direction;
            k = ref_voc;
        }
        */
        // if (if_found == true)
        // {
        // pslip->ddisp = vs - normal_direction * (pslip->disp *
        // dv_ji) / abs_dr;
        // }
        // else
        // {
        // if (if_voc == false)
        // {
        // if (pslip != NULL)
        // {
        // cout << "error at searching for matching pair!" << endl;
        // assert(0);
        // }
        // pslip = (MDSlipDisp *) malloc(sizeof (struct MDSlipDisp));
        // if (slip_array_[cellnum] != NULL)
        // previous->next = pslip;
        // else
        // {
        // slip_array_[cellnum] = pslip;
        // }
        // pslip->next = NULL;
        // pslip->ref_num1 = i;
        // pslip->ref_num2 = j;
        // pslip->disp = 0;
        // pslip->ddisp =
        // vs - normal_direction * (pslip->disp * dv_ji) / abs_dr;
        // }
        // else
        // {
        // pslip = pvoc;
        // pslip->ref_num1 = i;
        // pslip->ref_num2 = j;
        // pslip->disp = 0;
        // pslip->ddisp = vs - normal_direction *
        // (pslip->disp * dv_ji) / abs_dr;
        // }
        // }
        // abs_Fss = kt_ * sqrt(slip_array_[i][k].disp *
        // slip_array_[i][k].disp);
        // cout << pslip->disp.x/abs_Fss << "\t" << pslip->disp.y/abs_Fss <<
        // "\t" << pslip->disp.z/abs_Fss << endl;
        // cout << vs_direction.x << "\t" << vs_direction.y << "\t" <<
        // vs_direction.z << endl;
        // cout << endl;
        abs_Fn = sqrt(overlap) * (kn_ * overlap - gn_ * abs_vn);
        // abs_Fs = min(abs_Fss, mu_pp_ * abs_Fn);
        abs_Fs = mu_pp_ * abs_Fn;
        Fn = -abs_Fn * normal_direction;
        Fs = abs_Fs * vs_direction;
        sand_array_[i].F = sand_array_[i].F + Fn + Fs;
        sand_array_[i].abF = sand_array_[i].abF + abs(abs_Fn);
        //sand_array_[i].moment =
        //    sand_array_[i].moment + CrossProduct(dr_ji * sand_array_[i].r / abs_dr, Fs);
        sand_array_[j].F = sand_array_[j].F - Fn - Fs;
        sand_array_[j].abF = sand_array_[j].abF + abs(abs_Fn);
        //sand_array_[j].moment =
        //    sand_array_[j].moment + CrossProduct(dr_ji * sand_array_[j].r / abs_dr, Fs);
        if (isnan(abs_Fn))
        {
            cout << sand_array_[i].pos.x << "\t" << sand_array_[i].pos.y << "\t" << sand_array_[i].pos.z << "\n";
            cout << sand_array_[j].pos.x << "\t" << sand_array_[j].pos.y << "\t" << sand_array_[j].pos.z << "\n";
            assert(!isnan(abs_Fn) && "Large Overlap");
        }
    }
    /*
    else
    {
        if (if_found == true)
        {
            slip_array_[i][k].ref_num = 0;
        }
    }*/
}

void
MDSand::SwitchOutput(bool output)
{
    if_output_ = output;
}

void
MDSand::OutputState()
{
    double temp = sys_time_;
    stringstream statename;

    statename << ".//timesnap=" << temp << "_statesave";
    para_file_.open(statename.str().c_str(), ios::out | ios::binary);
    para_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
    int count = 0;

    for (int i = 0; i < num_particles_; i++)
    {
        if (sand_array_[i].remove == 1)
        {
            count++;
            continue;
        }
        temp = (sand_array_[i].pos.x);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));

        temp = (sand_array_[i].pos.y);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));

        temp = (sand_array_[i].pos.z);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));

        temp = (sand_array_[i].r);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));

        temp = (sand_array_[i].v.x);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));

        temp = (sand_array_[i].v.y);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));

        temp = (sand_array_[i].v.z);
        para_file_.write(reinterpret_cast<char *> (&temp), sizeof (double));
    }
    para_file_.close();
}

void
MDSand::OutputFile()
{
    // Write to binary file with time. this record every balls position and
    // radius and abF. Big file, to make animation mainly.
    if (if_output_)
    {
        energy_file_ << sys_time_ << "\t";
        float temp = float (sys_time_);
        sand_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
        int count = 0;

        for (int i = 0; i < num_particles_; i++)
        {
            if (sand_array_[i].remove == 1)
            {
                count++;
                continue;
            }
            temp = float (sand_array_[i].r);
            sand_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
            temp = float (sand_array_[i].pos.x);
            sand_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
            temp = float (sand_array_[i].pos.y);
            sand_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
            temp = float (sand_array_[i].pos.z);
            sand_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
            //temp = float (sqrt(sand_array_[i].v * sand_array_[i].v));
            //temp = float (sand_array_[i].abF);
            //sand_file_.write(reinterpret_cast<char *> (&temp), sizeof (float));
        }
        sand_file_.flush();
        if (count > 0)
        {
            cout << num_particles_ - count << endl;
            sand_file_.close();
            exit(-1);
        }
        energy_file_ << sys_kinematic_energy_ << "\t" << sys_potential_Energy_
            << "\t" << sys_dissip_energy_ << "\t";
        energy_file_ << sys_energy_ << endl;
    }
}

void
MDSand::InteractExt(const MDShapeParser & shape_parser)
{
    bool collision;
    RealVec shape_force, shape_moment, virtual_pos, virtual_vel;
    RealVec tmp_force, tmp_pos, tmp_moment;
    double shape_mass, virtual_radius;

    for (int inode = 0; inode < shape_parser.n; inode++)
    {
        shape_parser.pshape[inode]->GetNodeMass(shape_mass);
        shape_parser.pshape[inode]->GetNodePosition(tmp_pos);
        shape_force = 0;
        shape_moment = 0;
        for (int i = 0; i < num_particles_; i++)
        {
            bool if_inside = 0;

            collision =
                shape_parser.pshape[inode]->GetSurfGeo(sand_array_[i].pos,
                sand_array_[i].r,
                virtual_pos,
                virtual_vel, virtual_radius, if_inside);
            if (if_inside == 1)
            {
                // some particles inside the shape
#if MD_GEN_SAND == 1
                sand_array_[i].remove = 1;
#else
                OutputState();
                assert(0 && (cout << "particles inside shape!" << endl));
#endif
            }
            if (if_inside == 0 && collision)
            {
                tmp_force = 0;
                InteractVirtualBall(i, virtual_pos, virtual_vel, tmp_force,
                    virtual_radius, shape_mass);
                shape_force = shape_force + tmp_force;
                tmp_moment = CrossProduct(virtual_pos - tmp_pos, tmp_force);
                shape_moment = shape_moment + tmp_moment;
                /*
                 * if (sqrt(tmp_force * tmp_force) > 1e5) { cout << sand_array_[i].r << endl; cout << sand_array_[i].pos.x
                 * <<"\t"<< sand_array_[i].pos.y <<"\t"<< sand_array_[i].pos.z << endl; cout << sand_array_[i].v.x <<"\t"<<
                 * sand_array_[i].v.y <<"\t"<< sand_array_[i].v.z << endl; cout << virtual_pos.x <<"\t"<< virtual_pos.y <<"\t"<<
                 * virtual_pos.z << endl; cout << virtual_vel.x <<"\t"<< virtual_vel.y <<"\t"<< virtual_vel.z << endl; cout <<
                 * tmp_pos.x <<"\t"<< tmp_pos.y <<"\t"<< tmp_pos.z << endl; exit(0); }
                 */
#if MD_GEN_SAND == 1
                if (sqrt(tmp_force * tmp_force) > 1e-8)
                {
                    sand_array_[i].remove = 1;
                }
#endif
            }
        }
        // if(shape_force.y!=0) {
        // cout << shape_force.y << endl;
        // cout << sys_step_ << endl;
        // cout << tmp_pos.y << endl;
        // assert(0);
        // }
        shape_parser.pshape[inode]->SetNodeForce(shape_force);
        shape_parser.pshape[inode]->SetNodeMoment(shape_moment);
    }
}

void
MDSand::ExForce()
{
    //    for (int j = 0; j < n_boulder_; j++)
    //    {
    //        int address = boulder_array_[j] - sand_array_;
    //
    //        // if (address != 72174 && address != 50542) {
    //        if (true)
    //        {
    //            for (int i = 0; i < num_particles_; i++)
    //            {
    //                if (i != address)
    //                    Interact(address, i);
    //            }
    //        }
    //    }
}

/*
 * run the system for one time step , update the force on each sand ball
 */
void
MDSand::RunOneStep()
{

    int i, k;
    IntVec ci, cu;
    MDCell *pair1, *pair2;

    // this part is to calculate the transition of particles between cells
    for (i = 0; i < cell_num_x_ * cell_num_y_ * cell_num_z_; i++)
    {
        pair1 = box_cell_[i];
        CellIndexToDim23(i, ci);
        while (pair1 != NULL)
        {
            int tempnum;

            tempnum = (*pair1).ref_num;
            PositionToCell(sand_array_[tempnum].pos, cu);
            if (!equal(ci, cu))
            {
                TransCell(tempnum, ci, cu);
                pair1 = box_cell_[i];
                continue;
            }
            pair1 = (*pair1).next;
        }
    }

    // now calculate the interaction pairs
    for (i = 0; i < cell_num_x_ * cell_num_y_ * cell_num_z_; i++)
        // for (ci = 0; !BoundaryEnd(ci); ++ci)
    {
        pair1 = box_cell_[i];
        CellIndexToDim23(i, ci);
        while (pair1 != NULL)
        {
            for (k = 0; k < MD_N_NBR; k++)
            {
                cu = ci + md_cell_nbr_[k];
                // jump when it is out of boundary
                if (RegionOut(cu))
                    continue;
                pair2 = box_cell_[CellIndexToDimOne(cu)];
                while (pair2 != NULL)
                {
                    // k==0 means the original cell(same cell) so we want
                    // pnumB > pnumA to avoid repeated pair
                    if ((k == 0 && (*pair1).ref_num < (*pair2).ref_num) || k != 0)
                    {
                        if ((*pair1).ref_num < (*pair2).ref_num)
                        {
                            Interact((*pair1).ref_num, (*pair2).ref_num);
                        }
                        else
                        {
                            Interact((*pair2).ref_num, (*pair1).ref_num);
                        }
                    }
                    pair2 = (*pair2).next;
                }
            }
            pair1 = (*pair1).next;
        }

    }

}

void
MDSand::BoundaryCondition(unsigned int controlnumber)
{

    int i;

    if (controlnumber == 0)
    {
        for (i = 0; i < num_particles_; i++)
        {
            // if (sandDisk[i].Position.x > Container_width) {
            // overlap = sandDisk[i].Position.x - Container_width;
            // Fn = (K * overlap +
            // mr * GN * fabs(sandDisk[i].Vx) *
            // Heaviside(sandDisk[i].Vx));
            // Fs = min(mr * GS * fabs(sandDisk[i].Vy), MU * fabs(Fn));
            // if (sandDisk[i].Vy > 0)
            // Fs = -Fs;
            // sandDisk[i].Fx += (-Fn);
            // sandDisk[i].Fy += Fs;
            // }
            // if (sandDisk[i].Position.x < 0) {
            // overlap = -sandDisk[i].Position.x;
            // Fn = (K * overlap +
            // mr * GN * fabs(sandDisk[i].Vx) *
            // Heaviside(-sandDisk[i].Vx));
            // Fs = min(mr * GS * fabs(sandDisk[i].Vy), MU * fabs(Fn));
            // if (sandDisk[i].Vy > 0)
            // Fs = -Fs;
            // sandDisk[i].Fx += Fn;
            // sandDisk[i].Fy += Fs;
            // }
            // if (sandDisk[i].Position.y < 0) {
            // overlap = -sandDisk[i].Position.y;
            // Fn = (K * overlap +
            // mr * GN * fabs(sandDisk[i].Vy) *
            // Heaviside(-sandDisk[i].Vy));
            // Fs = min(mr * GS * fabs(sandDisk[i].Vx), MU * fabs(Fn));
            // if (sandDisk[i].Vx > 0)
            // Fs = -Fs;
            // sandDisk[i].Fx += Fs;
            // sandDisk[i].Fy += Fn;
            // }
            // "sticky" wall condition:
            if (sand_array_[i].pos.x + sand_array_[i].r > box_length_x_)
            {
                sand_array_[i].pos.x = box_length_x_ - sand_array_[i].r;
                sys_dissip_energy_ =
                    sys_dissip_energy_ +
                    0.5 * sand_array_[i].m * (sand_array_[i].v.x * sand_array_[i].v.x);
                sand_array_[i].v.x = 0;
            }
            if (sand_array_[i].pos.x - sand_array_[i].r < 0)
            {
                sand_array_[i].pos.x = sand_array_[i].r;
                sys_dissip_energy_ =
                    sys_dissip_energy_ +
                    0.5 * sand_array_[i].m * (sand_array_[i].v.x * sand_array_[i].v.x);
                sand_array_[i].v.x = 0;
            }
            if (sand_array_[i].pos.y - sand_array_[i].r < 0)
            {
                sand_array_[i].pos.y = sand_array_[i].r;
                sys_dissip_energy_ =
                    sys_dissip_energy_ +
                    0.5 * sand_array_[i].m * (sand_array_[i].v.y * sand_array_[i].v.y);
                sand_array_[i].v.y = 0;
            }
            if (sand_array_[i].pos.y + sand_array_[i].r > box_length_y_)
            {
                sand_array_[i].pos.y = box_length_y_ - sand_array_[i].r;
                sys_dissip_energy_ =
                    sys_dissip_energy_ +
                    0.5 * sand_array_[i].m * (sand_array_[i].v.y * sand_array_[i].v.y);
                sand_array_[i].v.y = 0;
            }
            if (sand_array_[i].pos.z - sand_array_[i].r < 0)
            {
                sand_array_[i].pos.z = sand_array_[i].r;
                sys_dissip_energy_ =
                    sys_dissip_energy_ +
                    0.5 * sand_array_[i].m * (sand_array_[i].v.z * sand_array_[i].v.z);
                sand_array_[i].v.z = 0;
            }
            if (sand_array_[i].pos.z + sand_array_[i].r > box_length_z_)
            {
                sand_array_[i].pos.z = box_length_z_ - sand_array_[i].r;
                sys_dissip_energy_ =
                    sys_dissip_energy_ +
                    0.5 * sand_array_[i].m * (sand_array_[i].v.z * sand_array_[i].v.z);
                sand_array_[i].v.z = 0;
            }
            sys_kinematic_energy_ =
                sys_kinematic_energy_ +
                0.5 * sand_array_[i].m * (sand_array_[i].v * sand_array_[i].v);
            sys_potential_Energy_ =
                sys_potential_Energy_ + sand_array_[i].m * GRAVITY * sand_array_[i].pos.z;
        }
    }
    //    if (controlnumber == 1)
    //    {
    //        for (i = 0; i < num_particles_; i++)
    //        {
    //            if (sand_array_[i].pos.z - sand_array_[i].r < 0)
    //            {
    //                sand_array_[i].pos.z = sand_array_[i].r;
    //                sand_array_[i].v.z = 0;
    //            }
    //            if (sand_array_[i].pos.z + sand_array_[i].r > box_length_z_)
    //            {
    //                sand_array_[i].pos.z = box_length_z_ - sand_array_[i].r;
    //                sand_array_[i].v.z = 0;
    //            }
    //            double
    //                rho =
    //                sqrt ((sand_array_[i].pos.x -
    //                       3) * (sand_array_[i].pos.x - 3) +
    //                      (sand_array_[i].pos.y - 3) * (sand_array_[i].pos.y - 3));
    //            RealVec dr;
    //
    //            if (rho + sand_array_[i].r > 3)
    //            {
    //                sand_array_[i].pos.x =
    //                    3 + (sand_array_[i].pos.x - 3) * (3 - sand_array_[i].r) / rho;
    //                sand_array_[i].pos.y =
    //                    3 + (sand_array_[i].pos.y - 3) * (3 - sand_array_[i].r) / rho;
    //                dr.x = sand_array_[i].pos.x - 3;
    //                dr.z = 0;
    //                dr.y = sand_array_[i].pos.y - 3;
    //                sand_array_[i].v = sand_array_[i].v - dr * (sand_array_[i].v * dr) / (rho * rho);
    //            }
    //            if (sand_array_[i].pos.z > 3.5 && sand_array_[i].pos.z < 6 + 0.4142 * sand_array_[i].r)
    //            {
    //                double tmpr = sand_array_[i].pos.z - 3;
    //
    //                if (rho < tmpr)
    //                {
    //                    if (rho + sand_array_[i].r * 1.4142 > tmpr)
    //                    {
    //                        sand_array_[i].pos.x =
    //                            3 + (sand_array_[i].pos.x - 3) * (tmpr -
    //                                                              1.4142 * sand_array_[i].r) / rho;
    //                        sand_array_[i].pos.y =
    //                            3 + (sand_array_[i].pos.y - 3) * (tmpr -
    //                                                              1.4142 * sand_array_[i].r) / rho;
    //                        dr.x = sand_array_[i].pos.x - 3;
    //                        dr.z = -(tmpr - 1.4142 * sand_array_[i].r);
    //                        dr.y = sand_array_[i].pos.y - 3;
    //                        sand_array_[i].v =
    //                            sand_array_[i].v - dr * (sand_array_[i].v * dr) / (rho * rho * 2);
    //                    }
    //                }
    //                else
    //                {
    //                    if (rho - sand_array_[i].r * 1.4142 < tmpr)
    //                    {
    //                        sand_array_[i].pos.x =
    //                            3 + (sand_array_[i].pos.x - 3) * (tmpr +
    //                                                              1.4142 * sand_array_[i].r) / rho;
    //                        sand_array_[i].pos.y =
    //                            3 + (sand_array_[i].pos.y - 3) * (tmpr +
    //                                                              1.4142 * sand_array_[i].r) / rho;
    //                        dr.x = sand_array_[i].pos.x - 3;
    //                        dr.z = -(tmpr + 1.4142 * sand_array_[i].r);
    //                        dr.y = sand_array_[i].pos.y - 3;
    //                        sand_array_[i].v =
    //                            sand_array_[i].v - dr * (sand_array_[i].v * dr) / (rho * rho * 2);
    //                    }
    //                }
    //            }
    //            sys_kinematic_energy_ =
    //                sys_kinematic_energy_ +
    //                0.5 * sand_array_[i].m * (sand_array_[i].v * sand_array_[i].v);
    //            sys_potential_Energy_ =
    //                sys_potential_Energy_ + sand_array_[i].m * GRAVITY * sand_array_[i].pos.z;
    //        }
    //    }
}

void
MDSand::Integration()
{
    int i = 0;

    for (i = 0; i < num_particles_; i++)
    {
        //sand_array_[i].omega = sand_array_[i].omega + step_tau_
        //    * sand_array_[i].moment / (0.4 * sand_array_[i].m *
        //    sand_array_[i].r * sand_array_[i].r);
        sand_array_[i].F = sand_array_[i].F + sand_array_[i].m * gravity_vector;
        sand_array_[i].v = sand_array_[i].v + step_tau_ * sand_array_[i].F / sand_array_[i].m;
        sand_array_[i].pos = sand_array_[i].pos + sand_array_[i].v * step_tau_;
    }
    /*
     * MDSlipDisp* pslip, *previous; for (i = 0; i < MD_CELL_NUM_X * MD_CELL_NUM_Y * MD_CELL_NUM_Z; i++) { pslip =
     * slip_array_[i]; previous = pslip; while (pslip != NULL) { pslip->disp = pslip->disp + pslip->ddisp * step_tau_; previous =
     * pslip; pslip = pslip->next; } }
     */
}

void
MDSand::Run(int timeSteps, const MDShapeParser & shape_parser, bool flag)
{
    int i;
    int endTimeStep;

    endTimeStep = sys_step_ + timeSteps;
    for (; sys_step_ < endTimeStep;)
    {

        // clear the force on all the sand objects, to re-calculate
        for (i = 0; i < num_particles_; i++)
        {
            sand_array_[i].F = 0;
            sand_array_[i].moment = 0;
            sand_array_[i].abF = 0;
        }
        /*
         * run the system for one time step tau, update the force on each sand ball; If flag is set to true, also update the
         * force on the external nodes
         */
        RunOneStep();
        // ExForce();
        if (flag)
        {
            InteractExt(shape_parser);
#if MD_GEN_SAND == 1
            OutputFile();
            exit(0);
#endif
        }
        Integration();
        sys_step_++;
        sys_time_ = sys_step_ * step_tau_;
        sys_kinematic_energy_ = 0;
        sys_potential_Energy_ = 0;
        BoundaryCondition(0);
        sys_energy_ = sys_kinematic_energy_ + sys_potential_Energy_ + sys_dissip_energy_;
        if (sys_step_ % sav_intval_step_ == 0)
        {
            cout.flush();
            cout << "time=" << fixed << setprecision(5) << sys_time_ << "\r";
            OutputFile();
        }
    }
}

MDSand::~MDSand()
{
    MDCell *pointa, *pointb;

    for (int i = 0; i < cell_num_x_ * cell_num_y_ * cell_num_z_; i++)
    {
        pointa = box_cell_[i];
        while (pointa != NULL)
        {
            pointb = (*pointa).next;
            free(pointa);
            pointa = pointb;
        }
    }
    delete[]sand_array_;
    delete[]box_cell_;
    delete[]md_cell_nbr_;
    for (int i = 0; i < num_particles_; i++) {
	    delete slip_array_[i];
    }
    delete[]slip_array_;
}
