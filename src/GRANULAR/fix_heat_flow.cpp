/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_heat_flow.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "input.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, CONSTANT, TYPE };
enum { NO_SOURCE, CONSTANT_SOURCE, EQUAL_SOURCE, ATOM_SOURCE };
enum { NO_RXN, CONSTANT_RXN, EQUAL_RXN, ATOM_RXN };

/* ---------------------------------------------------------------------- */

FixHeatFlow::FixHeatFlow(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), sourceval(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix heat/flow", error);

  cp_style = NONE;
  reaction = 0;
  nreaction = 0;
  comm_forward = 1;
  comm_reverse = 1;

  int ntypes = atom->ntypes;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "constant") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix command");
      cp_style = CONSTANT;
      cp = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      if (cp < 0.0) error->all(FLERR, "Illegal fix command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "type") == 0) {
      if (iarg + ntypes + 1 > narg) error->all(FLERR, "Illegal fix command");
      cp_style = TYPE;
      memory->create(cp_type, ntypes + 1, "fix/temp/integrate:cp_type");
      for (int i = 1; i <= ntypes; i++) {
        cp_type[i] = utils::numeric(FLERR, arg[iarg + i], false, lmp);
        if (cp_type[i] < 0.0) error->all(FLERR, "Illegal fix command");
      }
      iarg += ntypes + 1;
    } else if (strcmp(arg[iarg], "source") == 0) {
      source = 1;
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix heat/flow source", error);
      if (utils::strmatch(arg[iarg + 1], "^v_")) {
        source_str = utils::strdup(arg[iarg + 1] + 2);
      } else {
        source_value = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
        source_style = CONSTANT_SOURCE;
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "reaction") == 0) {
      //Args are :
      // - number of reactions (nreact)
      // - name of custom property/atom vector (or array if nreact > 1)
      // - per-atom or equal-style variable names or constants that define reaction rates (nreact repeats)
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix heat/flow reaction", error);
      reaction = 1;
      nreaction = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (iarg + nreaction + 3 > narg)
        utils::missing_cmd_args(FLERR, "fix heat/flow reaction", error);     
      
      int is_double, cols;
      reaction_index = atom->find_custom(arg[iarg + 2], is_double, cols);
      if (reaction_index == -1) {
        error->all(FLERR, "Fix heat/flow reaction requires previously defined property/atom");
      }
      if (cols != nreaction) {
        error->all(FLERR,
                   "Fix heat/flow reaction requires property/atom array with number of columns "
                   "that match number of reactions");
      }     

      reaction_str = new char *[nreaction];
      reaction_style = new int[nreaction];
      reaction_var = new int[nreaction];
      reaction_value = new double[nreaction];
      for (int i = 0; i < nreaction; i++) {
        if (utils::strmatch(arg[iarg + 3 + i], "^v_")) {
          reaction_str[i] = utils::strdup(arg[iarg + 3 + i] + 2);
        } else {
          reaction_value[i] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
          reaction_style[i] = CONSTANT_RXN;
        }
      }
      iarg += 3+nreaction;
    } else {
      error->all(FLERR, "Unknown fix heat/flow keyword {}", arg[3]);
    }
  }

  if (cp_style == NONE) error->all(FLERR, "Must specify specific heat in fix heat/flow");
  dynamic_group_allow = 1;

  maxatom = 1;
  memory->create(sourceval, maxatom, "heatflow:sourceval");
  memory->create(rxn_array, maxatom, nreaction, "heatflow:rxn_array");
}

/* ---------------------------------------------------------------------- */

int FixHeatFlow::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::init()
{
  dt = update->dt;

  if (!atom->temperature_flag)
    error->all(FLERR, "Fix heat/flow requires atom style with temperature property");
  if (!atom->heatflow_flag)
    error->all(FLERR, "Fix temp/integrate requires atom style with heatflow property");

  if (source_str) {
    source_var = input->variable->find(source_str);
    if (source_var < 0)
      error->all(FLERR, "Variable {} for fix heat/flow source does not exist", source_str);
    if (input->variable->equalstyle(source_var))
      source_style = EQUAL_SOURCE;
    else if (input->variable->atomstyle(source_var)) {
      source_style = ATOM_SOURCE;
    } else
      error->all(FLERR, "Variable {} for fix heat/flow source is invalid style", source_str);
  }
  if (reaction) {
    reaction_atom = 0;
    for (int i = 0; i < nreaction; i++) {
      if (reaction_str[i]) {
        reaction_var[i] = input->variable->find(reaction_str[i]);
        if (reaction_var[i] < 0)
          error->all(FLERR, "Variable {} for fix heat/flow reaction does not exist",
                     reaction_str[i]);
        if (input->variable->equalstyle(reaction_var[i]))
          reaction_style[i] = EQUAL_RXN;
        else if (input->variable->atomstyle(reaction_var[i])) {
          reaction_style[i] = ATOM_RXN;
          reaction_atom = 1;
        } else
          error->all(FLERR, "Variable {} for fix heat/flow reaction is invalid style",
                     reaction_str[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::setup_pre_force(int /*vflag*/)
{
  // Identify whether this is the first instance of fix heat/flow
  first_flag = 0;
  int i = 0;
  auto fixlist = modify->get_fix_by_style("heat/flow");
  for (const auto &ifix : fixlist) {
    if (strcmp(ifix->id, id) == 0) break;
    i++;
  }
  if (i == 0) first_flag = 1;
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::pre_force(int /*vflag*/)
{
  // send updated temperatures to ghosts if first instance of fix
  // then clear heatflow for next force calculation
  double *heatflow = atom->heatflow;
  if (first_flag) {
    comm->forward_comm(this);
    for (int i = 0; i < atom->nmax; i++) heatflow[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::final_integrate()
{
  // update temperature of atoms in group

  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // add ghost contributions to heatflow if first instance of fix
  if (first_flag) comm->reverse_comm(this);

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) temperature[i] += dt * heatflow[i] / (calc_cp(i) * rmass[i]);
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) temperature[i] += dt * heatflow[i] / (calc_cp(i) * mass[type[i]]);
  }

  if (source) {
    // Reallocate sourceval array if necessary
    if ((source_style == ATOM_SOURCE) && (atom->nmax > maxatom)) {
      maxatom = atom->nmax;
      memory->destroy(sourceval);
      memory->create(sourceval, maxatom, "heatflow:sourceval");
    }

    if (source_style == EQUAL_SOURCE) {
      source_value = input->variable->compute_equal(source_var);
    } else if (source_style == ATOM_SOURCE) {
      input->variable->compute_atom(source_var, igroup, sourceval, 1, 0);
    }
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (source_style == ATOM_SOURCE) source_value = sourceval[i];
        temperature[i] += dt * source_value;
      }
    }
  }

  if (reaction) {
    // Rellocate reaction rate arrays if necessary
    if ((reaction_atom) && (atom->nmax > maxatom)) {
      maxatom = atom->nmax;
      memory->destroy(rxn_array);
      memory->create(rxn_array, maxatom, nreaction, "heatflow:reaction");
    }

    for (int j = 0; j < nreaction; j++) {
      if (reaction_style[j] == ATOM_RXN) {
        input->variable->compute_atom(reaction_var[j], igroup, &rxn_array[0][j], nreaction, 0);
      }
      else if (reaction_style[j] == EQUAL_RXN){
        reaction_value[j] = input->variable->compute_equal(reaction_var[j]);
      }
    }

    // Update extents of reaction
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        for (int j = 0; j < nreaction; j++) {
          if (reaction_style[j] == ATOM_RXN) reaction_value[j] = rxn_array[i][j];        
<<<<<<< HEAD
          atom->darray[reaction_index][i][j] += dt * reaction_value[j];
=======
          reaction_extents[i][j] += dt * reaction_value[j];
>>>>>>> d21518e47d727d671f29ea5c45f892220c2a37e6
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::final_integrate_respa(int /*ilevel*/, int /*iloop*/)
{
  dt = update->dt;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::reset_dt()
{
  dt = update->dt;
}

/* ---------------------------------------------------------------------- */

double FixHeatFlow::calc_cp(int i)
{
  if (cp_style == TYPE) {
    return cp_type[atom->type[i]];
  } else {
    return cp;
  }
}

/* ---------------------------------------------------------------------- */

int FixHeatFlow::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  double *temperature = atom->temperature;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = temperature[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;

  double *temperature = atom->temperature;

  for (i = first; i < last; i++) temperature[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int FixHeatFlow::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  double *heatflow = atom->heatflow;

  for (int i = first; i < last; i++) { buf[m++] = heatflow[i]; }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  double *heatflow = atom->heatflow;

  for (int i = 0; i < n; i++) heatflow[list[i]] += buf[m++];
}
