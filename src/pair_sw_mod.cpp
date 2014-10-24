/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Modified for MoS2 according to JW Jiang, J Appl Phys 114, 064307 (2013)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_sw_mod.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWMod::PairSWMod(LAMMPS *lmp) : PairSW(lmp)
{};

PairSWMod::~PairSWMod()
{};

/* ---------------------------------------------------------------------- */

void PairSWMod::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2,rsq3;
  double delr1[3],delr2[3],delr3[3],fj[3],fk[3];
  double bondss = 3.78;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtag = tag[j];

      if (itag > jtag) {
        if ((itag+jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag+jtag) % 2 == 1) continue;
      } else {
        if (x[j][2] < ztmp) continue;
        if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
        if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = map[type[j]];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      ijparam = elem2param[itype][jtype][jtype];
      if (rsq > params[ijparam].cutsq) continue;

      twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    jnumm1 = jnum - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if (rsq1 > params[ijparam].cutsq) continue;

      for (kk = jj+1; kk < jnum; kk++) {
        k = jlist[kk];
        k &= NEIGHMASK;
        ktype = map[type[k]];
        ikparam = elem2param[itype][ktype][ktype];
        ijkparam = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        if (rsq2 > params[ikparam].cutsq) continue;

        /* For angle 123, r12 and r13 are two edges. I add an additional 
           requirement for r23<=bondss. This is necessary when there are
           more than 4 first-nearest-neighboring atoms for atom 1. For
           example, in MoS2, there are six S atoms around the Mo atom.
           Two different angles are formed by the Mo atom and its 6 
           neighboring S atoms. I only want to add angle interaction
           for one of these two angle. The other angle will be excluded
           by the additional requirement.(29/05/13/Wed/jjw)*/
        delr3[0] = x[k][0] - x[j][0];
        delr3[1] = x[k][1] - x[j][1];
        delr3[2] = x[k][2] - x[j][2];
        rsq3 = delr3[0]*delr3[0] + delr3[1]*delr3[1] + delr3[2]*delr3[2];
        if (rsq3 > bondss) continue;

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        f[i][0] -= fj[0] + fk[0];
        f[i][1] -= fj[1] + fk[1];
        f[i][2] -= fj[2] + fk[2];
        f[j][0] += fj[0];
        f[j][1] += fj[1];
        f[j][2] += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

