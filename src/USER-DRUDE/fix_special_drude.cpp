//#include "string.h"
//#include "stdlib.h"
#include "math.h"
#include "fix_special_drude.h"
#include "atom.h"
//#include "force.h"
#include "comm.h"
//#include "input.h"
//#include "variable.h"
//#include "random_mars.h"
//#include "group.h"
#include "update.h"
//#include "modify.h"
#include "memory.h"
//#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixSpecialDrude::FixSpecialDrude(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), done(false)
{}

/* ---------------------------------------------------------------------- */

void FixSpecialDrude::init()
{
  int dummy;
  index_drude = atom->find_custom("drude", dummy);
  if (index_drude == -1) 
    error->all(FLERR,"Unable to get DRUDE atom property");
  index_ref = atom->find_custom("ref", dummy);
  if (index_ref == -1) 
    error->all(FLERR,"Unable to get REF atom property");
}

/* ---------------------------------------------------------------------- */

int FixSpecialDrude::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpecialDrude::pre_force(int)
{
    
    //if (done) return;
    // Redefine lists of special neighbors so that Drudes are equivalent to their core
    // Need to call this after the standard list of special neighbors is built
    int nlocal = atom->nlocal, nall = atom->nlocal + atom->nghost;
    int *drude = atom->ivector[index_drude];
    int *ref = atom->ivector[index_ref];
    int **nspecial = atom->nspecial;
    tagint **special = atom->special;
    int j;
    int maxspecial = atom->maxspecial; // Normally there is no need to grow special if the special_bond extra keyword is used correctly
    /*std::cout << comm->me << " : " << maxspecial << std::endl;
    for (int i = 0; i < nall; i++) {
        std::cout << comm->me << " (avant) : " << atom->tag[i] << " (";
        if (i < nlocal) std::cout << "owned ";
        else std::cout << "ghost ";
        std::cout << i << ") ";
        for (int k = 0; k<atom->nspecial[i][2]; k++){
            if (k == atom->nspecial[i][0]) std::cout << "| "; 
            if (k == atom->nspecial[i][1]) std::cout << "| "; 
            std::cout << atom->special[i][k] << " ";
        } std::cout << std::endl;
    }*/
    
    
    /*for (int i = 0; i < nall; i++) {
        if (!ref[i]) continue;
        if (i < nlocal && atom->map(ref[i]) < nlocal) continue;
        if (i >= nlocal && atom->map(ref[i]) >= nlocal) continue;
        std::cout << comm->me << " : " << atom->tag[i] << " (";
        if (i < nlocal) std::cout << "owned ";
        else std::cout << "ghost ";
        std::cout << i << ") ";
        for (int k = 0; k<atom->nspecial[i][2]; k++){
            if (k == atom->nspecial[i][0]) std::cout << "| "; 
            if (k == atom->nspecial[i][1]) std::cout << "| "; 
            std::cout << atom->special[i][k] << " ";
        } std::cout << std::endl;
    }*/
    if (done) return;
    // All atoms need to know their special info.    
    // All atoms need to know their drudeness and the ref of their partner.
    // Loop on the core atoms
    for (int i=0; i < nall; i++){
        if (!drude[i]) { // Core atom (polarizable or not)
            // Remove Drude particles from the list because some are missing
            // and all of them are 1 neighbor too far away
            if (i >= nlocal) { // This is a ghost
                if (ref[i] == 0 
                  || atom->map(ref[i]) >= nlocal
                  || atom->map(ref[i]) < 0) continue;
                // This is either a non-polarizable atom
                // or a core whose Drude is also not owned
                // (ghost or unknown).
            }
            for (int k=0; k<nspecial[i][2]; k++) {
                j = atom->map(special[i][k]);
                if (drude[j]) { // if not j > 0, error "missing atom"
                    for (int l=k; l<nspecial[i][2]-1; l++) { // left shift
                        special[i][l] = special[i][l+1];
                    }
                    nspecial[i][2]--;
                    if (k < nspecial[i][1]){
                        nspecial[i][1]--;
                        if (k < nspecial[i][0]) {
                            nspecial[i][0]--;
                        }
                    }
                    k--;
                }
            }
        /*std::cout << comm->me << " (avant ajout) : " << atom->tag[i] << " (";
        if (i < nlocal) std::cout << "owned ";
        else std::cout << "ghost ";
        std::cout << i << ") ";
        for (int k = 0; k<atom->nspecial[i][2]; k++){
            if (k == atom->nspecial[i][0]) std::cout << "| "; 
            if (k == atom->nspecial[i][1]) std::cout << "| "; 
            std::cout << atom->special[i][k] << " ";
        } std::cout << std::endl;*/
            // Add back all Drude particles at the correct positions
            for (int k=0; k<nspecial[i][2]; k++) {
                j = atom->map(special[i][k]);
                if (ref[j] > 0) { // if not j > 0, error "missing atom"
                    // This neighbor is the core of a polarizable atom.
                    // Check size: need to grow special?
                    if (nspecial[i][2] == maxspecial) {
                        maxspecial++;
                        memory->grow(atom->special, atom->nmax, maxspecial, "atom:special");
                    }
                    for (int l=nspecial[i][2]; l>k+1; l--) { // right shift
                        special[i][l] = special[i][l-1];
                    }
                    special[i][k+1] = ref[j];
                    nspecial[i][2]++;
                    if (k < nspecial[i][1]) {
                        nspecial[i][1]++;
                        if (k < nspecial[i][0]) {
                            nspecial[i][0]++;
                        }
                    }
                    k++; // Must jump the new Drude particle
                }
            }
        /*std::cout << comm->me << " (avant self) : " << atom->tag[i] << " (";
        if (i < nlocal) std::cout << "owned ";
        else std::cout << "ghost ";
        std::cout << i << ") ";
        for (int k = 0; k<atom->nspecial[i][2]; k++){
            if (k == atom->nspecial[i][0]) std::cout << "| "; 
            if (k == atom->nspecial[i][1]) std::cout << "| "; 
            std::cout << atom->special[i][k] << " ";
        } std::cout << std::endl;*/
            // Also adds the own Drude partner at the beginning (if any)
            if (ref[i] > 0){
                // Check size: need to grow special?
    //std::cout << comm->me << " " << __LINE__ << "  " << nspecial[i][2] << std::endl;
                if (nspecial[i][2] == maxspecial) {
                    maxspecial++;
                    memory->grow(atom->special, atom->nmax, maxspecial, "atom:special");
                }
                for (int l=nspecial[i][2]; l>0; l--) { // right shift
    //std::cout << comm->me << " " << __LINE__ << "  " << nspecial[i][2] << std::endl;
                    special[i][l] = special[i][l-1];
                }
                special[i][0] = ref[i];
                nspecial[i][0]++;
                nspecial[i][1]++;
                nspecial[i][2]++;
            }
        }
    }
    
    /*std::cout << comm->me << " : " << maxspecial << std::endl;
    for (int i = 0; i < nall; i++) {
        std::cout << comm->me << " (avant comm) : " << atom->tag[i] << " (";
        if (i < nlocal) std::cout << "owned ";
        else std::cout << "ghost ";
        std::cout << i << ") ";
        for (int k = 0; k<atom->nspecial[i][2]; k++){
            if (k == atom->nspecial[i][0]) std::cout << "| "; 
            if (k == atom->nspecial[i][1]) std::cout << "| "; 
            std::cout << atom->special[i][k] << " ";
        } std::cout << std::endl;
    }*/
    // Communicate the new size of special.
    MPI_Allreduce(&maxspecial, &atom->maxspecial, 1, MPI_INT, MPI_MAX, world);
    if (maxspecial != atom->maxspecial) {
      maxspecial = atom->maxspecial;
      memory->grow(atom->special, atom->nmax, maxspecial, "atom:special");
    }
    // Let the ghost cores whose Drude is owned know their special lists. 
    comm_forward = maxspecial + 3;
    comm->forward_comm_fix(this);
    
    for (int i=0; i<nlocal; i++) {
        // Copy the Drude particle list of special neighbors from its core.
        if (drude[i]) {
            j = atom->map(ref[i]);
            nspecial[i][0] = nspecial[j][0];
            nspecial[i][1] = nspecial[j][1];
            nspecial[i][2] = nspecial[j][2];
            special[i][0] = ref[i]; // First is the core, instead of itself
            for (int k=1; k<nspecial[i][2]; k++){
                special[i][k] = special[j][k];
            }
        }
    } 
    // reset ghosts
    for (int i=nlocal; i<nall; i++) nspecial[i][2] = 0; 
    
    /*std::cout << comm->me << " : " << maxspecial << std::endl;
    for (int i = 0; i < nall; i++) {
        std::cout << comm->me << " (apres) : " << atom->tag[i] << " (";
        if (i < nlocal) std::cout << "owned ";
        else std::cout << "ghost ";
        std::cout << i << ") ";
        for (int k = 0; k<atom->nspecial[i][2]; k++){
            if (k == atom->nspecial[i][0]) std::cout << "| "; 
            if (k == atom->nspecial[i][1]) std::cout << "| "; 
            std::cout << atom->special[i][k] << " ";
        } std::cout << std::endl;
    }*/
    done = true;
}

/* ---------------------------------------------------------------------- */

int FixSpecialDrude::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int m = 0;
  for (int i=0; i<n; i++) {
    int j = list[i];
    buf[m++] = ubuf(atom->nspecial[j][0]).d;
    buf[m++] = ubuf(atom->nspecial[j][1]).d;
    buf[m++] = ubuf(atom->nspecial[j][2]).d;
    for (int k=0; k<atom->nspecial[j][2]; k++) {
      buf[m++] = ubuf(atom->special[j][k]).d;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixSpecialDrude::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i=first; i<last; i++) {
    atom->nspecial[i][0] = (int) ubuf(buf[m++]).i;
    atom->nspecial[i][1] = (int) ubuf(buf[m++]).i;
    atom->nspecial[i][2] = (int) ubuf(buf[m++]).i;
    for (int k=0; k<atom->nspecial[i][2]; k++) {
      atom->special[i][k] = (int) ubuf(buf[m++]).i;
    }
  }
}




