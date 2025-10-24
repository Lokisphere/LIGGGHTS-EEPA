/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

   Contributing authors:
   Rahul Mohanty (University of Edinburgh, P&G)
   Tomaz M. Zorec (University of Ljubljana)

-------------------------------------------------------------------------- 
   Modified and validated by Lokeshwar Mahto (Lokeshwar Mahto)
------------------------------------------------------------------------- */

#ifdef NORMAL_MODEL
NORMAL_MODEL(EEPA,eepa,15)
#else
#ifndef NORMAL_MODEL_EEPA_H_
#define NORMAL_MODEL_EEPA_H_
#include "contact_models.h"
#include "normal_model_base.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<EEPA> : public NormalModelBase
  {
  public:
    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase *c) :
      NormalModelBase(lmp, hsetup, c),
      Yeff(NULL),
      Geff(NULL),
      CoeffRestLog(NULL),
      betaeff(NULL),
      kn2kc(NULL),
      kn2k1(NULL),
      cex(0.0),
      dex(0.0),
      f_adh(NULL),
      gamma_surf(NULL),
      history_offset(0),
      tangential_damping(false),
      limitForce(false),
      cmb(c),
      ktm(NULL)
 	{
 	}
    inline void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce,true);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0");
      hsetup->add_history_value("old_delta", "0");
    }

    void connectToProperties(PropertyRegistry & registry) {

      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model eepa");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model eepa"); 
      registry.registerProperty("CoeffRestLog", &MODEL_PARAMS::createCoeffRestLog, "model eepa");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model eepa");
      registry.registerProperty("kn2k1", &MODEL_PARAMS::createUnloadingStiffness, "model eepa");
      registry.registerProperty("cex", &MODEL_PARAMS::createAdhesionExponent, "model eepa");
      registry.registerProperty("dex", &MODEL_PARAMS::createOverlapExponent, "model eepa");
      registry.registerProperty("f_adh", &MODEL_PARAMS::createPullOffForce, "model eepa");
      registry.registerProperty("gamma_surf", &MODEL_PARAMS::createSurfaceEnergy, "model eepa");
      registry.registerProperty("ktm", &MODEL_PARAMS::createTangentialMultiplier, "model eepa");   
      
      registry.connect("Yeff", Yeff,"model eepa");
      registry.connect("Geff", Geff,"model eepa");
      registry.connect("CoeffRestLog", CoeffRestLog, "model eepa");
      registry.connect("betaeff", betaeff,"model eepa");
      registry.connect("kn2k1", kn2k1,"model eepa");
      registry.connect("cex", cex,"model eepa");
      registry.connect("dex", dex,"model eepa");
      registry.connect("f_adh", f_adh,"model eepa");
      registry.connect("gamma_surf", gamma_surf,"model eepa");
      registry.connect("ktm", ktm,"model eepa");  
    }

    // effective exponent for stress-strain relationship

    inline double stressStrainExponent()
    {
      return dex;
    }

    inline double calculate_deltan_p_max (double deltan_p, double * const history, int count_flag, const double k2, double dex, double dex_i, double k_adh)
    {
      //calculating the maximum overlap
      double deltan_p_max;

      if (count_flag == 0 ) {
           if (deltan_p > history[0]) {
            deltan_p_max = deltan_p;
            history[0] = deltan_p_max;
           } else {
            deltan_p_max = history[0];
          }
      } else {
        deltan_p_max = pow((pow(history[1], dex) + ((k_adh/k2)*pow(history[1], cex))), dex_i);
        history[0] = deltan_p_max;
      }
      return history[0];
    }

    inline double calc_f_min(double k2, double f_min_lim, const double g_surf, const double f_0, double deltan_p_max, double deltan_pe_max, double reff)
    {
    	double a = sqrt(2 * deltan_p_max * reff);
    	double f_min = (1.5 * M_PI * g_surf * a) - f_0;
    	if (f_min > f_min_lim)
    	{
    		f_min = (f_min_lim - f_0) * 0.5;
    	}
    	return f_min;
    }

    inline double calculate_k_adh(double f_min, double deltan_pe_max, const double k2, const double dex_i, const double cex, const double f_0)
    {
	double new_k_c, delta_min;
	
       delta_min = pow((-f_min - f_0 + k2 * deltan_pe_max)/k2, dex_i);
       if (delta_min != 0) 
       {
       	new_k_c = (f_min + f_0) / pow(delta_min, cex);
       } else {
       	new_k_c = 0;
       	}
      return new_k_c;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = sidata.itype;
      const int jtype = sidata.jtype;
      const double deltan = sidata.deltan;
      double ri = sidata.radi;
      double rj = sidata.radj;
      double reff=sidata.is_wall ? sidata.radi : ((ri*rj)/(ri+rj));
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(sidata.is_non_spherical && atom->superquadric_flag)
            reff = sidata.reff;
#endif
      double meff=sidata.meff;
      const double f_0 = f_adh[itype][jtype];
      double sqrtval = sqrt(reff*sidata.deltan);
      double kn, kt;

      if(dex==1){
        kn = 2.0*Yeff[itype][jtype]*reff;
        kt = kn * ktm[itype][jtype];
      } else {
        kn = (4./3.)*Yeff[itype][jtype]*sqrt(reff);
        kt = 8 * Geff[itype][jtype] * sqrtval * ktm[itype][jtype];
        }

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double dex_i = 1./dex;
      const double k1 = kn;
      const double k2 = kn*kn2k1[itype][jtype];


      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_NORMAL_MODEL;
      double * const history = &sidata.contact_history[history_offset];

      double Fn_contact, deltan_e, deltan_ce, f_min_lim, fTmp, f_min;
      double k_adh = 0.0;

        const double lambda = pow((1. - k1/k2), dex_i);
        int count_flag = 0;

        double deltan_p = lambda*deltan;
        double deltan_p_max, deltan_pe_max;

        deltan_e = pow(deltan, dex);
        deltan_ce = pow(deltan, cex);


        deltan_p_max = calculate_deltan_p_max(deltan_p, &sidata.contact_history[history_offset], count_flag, k2, dex, dex_i, k_adh);

        // Normal force calculation for Edinburgh model

        const double g_surf = gamma_surf[itype][jtype];

        temp_calc:

        deltan_pe_max = pow(deltan_p_max, dex);
        f_min_lim = (k2 * deltan_pe_max) - f_0;
        f_min = calc_f_min(k2, f_min_lim, g_surf, f_0, deltan_p_max, deltan_pe_max, reff);
        k_adh = calculate_k_adh(f_min, deltan_pe_max, k2, dex_i, cex, f_0);
        
        fTmp=k2*(deltan_e - deltan_pe_max);

        if (fTmp >= k1 * deltan_e){ // loading
          Fn_contact = f_0 + k1 * deltan_e;
        }else{
         if (fTmp > (-k_adh * deltan_ce)){
            Fn_contact = f_0 + fTmp;
          }else{  // cohesion part

            if (deltan > history[1]){
              count_flag = 1;
              deltan_p_max = calculate_deltan_p_max(deltan_p, &sidata.contact_history[history_offset], count_flag, k2, dex, dex_i, k_adh);
              goto temp_calc;
            }
            Fn_contact = f_0 - k_adh * deltan_ce;
          }
        }
        history[1] = deltan;
      
      
      double gamman, gammat;
      gamman = sqrt(4.*meff*kn/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));
      gammat = sqrt(4.*meff*kt/(1.+(M_PI/CoeffRestLog[itype][jtype])*(M_PI/CoeffRestLog[itype][jtype])));

	if (!tangential_damping) gammat=0;
	
      const double Fn_damping = -gamman*sidata.vn;
      double Fn = Fn_damping + Fn_contact;

      sidata.Fn = Fn;
      sidata.Fn_contact = Fn_contact;
      sidata.kn = kn;
      sidata.kt = kt;
      sidata.gamman = gamman;
      sidata.gammat = gammat;
      sidata.f_min = f_min;

      // apply normal force
      if(sidata.is_wall) {
        const double Fn_ = Fn * sidata.area_ratio;
        i_forces.delta_F[0] = Fn_ * sidata.en[0];
        i_forces.delta_F[1] = Fn_ * sidata.en[1];
        i_forces.delta_F[2] = Fn_ * sidata.en[2];
      } else {
        i_forces.delta_F[0] = sidata.Fn * sidata.en[0];
        i_forces.delta_F[1] = sidata.Fn * sidata.en[1];
        i_forces.delta_F[2] = sidata.Fn * sidata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    inline void surfacesClose(SurfacesCloseData & scdata, ForceData&, ForceData&)
    {
      if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_NORMAL_MODEL;
      double * const history = &scdata.contact_history[history_offset];
      history[0] = 0.0;
      history[1] = 0.0;
    }

    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

  protected:
    double **Yeff;
    double **Geff;
    double **CoeffRestLog;
    double **betaeff;
    double **kn2kc;
    double **kn2k1;
    double cex;
    double dex;
    
    double **f_adh;
    double **gamma_surf;
    double **ktm;
    int history_offset;
    bool tangential_damping;
    bool limitForce;
    class ContactModelBase *cmb;
  };
}
}
#endif
#endif
