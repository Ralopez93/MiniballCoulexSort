#ifndef __ADDBACK_CXX
#define __ADDBACK_CXX

#define PBINS 800
#define PRANGE 800
#define PMIN -1.0 * PRANGE / PBINS
#define PMAX PRANGE + PMIN

#define GBINS 4000
#define GRANGE 4000
#define GMIN -1.0 * GRANGE / GBINS
#define GMAX GRANGE + GMIN

#define ELBINS 2000
#define ELRANGE 2000
#define ELMIN -1.0 * ELRANGE / ELBINS
#define ELMAX ELRANGE + ELMIN

#ifndef AddBack_hh
#include "AddBack.hh"
#endif

void AddBack::ClearEvt() {
  gen_array.clear();
  gtd_array.clear();
  clu_array.clear();
  cid_array.clear();
  sid_array.clear();
  sen_array.clear();
  ab_array.clear();

  vector<float>().swap(gen_array);
  vector<long long>().swap(gtd_array);
  vector<unsigned short>().swap(clu_array);
  vector<unsigned short>().swap(cid_array);
  vector<unsigned short>().swap(sid_array);
  vector<float>().swap(sen_array);
  vector<bool>().swap(ab_array);

  return;
}

void AddBack::MakeGammaRays(bool addback, bool reject, bool segsum) {
  // Reset some variables
  ab_mul = 0;

  // ------------------------------------------------------------------------ //
  // Only DGF events
  // ------------------------------------------------------------------------ //

  for (unsigned int j = 0; j < event->NumberOfDgfs(); j++) {
    reject_evt = false;
    ab_evt = false;
    veto_gamma = false;
    dgf_num = event->Dgf(j)->ModuleNumber();
    dgf_ch = event->Dgf(j)->Channel();
    dgf_en = event->Dgf(j)->Energy();
    dgf_t = event->Dgf(j)->Time() + Cal->CryTime(dgf_num / 2);

    // Miniball, this is where we do real gammas and addback
    if (0 <= dgf_num && dgf_num < 48 && 0 <= dgf_ch && dgf_ch < 4) {
      GammaEnergy = Cal->DgfEnergy(dgf_num, dgf_ch, dgf_en);

      if (dgf_num % 2 == 0 && dgf_ch < 3) {  // cores plus seg1 and seg2
        if (!(dgf_num / 2 == 15 && dgf_ch == 0) &&
            !(dgf_num / 2 == 17 &&
              dgf_ch == 0)) {  // pause filling crystal 15 core gammas
          E_gam_seg[dgf_num / 6][dgf_num % 6 / 2][dgf_ch]->Fill(dgf_en);
          E_gam_seg_cal[dgf_num / 6][dgf_num % 6 / 2][dgf_ch]->Fill(
              GammaEnergy);
        }
      }  // even DGF number

      else {
        E_gam_seg[dgf_num / 6][dgf_num % 6 / 2][dgf_ch + 3]->Fill(dgf_en);
        E_gam_seg_cal[dgf_num / 6][dgf_num % 6 / 2][dgf_ch + 3]->Fill(
            GammaEnergy);
      }  // odd DGF number

      // Look for a core event
      if (dgf_ch == 0 && dgf_num % 2 == 0) {
        MaxSegClu = dgf_num / 6;
        MaxSegCore = dgf_num % 6 / 2;
        if (dead_segments[MaxSegClu][MaxSegCore].size() == 7) continue;

        MaxSegId =
            0;  // initialise as core (if no segment hit (dead), use core!)
        SegSumEnergy = 0;  // add segment energies
        c15SegSumDgf_s12 = 0;
        c15SegSumDgf_s3456 = 0;
        c17_s4 = 0;
        c17_s5 = 0;
        c17_s45 = 0;
        c17_s1236 = 0;

        seg_mul = 0;  // segment multiplicity

        // Check for highest energy segment in same detector
        for (unsigned int k = 0; k < event->NumberOfDgfs(); k++) {
          dgf_num2 = event->Dgf(k)->ModuleNumber();
          dgf_ch2 = event->Dgf(k)->Channel();
          dgf_en2 = event->Dgf(k)->Energy();
          dgf_t2 = event->Dgf(k)->Time() + Cal->CryTime(dgf_num2 / 2);

          // We don't care if it's the same event, i.e. core event
          if (k == j) continue;

          // Make sure it's a Miniball DGF
          if (dgf_num2 < 0 || dgf_num2 > 47) continue;

          // Get global segment number (0-167)
          gSeg = (dgf_num2 / 2) * 7;
          if (dgf_num2 % 2 == 0)
            gSeg += dgf_ch2;
          else
            gSeg += dgf_ch2 + 3;

          // Skip if a different detector
          if (dgf_num2 != dgf_num && dgf_num2 != dgf_num + 1) continue;

          // Is it a crap segment?
          // for( unsigned int ds = 0; ds < dead_segments.size(); ds++ ) {

          //   if( gSeg == dead_segments[ds] && dgf_en2 > 20 ) veto_gamma =
          //   true;

          // }

          // Check if it should be swapped
          for (unsigned int ss = 0; ss < swap_segments.size(); ss++) {
            if (gSeg == swap_segments[ss][0]) gSeg = swap_segments[ss][1];
            if (gSeg == swap_segments[ss][1]) gSeg = swap_segments[ss][0];
          }

          GammaEnergy2 = Cal->DgfEnergy(dgf_num2, dgf_ch2, dgf_en2);

          // Calculate sum of segments if energy is deposited (i.e. >1 keV)
          if (GammaEnergy2 > 1.) {
            SegSumEnergy += GammaEnergy2;
            seg_mul++;
          }
          if (dgf_num2 / 2 == 15) {  // crystal 15, needs crosstalk correction
            // segments 1-2 vs 3-6 correspond to different eCore
            if (dgf_num2 % 2 == 0)
              c15SegSumDgf_s12 += dgf_en2;
            else
              c15SegSumDgf_s3456 += dgf_en2;
          } else if (dgf_num2 / 2 ==
                     17) {  // crystal 17, needs crosstalk correction
            if (dgf_num2 % 2 == 1) {
              if (dgf_ch2 == 1) {
                c17_s4 = dgf_en2;
                c17_s45 += dgf_en2;
              } else if (dgf_ch2 == 2) {
                c17_s5 = dgf_en2;
                c17_s45 += dgf_en2;
              } else
                c17_s1236 += dgf_en2;
            } else
              c17_s1236 += dgf_en2;
          }
          // Test maximum energy segment
          if (GammaEnergy2 < MaxSegEnergy) continue;

          // Reassign energy and id
          MaxSegEnergy = GammaEnergy2;
          MaxSegId = gSeg % 7;

        }  // k

        // Apply c15 correction
        if (c15SegSumDgf_s12 > 0 || c15SegSumDgf_s3456 > 0) {
          double fc15 = C15XTALK * ((c15SegSumDgf_s12 - c15SegSumDgf_s3456) /
                                    (c15SegSumDgf_s12 + c15SegSumDgf_s3456));
          //	  cout<<c15SegSumDgf_s12<<" "<<c15SegSumDgf_s3456<<"
          //"<<fc15<<endl; 	  cout<<"previous c15 dgf_en: "<<dgf_en<<endl;

          dgf_en = dgf_en / (1. + fc15);
          //	  cout<<"new c15 dgf_en: "<<dgf_en<<endl;

          E_gam_seg[5][0][0]->Fill(dgf_en);
          GammaEnergy = Cal->DgfEnergy(30, 0, dgf_en);
          E_gam_seg_cal[5][0][0]->Fill(GammaEnergy);

        } else if (MaxSegClu == 5 && MaxSegCore == 0) {
          E_gam_seg[5][0][0]->Fill(dgf_en);
          GammaEnergy = Cal->DgfEnergy(30, 0, dgf_en);
          E_gam_seg_cal[5][0][0]->Fill(GammaEnergy);
        }
        // Apply c17 cuts and corrections
        if (MaxSegClu == 5 && MaxSegCore == 2 && SegSumEnergy == 0)
          continue;  // skip core-only events

        if (c17_s4 > 0 || c17_s5 > 0) {
          if (c17_s4 > 0 && c17_s5 == 0) {
            double frac = (double)(c17_s4 + c17_s1236) / dgf_en;
            if (frac > C17S4CUT)
              dgf_en = dgf_en * (1 + C17S4XTALK * (1. - frac));

            else
              dgf_en =
                  dgf_en * (1 + C17S4XTALK * ((1. - C17S4CUT) +
                                              (frac - C17S4CUT) / C17S4CUT *
                                                  (1. - C17S4CUT)));

          } else if (c17_s5 > 0 && c17_s4 == 0) {
            double frac = (double)(c17_s5 + c17_s1236) / dgf_en;
            if (frac > C17S5CUT)
              dgf_en = dgf_en * (1 + C17S5XTALK * (1. - frac));
            else
              dgf_en =
                  dgf_en * (1 + C17S5XTALK * ((1. - C17S5CUT) +
                                              (frac - C17S5CUT) / C17S5CUT *
                                                  (1. - C17S5CUT)));
          } else {
            double frac = ((double)c17_s45 + c17_s1236) / dgf_en;
            dgf_en = dgf_en * (1 + C17S45XTALK * (1. - frac));
          }
          E_gam_seg[5][2][0]->Fill(dgf_en);
          GammaEnergy = Cal->DgfEnergy(34, 0, dgf_en);
          E_gam_seg_cal[5][2][0]->Fill(GammaEnergy);
        }  // segments 4 or 5 hit
        else if (MaxSegClu == 5 && MaxSegCore == 2) {
          E_gam_seg[5][2][0]->Fill(dgf_en);
          GammaEnergy = Cal->DgfEnergy(34, 0, dgf_en);
          E_gam_seg_cal[5][2][0]->Fill(GammaEnergy);
        }
        // Found highest energy segment //

        // Do the veto of crap segments
        //	if( veto_gamma ){
        //  cout<<"Check veto_gamma"<<endl;
        //   continue;
        // }
        //	cout<<"check here"<<endl;

        // check for reassigning dead-segment events
        if (dead_segments[MaxSegClu][MaxSegCore].size() ==
            1) {  // one dead segment, check if missing energy is greater
          if (MaxSegEnergy <
              GammaEnergy /
                  (double)(seg_mul +
                           1.)) {  // more energy in missing segment suspected;
                                   // check 40%(dead)/30%/30% distributions
            MaxSegEnergy = GammaEnergy - SegSumEnergy;
            MaxSegId = dead_segments[MaxSegClu][MaxSegCore][0];
          }
          SegSumEnergy = GammaEnergy;  // energy "recovered" anyhow
        } else if (dead_segments[MaxSegClu][MaxSegCore].size() >=
                   2) {  // ambiguity in dead segments, check energy
          if (MaxSegEnergy <
              GammaEnergy /
                  (double)(seg_mul +
                           1.)) {  // more energy in missing segments suspected
            MaxSegId = 0;          // resort to core as max. energy segment;
          }
          SegSumEnergy = GammaEnergy;
        }

        // Compare segment sum energy and core energy - less than 0.5%
        // difference or 2 keV
        if ((TMath::Abs(SegSumEnergy - GammaEnergy) / GammaEnergy < 0.005 ||
             TMath::Abs(SegSumEnergy - GammaEnergy) < 2.) &&
            segsum)
          GammaEnergy = SegSumEnergy;  // Overwrite with segment sum energy
        // else
        //   continue; // really throw it out

        // If energy is outside of 1% window and not using reject, throw it away
        // (comment out the line below if you want to use core energy instead)
        //	else if( segsum && !reject ) continue; // assuming energy
        //imbalance due to dead segment and not random signal

        // If using reject and segsum in combination, take only single-segment
        // hits
        if (segsum && reject) {
          if (seg_mul == 1)
            GammaEnergy = SegSumEnergy;
          else
            continue;
        }

        // Check <previous> gammas for addback (for a given crystal)
        for (unsigned int l = 0; l < gen_array.size(); l++) {
          if (gen_array[l] < GammaEnergy + 0.2 &&
              gen_array[l] > GammaEnergy - 0.2 && sid_array[l] == MaxSegId)
            continue;  // same event

          // Do the addback with the 3 conditions:
          // for the given gamma ray with same cluster
          // addback flag is on
          // ab_array[l] flag is off (not added back)
          if (clu_array[l] == MaxSegClu && addback && !ab_array[l]) {
            gen_array[l] += GammaEnergy;  // apply energy addback
            ab_evt = true;
            ab_mul++;
            ab_array.at(l) = true;

            if (sen_array[l] <
                MaxSegEnergy) {  // update with higher-energy segment info

              gtd_array[l] = dgf_t;
              clu_array[l] = MaxSegClu;                   // cluster number
              cid_array[l] = MaxSegClu * 3 + MaxSegCore;  // crystal number
              sid_array[l] = MaxSegId;                    // segment number
              sen_array[l] = MaxSegEnergy;
            }

          }  // addback

          // Reject/suppress if same cluster (if not doing single crystal
          // reject)
          if (clu_array[l] == MaxSegClu && reject && !segsum) {
            gen_array.erase(gen_array.begin() + l);
            gtd_array.erase(gtd_array.begin() + l);
            clu_array.erase(clu_array.begin() + l);
            cid_array.erase(cid_array.begin() + l);
            sid_array.erase(sid_array.begin() + l);
            sen_array.erase(sen_array.begin() + l);

            reject_evt = true;
            l--;

          }  // reject

        }  // previous gammas' loop finished

        if (ab_evt)
          continue;  // get next gamma (this one (j) has been added back to the
                     // previous gamma with same cluster)
        if (reject_evt) continue;  // get next gamma (this one is rejected)
        hABmult->Fill(ab_mul);

        gen_array.push_back(GammaEnergy);
        gtd_array.push_back(dgf_t);
        clu_array.push_back(MaxSegClu);
        cid_array.push_back(MaxSegClu * 3 + MaxSegCore);
        sid_array.push_back(MaxSegId);
        sen_array.push_back(MaxSegEnergy);
        if (ab_mul == 0)
          ab_array.push_back(false);
        else
          ab_array.push_back(true);

      }  // core event

      // Reset
      MaxSegEnergy = -99.;
      MaxSegClu = -1;
      MaxSegCore = -1;
      MaxSegId = -1;
    }
    // beam dump gamma energies
    if (dgf_num == 53 && dgf_ch >= 0 && dgf_ch < 2) {
      GammaEnergy = Cal->DgfEnergy(dgf_num, dgf_ch, dgf_en);
      E_BeamDump[dgf_ch]->Fill(GammaEnergy);
      T_BeamDump[dgf_ch]->Fill((dgf_t) / 40000000);

      for (unsigned int k = 0; k < event->NumberOfDgfs(); k++) {
        if (k == j) continue;
        dgf_num2 = event->Dgf(k)->ModuleNumber();
        dgf_ch2 = event->Dgf(k)->Channel();
        dgf_en2 = event->Dgf(k)->Energy();
        dgf_t2 = event->Dgf(k)->Time() + Cal->CryTime(dgf_num2 / 2);

        if (dgf_num2 != 53) continue;
        if (dgf_ch != 0 || dgf_ch2 != 1) continue;

        GammaEnergy2 = Cal->DgfEnergy(dgf_num2, dgf_ch2, dgf_en2);

        tdiff_BD->Fill(dgf_t - dgf_t2);

        if (TMath::Abs(dgf_t - dgf_t2) < 999999999.)

          bd_bd->Fill(GammaEnergy, GammaEnergy2);
      }
    }
  }

  // Fill the gamma-ray singles histogram
  for (unsigned int i = 0; i < gen_array.size(); i++) {
    E_gam_tot->Fill(gen_array[i]);

    E_gam_vs_seg->Fill(cid_array[i] * 7 + sid_array[i],
                       gen_array[i]);                 // segment
    E_gam_vs_core->Fill(cid_array[i], gen_array[i]);  // core
  }

  return;
}

void AddBack::MakeElectrons() {
  // SPEDE is in the 5th ADC
  for (unsigned int i = 0; i < subevent->Size(); i++) {
    adc_ch = subevent->AdcChannel(i);
    adc_en = subevent->AdcValue(i);

    ElectronEnergy = Cal->AdcEnergy(4, adc_ch, adc_en);

    // STM-16 one
    if (adc_ch < 12) {
      E_spede_seg[adc_ch]->Fill(adc_en);
      E_spede_seg_cal[adc_ch]->Fill(ElectronEnergy);
      E_spede->Fill(ElectronEnergy);

      gen_array.push_back(ElectronEnergy);
      gtd_array.push_back(adc_t);
      clu_array.push_back(8);
      cid_array.push_back(0);
      sid_array.push_back(adc_ch);
      sen_array.push_back(0);

    }

    // STM-16 two
    else if (adc_ch > 15 && adc_ch < 28) {
      E_spede_seg[adc_ch - 4]->Fill(adc_en);
      E_spede_seg_cal[adc_ch - 4]->Fill(ElectronEnergy);
      E_spede->Fill(ElectronEnergy);

      gen_array.push_back(ElectronEnergy);
      gtd_array.push_back(adc_t);
      clu_array.push_back(8);
      cid_array.push_back(0);
      sid_array.push_back(adc_ch - 4);
      sen_array.push_back(0);
    }

  }  // k

  return;
}

#endif
