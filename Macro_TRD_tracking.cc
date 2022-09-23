
//#include "TRD_TFTracker.h"
//#include "TRD_TFTracker.cxx"

R__LOAD_LIBRARY(TRD_Kalman_Tracking_cxx.so);
R__LOAD_LIBRARY(TRD_ST_Analyze_tracklets_cxx.so);

// Environment variables
//#define ENV_PI
#define ENV_ALEX
//#define ENV_PI_SVEN

void Macro_TRD_tracking(TString input_list = "run0_test.txt", Int_t event_plot = -1, Double_t track_path = 1000.0, Double_t beam_path = 0.0)
{
    // event_plot: -1 -> loop over all events, other wise plot or loop over single event

    // First compile type root and then:
    // .L TRD_ST_Analyze_tracklets.cxx++
    // .L TRD_Kalman_Tracking.cxx++
    // .x Macro_TRD_tracking.cc("Split_tracklets_vD_1.546_LA_0.16133_V3_3-3.txt",-1,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("run0_test.txt",-1,1000.0,0.0) -> calibrated run0 files

    // .x Macro_TRD_tracking.cc("Split_tracklets_calib_A_V0_1681-1700.txt",-1,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("new_trkl_test.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_PbPb2018_test.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_pPb2016_test.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_ppSim_test.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_Martin.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_test_new_align.txt",0,1000.0,0.0)

    // Run3 test
    // .x Macro_TRD_tracking.cc("List_Run3_test.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_Run3_test_cosmics.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_Run3_testB_cosmics.txt",0,1000.0,0.0)

    // Pb-Pb
    // .x Macro_TRD_tracking.cc("List_PbPb2018_testB.txt",0,1000.0,0.0)

    // pp 900 GeV
    // .x Macro_TRD_tracking.cc("List_pp900GeV.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_pp900GeV_corr.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_pp900GeV_corr_B02.txt",0,1000.0,0.0)

    // pp 900 GeV B = 0.2T
    // .x Macro_TRD_tracking.cc("List_00505673.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_test_00505673.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_single_137_00505673_V3.txt",0,1000.0,0.0)
    // .x Macro_TRD_tracking.cc("List_single_137TS_00505673_V3.txt",0,1000.0,0.0)


    // .x Macro_TRD_tracking.cc("List_S_particle.txt",0,1000.0,0.0)

    // 2069, 1531
    // opposite 1507


    gROOT->SetStyle("Plain");
    gStyle->SetOptFit(11);
    //gStyle->SetOptFit(1111);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFillColor(0);
    gStyle->SetPalette(27);

    gStyle->SetPalette( kDarkBodyRadiator);



    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasBorderMode(0);


    //gSystem ->Load("TRD_Kalman_Tracking_cxx.so");
    gSystem ->Load("TRD_ST_Analyze_tracklets_cxx.so");

    //------------------------------------
    // Define output file name and directory
    TString out_file_name = input_list;
    out_file_name += "_out.root";

    TString out_file_name_calib = input_list;
    out_file_name_calib += "_out_calib.root";

    TString input_dir  = "./Data_bkg/";
    TString output_dir = "./";
    TString inlists_dir = "./Lists_tracklets/";
#if defined(ENV_PI)
    inlists_dir = "/home/ceres/schmah/ALICE/TRD_self_tracking/Lists_tracklets/";
    input_dir   = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/Calib_tracklets/";
    output_dir  = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/ST_out/";
#endif

#if defined(ENV_PI_SVEN)
    inlists_dir = "/home/ceres/schmah/ALICE/TRD_self_tracking/Lists_tracklets/";
    input_dir   = "/misc/alidata120/alice_u/schmah/TRD_self_tracking/Calib_tracklets/";
    output_dir  = "/misc/alidata120/alice_u/hoppner/TRD_self_tracking/ST_out/";
#endif

#if defined(ENV_ALEX)
    inlists_dir = "./Lists_tracklets/";
    input_dir  = "./Data_bkg/";
    output_dir = "./ST_out/";
#endif
    Int_t use_prim_vertex           = -1; // 0 = no primary vertex, 1 = primary vertex used
    Int_t KF_tracker                = 0; // Kalman filter tracker
    Int_t TF_tracker                = 0; // Tensorflow tracker

    Int_t graphics                  = 0; // 0 = no 3D graphics, 1 = 3D graphics (#define USEEVE in TRD_ST_Analyze_tracklets needs to be defined too)
    Int_t draw_tracklets_TPC_match  = 0; // Draw tracklets matched with TPC tracks
    Int_t draw_all_TPC_tracks       = 1; // Draw all TPC tracks
    Int_t draw_all_TRD_tracks       = 0; // Draw all TRD tracks ->
    Int_t draw_all_tracklets        = 1; // Draw all TRD tracklets
    Int_t draw_found_tracklets      = 1; // Draws tracklets found by tracker
    Int_t draw_matched_TPC_track    = 0; // Draw TPC to TRD matched TPC track
    Int_t draw_matched_TRD_track    = 0; // Draw TPC to TRD matched Kalman/TF track
    Int_t draw_secondary_vertices   = 0; // Draws tracks and secondary vertices
    Int_t calibrate_gain            = 0; // Calibrate TRD chamber gain
    Int_t Bethe_flag                = 1; // 1: Bethe Bloch for gain calib || 0: TPC dEdx for gain calib
    Int_t animate_beams             = 0; // 1 for beam animation
    Int_t calibrate_vD              = 0;
    Int_t is_MC_event               = 1;
    Int_t draw_MC_tracks            = 1; // Draw Monte Carlo event
    Int_t draw_TRD_digits           = 0; // Draw TRD digits (old versions don't have the tree)
    Int_t bool_make_invariant_mass_hist_up_to_layer = 0; //calculate invariant mass of S particle up to generation <Int_t> and draws a histogram
    Int_t calc_reconstuction_efficiency_up_to_layer = 0;
    Bool_t draw_only_S_quark = 0;
    Bool_t draw_V0_radius_Lambda = 1;
    Bool_t calc_second_vertices = 0;
    //------------------------------------

    printf("TRD_ST_Analyze_tracklets started \n");
    Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze(output_dir,out_file_name,graphics);
    //Ali_TRD_ST_Analyze*  TRD_ST_Analyze = new Ali_TRD_ST_Analyze(output_dir,"test.root",graphics);

    TRD_ST_Analyze ->set_input_lists(inlists_dir);
    TRD_ST_Analyze ->set_input_dir(input_dir);
    TRD_ST_Analyze ->Init_tree(input_list.Data());

    Long64_t N_Events = TRD_ST_Analyze ->get_N_Events();
    TH1D* h_layer_radii = TRD_ST_Analyze ->get_layer_radii_hist();
    //Long64_t event = 10;

    TRD_Kalman_Trackfinder kalid;
    kalid.set_layer_radii_hist(h_layer_radii);

    //------------------------------------------------------------
    //Tensorflow tracker
  //  TRD_TFTrackMaker tftracker;
    //------------------------------------------------------------



    // photon events: 88, 285, 29, 273
    // photon events thermal shield: 47,
    // pi0 event: 378
    // nuclear interaction event: 158, 168(!), 3741, 92, 328(!)
    Int_t start_event = 0;
    Int_t stop_event  = (Int_t) N_Events;
    //Int_t stop_event  = 1;
    if(event_plot != -1)
    {
        start_event = event_plot;
        stop_event  = (event_plot+1);
    }


    for(Long64_t event = start_event; event < stop_event; event++)
    {

        if (event != 0  &&  event % 50 == 0)
            cout << "." << flush;
        if (event != 0  &&  event % 500 == 0)
        {
            printf("event: %lld out of %lld, %4.2f%% total done \n",event,N_Events,((Double_t)event/(Double_t)N_Events)*100.0);
        }

        vector< vector<Ali_TRD_ST_Tracklets*> > tracker_found_tracklets_KF;
        vector<vector<Double_t>>                mHelices_tracker_KF;
        vector<Double_t>                        mChi_2s_tracker_KF;
        vector< vector<Ali_TRD_ST_Tracklets*> > tracker_found_tracklets_TF;
        vector<vector<Double_t>>                mHelices_tracker_TF;

        TRD_ST_Analyze ->Loop_event(event,graphics);
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;
        
        //printf("new MC_event %d\n----------------------------------------\n",event);
        
        if(graphics) TRD_ST_Analyze ->Draw_event(event,graphics,draw_all_TPC_tracks,draw_all_tracklets,track_path,draw_TRD_digits);  // ->draws TPC tracks
        if(is_MC_event) TRD_ST_Analyze ->Scan_MC_Event(graphics,bool_make_invariant_mass_hist_up_to_layer);
        if(calc_reconstuction_efficiency_up_to_layer||draw_V0_radius_Lambda) TRD_ST_Analyze ->Check_possible_reconstruction(calc_reconstuction_efficiency_up_to_layer,event,draw_V0_radius_Lambda);
        if(graphics && draw_MC_tracks) TRD_ST_Analyze ->Draw_MC_event(event,graphics,draw_only_S_quark);
        //cout<<TRD_ST_Analyze->Tracklets[2]->get_TRD_index()<<endl;
        if(graphics && animate_beams) TRD_ST_Analyze ->Animate_beams(beam_path);

        if(draw_tracklets_TPC_match || KF_tracker) TRD_ST_Analyze ->Do_TPC_TRD_matching(event,3.0,10.0,graphics*draw_tracklets_TPC_match); // last one is graphics  --> draws kalman TRD tracklets

        if (calibrate_gain) TProfile* gain = TRD_ST_Analyze->Calibrate_gain(event,Bethe_flag);

        TRD_ST_Analyze ->set_TPC_helix_params(event);


        //------------------------------------------------------------
        //Kalman Filter tracker
        if(KF_tracker)
        {
            tracker_found_tracklets_KF = kalid.Kalman_Trackfind(TRD_ST_Analyze->Tracklets,TRD_ST_Analyze->Number_Tracklets,use_prim_vertex); // 0 = no primary vertex, 1 = primary vertex used
            vector< vector<Ali_TRD_ST_Tracklets*> > matched_tracks=TRD_ST_Analyze->matched_tracks; // TPC track matched tracklets
            vector< vector<Ali_TRD_ST_Tracklets*> > matched_beautiful_tracks;
            mHelices_tracker_KF = kalid.get_Kalman_helix_params();
            mChi_2s_tracker_KF = kalid.get_Kalman_chi_2();

            //if((Int_t)mHelices_tracker_KF.size() > 0) printf("event: %lld, size of mHelices_tracker: %d \n",event,(Int_t)mHelices_tracker_KF.size());

            //add here chi2 or use different new class already? Ali_Helix replace with new?

            TRD_ST_Analyze ->set_Kalman_helix_params(mHelices_tracker_KF);
            TRD_ST_Analyze ->set_Kalman_chi_2(mChi_2s_tracker_KF);
            TRD_ST_Analyze ->set_Kalman_TRD_tracklets(tracker_found_tracklets_KF);
            TRD_ST_Analyze ->flag_TRD_tracks_with_shared_tracklets();

            if(graphics && draw_found_tracklets) TRD_ST_Analyze ->Draw_Kalman_Tracklets(tracker_found_tracklets_KF); // Draws the Kalman matched TRD tracklets
            if(graphics && draw_all_TRD_tracks)  TRD_ST_Analyze ->Draw_Kalman_Helix_Tracks(-1,kGreen,5.0,500.0); // -1 -> all kalman tracks drawn, 3.0, 500.0
            TRD_ST_Analyze ->Match_kalman_tracks_to_TPC_tracks(graphics,draw_matched_TPC_track,draw_matched_TRD_track,kGreen+1);
            TRD_ST_Analyze ->Calc_Kalman_efficiency();

        }
        //------------------------------------------------------------


        if(calibrate_vD) TRD_ST_Analyze ->Calibrate(graphics);
        
        TRD_ST_Analyze ->set_self_event_info();
        if(calc_second_vertices){
            Int_t found_good_AP_vertex_TPC = TRD_ST_Analyze ->Calculate_secondary_vertices(graphics*draw_secondary_vertices,0,0); // (0 = no graphics), (0 = TRD, 1 = TPC), flag_fill_tree (0 don't fill)
            Int_t found_good_AP_vertex_TRD = TRD_ST_Analyze ->Calculate_secondary_vertices(graphics*draw_secondary_vertices,1,1); // (0 = no graphics), (0 = TRD, 1 = TPC), flag_fill_tree (0 don't fill)
        }//if(found_good_AP_vertex) printf(" ----> Good AP vertex found in event: %lld \n",event);

        //vector< vector<Ali_TRD_ST_Tracklets*> > tracker_found_tracklets=kalid.found_tracks;
        //vector<Double_t> track_accuracy;

        //Ali_TRD_ST_Tracklets* last_tracklet;
        TRD_ST_Analyze ->Write_Sexaquark_data_to_Tree(event);
    
    }
    TRD_ST_Analyze -> find_V0_pairs();
    if(bool_make_invariant_mass_hist_up_to_layer)   TRD_ST_Analyze ->Draw_Inv_Mass_histogram();
    if(draw_V0_radius_Lambda)                       TRD_ST_Analyze ->Draw_Count_vs_V0Lambda_histogram();
    if(calc_reconstuction_efficiency_up_to_layer)   TRD_ST_Analyze ->Count_reconstruction_efficiency();
    if(calibrate_vD)                                TRD_ST_Analyze ->Draw_n_Save_Calibration(output_dir,out_file_name_calib);
    if(calibrate_gain)                              TRD_ST_Analyze ->Draw_Save_Gain_calib();

    TRD_ST_Analyze ->Write();

}
