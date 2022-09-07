// vim:set ts=4 sw=4 fdm=marker et:
#include "Bfinder/Bfinder/interface/format.h"
#include "Bfinder/Bfinder/interface/Dntuple.h"
#include "Bfinder/Bfinder/interface/utilities.h"
#include "boost/tuple/tuple.hpp"
#include <chrono>
//
// class declaration
//
////add some class here.




class Dfinder : public edm::EDAnalyzer
{//{{
//{
    public:
        explicit Dfinder(const edm::ParameterSet&);
        ~Dfinder();
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		///here add some type.
		struct Track {
			double pt, eta, phi;
			int q;
			int index;

			Track(double ptp, double etap, double phip, int qp, int indexp) : pt(ptp), eta(etap), phi(phip), q(qp), index(indexp) { }

		};

		struct TrackXYZ {
			double px, py, pz;
			int q;
			int index;

			TrackXYZ(double pt, double eta, double phi, int qp, int indexp) : q(qp), index(indexp) { setup(pt, eta, phi); }

			TrackXYZ(const Track& tr) : q(tr.q) , index(tr.index) { setup(tr.pt, tr.eta, tr.phi); }

			void setup(double pt, double eta, double phi) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
			}

		};//TrackXYZ

		struct TrackXYZP2 {
			double px, py, pz, p2;
			int q;
			int index;

			TrackXYZP2(double pt, double eta, double phi, int qp, int indexp) : q(qp),index(indexp) { setup(pt, eta, phi); }

			TrackXYZP2(const Track& tr) : q(tr.q), index(tr.index) { setup(tr.pt, tr.eta, tr.phi); }

			TrackXYZP2(const TrackXYZ& tr) : px(tr.px), py(tr.py), pz(tr.pz), q(tr.q), index(tr.index) {
				p2 = px * px + py * py + pz * pz;
			}

			void setup(double pt, double eta, double phi) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
				p2 = pt * pt + pz * pz;
			}

		};//TrackXYZP2

		struct Triplet {
			int i1, i2, i3;
			//double m2; // invariant mass squared (can reconstruct p0 and y from it later)

			Triplet(int i1p, int i2p, int i3p) : i1(i1p), i2(i2p), i3(i3p) { }
		};



		struct TrackSelected{
			 int index_tk1, index_tk2, index_tk3, permutation_number;

			 TrackSelected(int index_tk1p, int index_tk2p, int index_tk3p, int permutation_numberp) : index_tk1(index_tk1p), index_tk2(index_tk2p), index_tk3(index_tk3p), permutation_number(permutation_numberp) { }
		};//TrackXYZ


		struct P3 {
			double px, py, pz;

			P3(double pt, double eta, double phi) {
				setup(pt, eta, phi);
			}

			P3(const Track& tr) {
				setup(tr.pt, tr.eta, tr.phi);
			}

			P3(const TrackXYZ& tr) {
				px = tr.px;
				py = tr.py;
				pz = tr.pz;
			}

			P3(const TrackXYZP2& tr) {
				px = tr.px;
				py = tr.py;
				pz = tr.pz;
			}

			void setup(double pt, double eta, double phi) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
			}

			P3& operator += (const P3& a) {
				px += a.px;
				py += a.py;
				pz += a.pz;
				return *this;
			}

		};//P3




		struct P4 {
			double px, py, pz, p0;

			P4(double pt, double eta, double phi, double m) {
				setup(pt, eta, phi, m);
			}

			P4(const Track& tr, double m) {
				setup(tr.pt, tr.eta, tr.phi, m);
			}

			P4(const TrackXYZ& tr, double m) {
				px = tr.px;
				py = tr.py;
				pz = tr.pz;
				double p2 = px * px + py * py + pz * pz;
				p0 = sqrt(p2 + m * m);
			}


			void setup(double pt, double eta, double phi, double m) {
				px = pt * cos(phi);
				py = pt * sin(phi);
				pz = pt * sinh(eta);
				double p2 = pt * pt + pz * pz;
				p0 = sqrt(p2 + m * m);
			}

			P4& operator += (const P4& a) {
				px += a.px;
				py += a.py;
				pz += a.pz;
				p0 += a.p0;
				return *this;
			}

		};//P4



    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
 
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        
        virtual std::vector< std::vector< std::pair<float, int> > > GetPermu(std::vector< std::pair<float, int> > InVec);
        virtual std::vector< std::vector< std::pair<float, int> > > DelDuplicate(std::vector< std::vector< std::pair<float, int> > > InVec);
        
		virtual void BranchOutNTk(
            DInfoBranches &DInfo, 
            std::vector<pat::GenericParticle> input_tracks, 
            reco::Vertex thePrimaryV,
			reco::Vertex theBeamSpotV,	
			vector<Track> lst,
			vector<TrackXYZP2> lstXYZP2,
            std::vector<int> &D_counter,
            float *mass_window,
            double tktkRes_mass,
            double tktkRes_mass_window,
            bool doConstrainFit,
            bool SequentialFit,
            int Dchannel_number,
			int TkCombinationMethod = -1
        );
		
		virtual void TkCombinationPermutation_Lc_v3(
            std::vector<pat::GenericParticle> input_tracks, 
            //std::vector<int> isNeededTrackIdx,
			vector<Track> lst,
			vector<TrackXYZP2> lstXYZP2,
            float *mass_window,
            //std::vector< std::pair<float, int> > TkMassCharge,
            double tktkRes_mass,
            double tktkRes_mass_window,
        	vector<TrackSelected> &selectedTkhidxSet,
			int Dchannel_number
        );
        // ----------member data ---------------------------
		edm::ESHandle<MagneticField> bField;
        edm::ParameterSet theConfig;

        bool detailMode_;
        bool dropUnusedTracks_;
        std::vector<int> Dchannel_;
        //edm::InputTag hltLabel_;
        edm::EDGetTokenT< reco::GenParticleCollection > genLabel_;
        edm::EDGetTokenT< std::vector<pat::GenericParticle> > trackLabel_;
        edm::EDGetTokenT< std::vector<PileupSummaryInfo> > puInfoLabel_;
        edm::EDGetTokenT< reco::BeamSpot > bsLabel_;
        edm::EDGetTokenT< reco::VertexCollection > pvLabel_;
        edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
        edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
        
        double tkPtCut_;
        double tkEtaCut_;
        std::vector<double> dCutSeparating_PtVal_;
        std::vector<double> dPtCut_;
        std::vector<double> dEtaCut_;
        std::vector<double> dRapidityCut_;
        std::vector<double> VtxChiProbCut_;
        std::vector<double> svpvDistanceCut_lowptD_;
        std::vector<double> svpvDistanceCut_highptD_;
        std::vector<double> alphaCut_;
        std::vector<double> MaxDocaCut_;
        std::vector<double> tktkRes_dCutSeparating_PtVal_;
        std::vector<double> tktkRes_dPtCut_;
        std::vector<double> tktkRes_dEtaCut_;
        std::vector<double> tktkRes_VtxChiProbCut_;
        std::vector<double> tktkRes_svpvDistanceCut_lowptD_;
        std::vector<double> tktkRes_svpvDistanceCut_highptD_;
		std::vector<double> tktkRes_alphaCut_;
        bool RunOnMC_;
        bool doTkPreCut_;
        bool makeDntuple_;
        bool doDntupleSkim_;
        bool printInfo_;
		bool printFill_Info_;
        bool readDedx_;
        edm::EDGetTokenT<edm::ValueMap<float> > MVAMapLabel_;
        edm::EDGetTokenT< std::vector<float> > MVAMapLabelpA_;
        edm::InputTag MVAMapLabelInputTag_;

        edm::Service<TFileService> fs;
        TTree *root;
        EvtInfoBranches     EvtInfo;
        VtxInfoBranches     VtxInfo;
        TrackInfoBranches   TrackInfo;
        DInfoBranches       DInfo;
        GenInfoBranches     GenInfo;
        CommonFuncts        Functs;
        DntupleBranches     *Dntuple = new DntupleBranches;
        TTree* ntD1; 
        TTree* ntD2;
        TTree* ntD3; 
        TTree* ntD4; 
        TTree* ntD5; 
        TTree* ntD6; 
        TTree* ntD7; 
        TTree* ntD8;
        TTree* ntGen;

        //histograms
        TH1F *TrackCutLevel;
        //How many channel
        static int const Nchannel = 20;
        std::vector<TH1D*> DMassCutLevel;
	
        
};//}}}

void Dfinder::beginJob()
{//{{{
    //auto t0 = std::chrono::high_resolution_clock::now();
    
	root = fs->make<TTree>("root","root");
    ntD1 = fs->make<TTree>("ntDkpi","");           Dntuple->buildDBranch(ntD1);
    ntD2 = fs->make<TTree>("ntDkpipi","");         Dntuple->buildDBranch(ntD2);
    ntD3 = fs->make<TTree>("ntDkpipipi","");       Dntuple->buildDBranch(ntD3);
    ntD4 = fs->make<TTree>("ntDPhikkpi","");       Dntuple->buildDBranch(ntD4);
    ntD5 = fs->make<TTree>("ntDD0kpipi","");       Dntuple->buildDBranch(ntD5);
    ntD6 = fs->make<TTree>("ntDD0kpipipipi","");   Dntuple->buildDBranch(ntD6);
    ntD7 = fs->make<TTree>("ntBptoD0pi","");       Dntuple->buildDBranch(ntD7);
    ntD8 = fs->make<TTree>("ntLambdaCtopkpi","");  Dntuple->buildDBranch(ntD8);
    ntGen = fs->make<TTree>("ntGen","");           Dntuple->buildGenBranch(ntGen);
    EvtInfo.regTree(root);
    VtxInfo.regTree(root);
    TrackInfo.regTree(root, detailMode_);
    DInfo.regTree(root, detailMode_);
    GenInfo.regTree(root);
/*	
	auto t1 = std::chrono::high_resolution_clock::now(); 
	double dt = 1e-3 * std::chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
	std::cout<<"duration: "<<dt<<std::endl;
*/
}//}}}

Dfinder::Dfinder(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
    //auto t0 = std::chrono::high_resolution_clock::now();
	detailMode_ = iConfig.getParameter<bool>("detailMode");
    dropUnusedTracks_ = iConfig.getParameter<bool>("dropUnusedTracks");

    Dchannel_ = iConfig.getParameter<std::vector<int> >("Dchannel");
    genLabel_           = consumes< reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>("GenLabel"));
    trackLabel_         = consumes< std::vector<pat::GenericParticle> >(iConfig.getParameter<edm::InputTag>("TrackLabel"));
    //hltLabel_           = iConfig.getParameter<edm::InputTag>("HLTLabel");
    puInfoLabel_    = consumes< std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfoLabel"));
    bsLabel_        = consumes< reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("BSLabel"));
    pvLabel_        = consumes< reco::VertexCollection >(iConfig.getParameter<edm::InputTag>("PVLabel"));

    tkPtCut_ = iConfig.getParameter<double>("tkPtCut");
    tkEtaCut_ = iConfig.getParameter<double>("tkEtaCut");
    dCutSeparating_PtVal_ = iConfig.getParameter<std::vector<double> >("dCutSeparating_PtVal");
    dPtCut_ = iConfig.getParameter<std::vector<double> >("dPtCut");
    dEtaCut_ = iConfig.getParameter<std::vector<double> >("dEtaCut");
    dRapidityCut_ = iConfig.getParameter<std::vector<double> >("dRapidityCut");
    VtxChiProbCut_ = iConfig.getParameter<std::vector<double> >("VtxChiProbCut");
    svpvDistanceCut_lowptD_ = iConfig.getParameter<std::vector<double> >("svpvDistanceCut_lowptD");
    svpvDistanceCut_highptD_ = iConfig.getParameter<std::vector<double> >("svpvDistanceCut_highptD");
    alphaCut_ = iConfig.getParameter<std::vector<double> >("alphaCut");
    MaxDocaCut_ = iConfig.getParameter<std::vector<double> >("MaxDocaCut");
    tktkRes_dCutSeparating_PtVal_ = iConfig.getParameter<std::vector<double> >("tktkRes_dCutSeparating_PtVal");
    tktkRes_dPtCut_ = iConfig.getParameter<std::vector<double> >("tktkRes_dPtCut");
    tktkRes_dEtaCut_ = iConfig.getParameter<std::vector<double> >("tktkRes_dEtaCut");
    tktkRes_VtxChiProbCut_ = iConfig.getParameter<std::vector<double> >("tktkRes_VtxChiProbCut");
    tktkRes_svpvDistanceCut_lowptD_ = iConfig.getParameter<std::vector<double> >("tktkRes_svpvDistanceCut_lowptD");
    tktkRes_svpvDistanceCut_highptD_ = iConfig.getParameter<std::vector<double> >("tktkRes_svpvDistanceCut_highptD");
	tktkRes_alphaCut_ = iConfig.getParameter<std::vector<double> >("tktkRes_alphaCut");
    RunOnMC_ = iConfig.getParameter<bool>("RunOnMC");
    doTkPreCut_ = iConfig.getParameter<bool>("doTkPreCut");
    makeDntuple_ = iConfig.getParameter<bool>("makeDntuple");
    doDntupleSkim_ = iConfig.getParameter<bool>("doDntupleSkim");
    printInfo_ = iConfig.getParameter<bool>("printInfo");
    printFill_Info_ = iConfig.getParameter<bool>("printFill");
	readDedx_ = iConfig.getParameter<bool>("readDedx");
    MVAMapLabelInputTag_ = iConfig.getParameter<edm::InputTag>("MVAMapLabel");
    MVAMapLabel_ = consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("MVAMapLabel"));
    MVAMapLabelpA_ = consumes< std::vector<float> >(iConfig.getParameter<edm::InputTag>("MVAMapLabel"));
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("Dedx_Token1"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("Dedx_Token2"));

    TrackCutLevel       = fs->make<TH1F>("TrackCutLevel"    , "TrackCutLevel"   , 10, 0, 10);
    for(unsigned int i = 0; i < Dchannel_.size(); i++){
        TH1D* DMassCutLevel_temp      = fs->make<TH1D>(TString::Format("DMassCutLevel_i")   ,TString::Format("DMassCutLevel_i")  , 10, 0, 10);
        DMassCutLevel.push_back(DMassCutLevel_temp);
    }
}//}}}


Dfinder::~Dfinder()
{//{{{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}//}}}

//
// member functions
//



// ------------ method called for each event  ------------
void Dfinder::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   
   
   //checking input parameter size
    if((Dchannel_.size() != dCutSeparating_PtVal_.size())
    || (Dchannel_.size() != dPtCut_.size())
    || (Dchannel_.size() != dEtaCut_.size())
    || (Dchannel_.size() != dRapidityCut_.size())
    || (Dchannel_.size() != VtxChiProbCut_.size())
    || (Dchannel_.size() != svpvDistanceCut_lowptD_.size())
    || (Dchannel_.size() != svpvDistanceCut_highptD_.size())
    || (Dchannel_.size() != alphaCut_.size())
    || (Dchannel_.size() != MaxDocaCut_.size())
    || (Dchannel_.size() != tktkRes_dCutSeparating_PtVal_.size())
    || (Dchannel_.size() != tktkRes_dPtCut_.size())
    || (Dchannel_.size() != tktkRes_dEtaCut_.size())
    || (Dchannel_.size() != tktkRes_VtxChiProbCut_.size())
    || (Dchannel_.size() != tktkRes_svpvDistanceCut_lowptD_.size())
    || (Dchannel_.size() != tktkRes_svpvDistanceCut_highptD_.size())
    || (Dchannel_.size() != tktkRes_alphaCut_.size())
    ){
        std::cout<<"Unmatched input parameter vector size, EXIT"<<std::endl;
        return;
	}
	

    //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
    using namespace edm;
    using namespace reco;
    //ESHandle<MagneticField> bField;
    iSetup.get<IdealMagneticFieldRecord>().get(bField);

    // Change used muon and track collections
    edm::Handle< std::vector<pat::GenericParticle> > tks;
    iEvent.getByToken(trackLabel_, tks);

    //CLEAN all memory
    memset(&EvtInfo     ,0x00,sizeof(EvtInfo)   );
    memset(&VtxInfo     ,0x00,sizeof(VtxInfo)   );
    memset(&TrackInfo   ,0x00,sizeof(TrackInfo) );
    memset(&DInfo       ,0x00,sizeof(DInfo)    );
    memset(&GenInfo     ,0x00,sizeof(GenInfo)   );
	
	//auto t0 = std::chrono::high_resolution_clock::now();

    // EvtInfo section{{{
    EvtInfo.RunNo   = iEvent.id().run();
    EvtInfo.EvtNo   = iEvent.id().event();
    //std::cout<<"(EvtInfo.EvtNo)"<<EvtInfo.EvtNo<<std::endl;
    EvtInfo.BxNo    = iEvent.bunchCrossing();
    EvtInfo.LumiNo  = iEvent.luminosityBlock();
    EvtInfo.Orbit   = iEvent.orbitNumber();
    EvtInfo.McFlag  = !iEvent.isRealData();
    //EvtInfo.hltnames->clear();
    //EvtInfo.nTrgBook= N_TRIGGER_BOOKINGS;

    //Using HI HLT analysis now
    /*
    //HLT{{{
    edm::Handle<TriggerResults> TrgResultsHandle; //catch triggerresults
    bool with_TriggerResults = iEvent.getByLabel(hltLabel_,TrgResultsHandle);
    if(!with_TriggerResults){//
        std::cout << "Sorry there is no TriggerResult in the file" << std::endl;
    }else{
        //get the names of the triggers
        const edm::TriggerNames &TrgNames = iEvent.triggerNames(*TrgResultsHandle);
        EvtInfo.trgCount = 0;
        for(int i=0; i< N_TRIGGER_BOOKINGS; i++){
            unsigned int TrgIndex = TrgNames.triggerIndex(TriggerBooking[i]);
            if (TrgIndex == TrgNames.size()) {
                EvtInfo.trgBook[i] = -4; // The trigger path is not known in this event.
            }else if ( !TrgResultsHandle->wasrun( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -3; // The trigger path was not included in this event.
            }else if ( !TrgResultsHandle->accept( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -2; // The trigger path was not accepted in this event.
            }else if (  TrgResultsHandle->error ( TrgIndex ) ) {
                EvtInfo.trgBook[i] = -1; // The trigger path has an error in this event.
            }else {
                EvtInfo.trgBook[i] = +1; // It's triggered.
                EvtInfo.trgCount++; 
            }
        }
        EvtInfo.nHLT = TrgNames.size();
        for(unsigned int i=0; i<TrgNames.size(); i++){
            EvtInfo.hltBits[i] = (TrgResultsHandle->accept(i) == true) ? 1:0;
        }
    }//end(!with_TriggerResults)}}}
    */

    // Handle primary vertex properties
    Vertex thePrimaryV;
    math::XYZPoint RefVtx;
    //get beamspot information
    Vertex theBeamSpotV;
    reco::BeamSpot beamSpot;
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByToken(bsLabel_, beamSpotHandle);
    if (beamSpotHandle.isValid()){
        beamSpot = *beamSpotHandle;
        theBeamSpotV = Vertex(beamSpot.position(), beamSpot.covariance3D());
        EvtInfo.BSx             = beamSpot.x0();
        EvtInfo.BSy             = beamSpot.y0();
        EvtInfo.BSz             = beamSpot.z0();
        EvtInfo.BSxErr          = beamSpot.x0Error();
        EvtInfo.BSyErr          = beamSpot.y0Error();
        EvtInfo.BSzErr          = beamSpot.z0Error();
        EvtInfo.BSdxdz          = beamSpot.dxdz();
        EvtInfo.BSdydz          = beamSpot.dydz();
        EvtInfo.BSdxdzErr       = beamSpot.dxdzError();
        EvtInfo.BSdydzErr       = beamSpot.dydzError();
        EvtInfo.BSWidthX        = beamSpot.BeamWidthX();
        EvtInfo.BSWidthXErr     = beamSpot.BeamWidthXError();
        EvtInfo.BSWidthY        = beamSpot.BeamWidthY();
        EvtInfo.BSWidthYErr     = beamSpot.BeamWidthYError();
    }else{
        std::cout<< "No beam spot available from EventSetup \n";
    }

    //get vertex informationa
    edm::Handle<reco::VertexCollection> VertexHandle;
    iEvent.getByToken(pvLabel_, VertexHandle);

    /*  
    if (!VertexHandle.failedToGet() && VertexHandle->size()>0){
        //int nVtxTrks = 0;//outdated PV definition
        double max_tkSt = 0;
        for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin(); it_vtx != VertexHandle->end(); it_vtx++){
            if (!it_vtx->isValid()) continue;
            //find primary vertex with largest St
            double tkSt = 0;
            for(std::vector<reco::TrackBaseRef>::const_iterator it_tk = it_vtx->tracks_begin();
                it_tk != it_vtx->tracks_end(); it_tk++){
                tkSt += it_tk->get()->pt();
            }
            if (tkSt > max_tkSt){
                max_tkSt = tkSt;
                thePrimaryV = Vertex(*it_vtx);
            }
        }
    }else{ 
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    RefVtx = thePrimaryV.position();
    */

    double PVBS_Pt_Max = -100.;
    reco::Vertex PVtx_BS;
    if( VertexHandle.isValid() && !VertexHandle.failedToGet() && VertexHandle->size() > 0) {
        //const vector<reco::Vertex> VerticesBS = *VertexHandle;
        for(std::vector<reco::Vertex>::const_iterator it_vtx = VertexHandle->begin();it_vtx != VertexHandle->end(); it_vtx++ ) {
        if (VtxInfo.Size>=MAX_Vertices) {
            std::cout << "PVBS " << VtxInfo.Size << std::endl;
            fprintf(stderr,"ERROR: number of  Vertices exceeds the size of array.\n");
            break;//exit(0);
        }
        VtxInfo.isValid[VtxInfo.Size] = it_vtx->isValid();
        VtxInfo.isFake[VtxInfo.Size] = it_vtx->isFake();
        VtxInfo.Ndof[VtxInfo.Size] = it_vtx->ndof();
        VtxInfo.NormalizedChi2[VtxInfo.Size] = it_vtx->normalizedChi2();
        VtxInfo.x[VtxInfo.Size] = it_vtx->x(); 
        VtxInfo.y[VtxInfo.Size] = it_vtx->y();
        VtxInfo.z[VtxInfo.Size] = it_vtx->z();
        VtxInfo.Pt_Sum[VtxInfo.Size] = 0.;
        VtxInfo.Pt_Sum2[VtxInfo.Size] = 0.;
        //if its hiSelectedVertex, then there will be only one vertex and will have no associated tracks
        if(int(VertexHandle->end()-VertexHandle->begin())==1){
            thePrimaryV = *it_vtx;
            VtxInfo.Size++;
            break;
        }
        for (reco::Vertex::trackRef_iterator it = it_vtx->tracks_begin(); it != it_vtx->tracks_end(); it++) {
           VtxInfo.Pt_Sum[VtxInfo.Size] += (*it)->pt();
           VtxInfo.Pt_Sum2[VtxInfo.Size] += ((*it)->pt() * (*it)->pt());
        }
        if( VtxInfo.Pt_Sum[VtxInfo.Size] >= PVBS_Pt_Max ){
            PVBS_Pt_Max = VtxInfo.Pt_Sum[VtxInfo.Size];
            thePrimaryV = *it_vtx;
        }            
        VtxInfo.Size++;
        }
    }else{ 
        thePrimaryV = Vertex(beamSpot.position(), beamSpot.covariance3D());
    }
    RefVtx = thePrimaryV.position();

    EvtInfo.PVx     = thePrimaryV.position().x();
    EvtInfo.PVy     = thePrimaryV.position().y();
    EvtInfo.PVz     = thePrimaryV.position().z();
    EvtInfo.PVxE    = thePrimaryV.xError();
    EvtInfo.PVyE    = thePrimaryV.yError();
    EvtInfo.PVzE    = thePrimaryV.zError();
    EvtInfo.PVnchi2 = thePrimaryV.normalizedChi2();
    EvtInfo.PVchi2  = thePrimaryV.chi2();

    // get pile-up information
    if (!iEvent.isRealData() && RunOnMC_){
        edm::Handle<std::vector<PileupSummaryInfo> >  PUHandle;
        iEvent.getByToken(puInfoLabel_, PUHandle);
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = PUHandle->begin(); PVI != PUHandle->end(); ++PVI) {
            EvtInfo.nPU[EvtInfo.nBX]   = PVI->getPU_NumInteractions();
            EvtInfo.BXPU[EvtInfo.nBX]  = PVI->getBunchCrossing();
            EvtInfo.trueIT[EvtInfo.nBX]= PVI->getTrueNumInteractions();
            EvtInfo.nBX += 1;
        }
    }else{
        EvtInfo.nBX = 0;
    }

    //}}}
    //printf("-----*****DEBUG:End of EvtInfo.\n");

	// Double check size=0.
    TrackInfo.size  = 0;
    DInfo.size      = 0;
    GenInfo.size    = 0;
    
    std::vector<int> D_counter;
    for(unsigned int i = 0; i < Dchannel_.size(); i++){
        D_counter.push_back(0);
    }

    std::vector<pat::GenericParticle>   input_tracks;
    input_tracks = *tks;
    try{
		const reco::GenParticle* genMuonPtr[MAX_MUON];
		memset(genMuonPtr,0x00,MAX_MUON);
		const reco::GenParticle* genTrackPtr[MAX_TRACK];
		memset(genTrackPtr,0x00,MAX_GEN);
		//standard check for validity of input data
		if (0){
			if (printInfo_) std::cout << "There's no muon : " << iEvent.id() << std::endl;
		}else{
			if (input_tracks.size() == 0){
				if (printInfo_) std::cout << "There's no track: " << iEvent.id() << std::endl;
			}else{
				if (printInfo_) std::cout << "Got " << input_tracks.size() << " tracks" << std::endl;
				if (input_tracks.size() > 0){

					//Preselect tracks{{{
					std::vector<bool> isNeededTrack;// Are the tracks redundant?
					vector<Track> lst;
					int PassedTrk = 0;
					for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
							tk_it != input_tracks.end(); tk_it++){
						if(PassedTrk >= MAX_TRACK){
							fprintf(stderr,"ERROR: number of tracks exceeds the size of array.\n");
							break;
						}
						isNeededTrack.push_back(false);
						if (printFill_Info_){
							TrackCutLevel->Fill(0);//number of all tracks
						}
						if (printFill_Info_){
							TrackCutLevel->Fill(1);//
						}
						if (tk_it->pt()<tkPtCut_)                           continue;
						if (printFill_Info_){
							TrackCutLevel->Fill(2);
						}
						if (fabs(tk_it->eta())>tkEtaCut_)                   continue;
						if (printFill_Info_){
							TrackCutLevel->Fill(3);
						}
						//if (fabs(tk_it->eta()) > 2.5)                       continue;
						if (printFill_Info_){
							TrackCutLevel->Fill(4);
						}
						if(doTkPreCut_){
							if( !(tk_it->track()->quality(reco::TrackBase::highPurity))) continue;
							//d0 analysis cuts
							//if(tk_it->track()->hitPattern().numberOfValidHits() < 12) continue;
							//if(tk_it->track()->ptError()/tk_it->track()->pt() > 0.075) continue;
							//if(tk_it->track()->normalizedChi2()/tk_it->track()->hitPattern().trackerLayersWithMeasurement() > 0.25) continue;
							//outdated selections
							//if (tk_it->track()->normalizedChi2()>5)             continue;
							//if (tk_it->p()>200 || tk_it->pt()>200)              continue;
							//if (tk_it->track()->hitPattern().numberOfValidStripHits()<10)continue;
							//if (tk_it->track()->hitPattern().numberOfValidPixelHits()<2) continue;
							if (printFill_Info_){
								TrackCutLevel->Fill(5);
							}
						}
						isNeededTrack[tk_it-input_tracks.begin()] = true;
						lst.push_back(
								Track(input_tracks[tk_it-input_tracks.begin()].pt(),
									input_tracks[tk_it-input_tracks.begin()].eta(),
									input_tracks[tk_it-input_tracks.begin()].phi(),
									input_tracks[tk_it-input_tracks.begin()].charge()+1,
									tk_it-input_tracks.begin()
									)
								);//here store charge+1




                        PassedTrk++;
                    }//end of track preselection}}}
                    //printf("-----*****DEBUG:End of track preselection.\n");
				if(printInfo_) std::cout<<"PassedTrk: "<<PassedTrk<<std::endl;
				
				vector<TrackXYZP2> lstXYZP2;
				int n = (int) lst.size();
				for (int i = 0; i < n; i ++) {
					TrackXYZP2 tr = lst[i];
					lstXYZP2.push_back(tr);
				}

                    //auto t0 = std::chrono::high_resolution_clock::now();
					// DInfo section{{{
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-
                    //////////////////////////////////////////////////////////////////////////
                    float d0_mass_window[2] = {D0_MASS-0.2,D0_MASS+0.2};
                    
                    if(Dchannel_[0] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, d0_mass_window, -1, -1, false, false, 1, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 1, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[1] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, d0_mass_window, -1, -1, false, false, 2, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 2, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
					float dplus_mass_window[2] = {DPLUS_MASS-0.2,DPLUS_MASS+0.2};
                    if(Dchannel_[2] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, dplus_mass_window, -1, -1, false, false, 3, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dplus_mass_window, InVec, -1, -1, false, false, 3, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[3] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, dplus_mass_window, -1, -1, false, false, 4, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dplus_mass_window, InVec, -1, -1, false, false, 4, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K-pi-pi+pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[4] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        std::pair<float, int> tk4 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, d0_mass_window, -1, -1, false, false, 5, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 5, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+pi+pi-pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[5] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk4 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, d0_mass_window, -1, -1, false, false, 6, 0);
                        }
                        //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, d0_mass_window, InVec, -1, -1, false, false, 6, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+K-(Phi)pi+
                    //////////////////////////////////////////////////////////////////////////
					float dsubs_mass_window[2] = {DSUBS_MASS-0.2,DSUBS_MASS+0.2};
                    if(Dchannel_[6] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dsubs_mass_window, PermuVec[i], PHI_MASS, 0.1, false, false, 7, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, dsubs_mass_window, PHI_MASS, 0.1, false, false, 7, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: K+K-(Phi)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[7] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dsubs_mass_window, PermuVec[i], PHI_MASS, 0.1, false, false, 8, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, dsubs_mass_window, PHI_MASS, 0.1, false, false, 8, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0(K-pi+)pi+
                    //////////////////////////////////////////////////////////////////////////
                    float dstar_mass_window[2] = {DSTAR_MASS-0.2,DSTAR_MASS+0.2};
                    if(Dchannel_[8] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 9, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, dstar_mass_window, D0_MASS, 0.1, false, true, 9, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0bar(K+pi-)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[9] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 10, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, theBeamSpotV,lst,lstXYZP2, D_counter, dstar_mass_window, D0_MASS, 0.1, false, true, 10, 1);
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0(K-pi-pi+pi+)pi+
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[10] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk4 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk5 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        InVec.push_back(tk5);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 11, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV,theBeamSpotV,lst,lstXYZP2, D_counter, dstar_mass_window, D0_MASS, 0.1, false, true, 11, 1);
                    }

                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0bar(K+pi+pi-pi-)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[11] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk4 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk5 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        InVec.push_back(tk4);
                        InVec.push_back(tk5);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                        //    Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV, isNeededTrackIdx, D_counter, dstar_mass_window, PermuVec[i], D0_MASS, 0.1, false, true, 12, 0);
                        //}
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV,theBeamSpotV,lst,lstXYZP2, D_counter, dstar_mass_window, D0_MASS, 0.1, false, true, 12, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0bar(K+pi-)pi+
                    //////////////////////////////////////////////////////////////////////////
                    //float bplus_mass_window[2] = {BPLUS_MASS-0.2,BPLUS_MASS+0.2};
                    float bplus_mass_window[2] = {4.5, 6.5};
                    if(Dchannel_[12] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(-PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV,theBeamSpotV,lst,lstXYZP2, D_counter, bplus_mass_window, D0_MASS, 0.1, false, true, 13, 1);
                    }
                    //////////////////////////////////////////////////////////////////////////
                    // RECONSTRUCTION: D0(K-pi+)pi-
                    //////////////////////////////////////////////////////////////////////////
                    if(Dchannel_[13] == 1){
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
                        std::pair<float, int> tk1 = std::make_pair(-KAON_MASS, 1);
                        std::pair<float, int> tk2 = std::make_pair(PION_MASS, 1);
                        std::pair<float, int> tk3 = std::make_pair(-PION_MASS, 0);
                        InVec.push_back(tk1);
                        InVec.push_back(tk2);
                        InVec.push_back(tk3);
                        PermuVec = GetPermu(InVec);
                        PermuVec = DelDuplicate(PermuVec);
                        Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV,theBeamSpotV,lst,lstXYZP2, D_counter, bplus_mass_window, D0_MASS, 0.1, false, true, 14, 1);
                    }
                   ///////////////////////////////////////////////////////////////////////////
                   //RECONSTRUCTION: pi+p+k-(for lambda_C)
                   ///////////////////////////////////////////////////////////////////////////
                  float lambdaC_mass_window[2] = {2.1,2.5};
				  int method = 4;
                 // float lambdaC_mass_window[2] = {0,10};
				  if(Dchannel_[14] == 1){
					    /*
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;
						std::pair<float, int> tk1 = std::make_pair(+PION_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(+PROTON_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(-KAON_MASS, 0);
                        */
						//InVec.push_back(tk1);
                        //InVec.push_back(tk2);
                        //InVec.push_back(tk3);
                        //PermuVec = GetPermu(InVec);
                        //PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                            Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV,theBeamSpotV,lst,lstXYZP2, D_counter, lambdaC_mass_window, -1, -1, false, false, 15, method);
                        //}
                    }
                   ///////////////////////////////////////////////////////////////////////////
                   //RECONSTRUCTION: pi-pbar-k+(for anti_lambda_C)
                   ///////////////////////////////////////////////////////////////////////////
                   if(Dchannel_[15] == 1){
					   /*
                        std::vector< std::vector< std::pair<float, int> > > PermuVec;
                        std::vector< std::pair<float, int> > InVec;   
						std::pair<float, int> tk1 = std::make_pair(-PION_MASS, 0);
                        std::pair<float, int> tk2 = std::make_pair(-PROTON_MASS, 0);
                        std::pair<float, int> tk3 = std::make_pair(+KAON_MASS, 0);
                        */
						//InVec.push_back(tk1);
                        //InVec.push_back(tk2);
                        //InVec.push_back(tk3);
                        //PermuVec = GetPermu(InVec);
                        //PermuVec = DelDuplicate(PermuVec);
                        //for(unsigned int i = 0; i < PermuVec.size(); i++){
                           //Dfinder::BranchOutNTk( DInfo, input_tracks, thePrimaryV,lst,lstXYZP2, D_counter, lambdaC_mass_window, -1, -1, false, false, 16, method);
                        //}
                     }
                  







                    if(printInfo_){
                        printf("D_counter: ");
                        for(unsigned int i = 0; i < Dchannel_.size(); i++){
                            printf("%d/", D_counter[i]);
                        }
                        printf("\n");

                    }//}}}
                    //printf("-----*****DEBUG:End of DInfo.\n");

                    
					// TrackInfo section {{{
                    // Setup MVA
                    Handle<edm::ValueMap<float> > mvaoutput;
                    Handle< std::vector<float> > mvaoutputpA;
                    std::vector<float>   mvavector;
                    if(MVAMapLabelInputTag_.instance() == "MVAVals") {
                        iEvent.getByToken(MVAMapLabel_, mvaoutput);
                    }
                    if(MVAMapLabelInputTag_.instance() == "MVAValues") {
                        iEvent.getByToken(MVAMapLabelpA_, mvaoutputpA);
                        mvavector = *mvaoutputpA;
                        assert(mvavector.size()==input_tracks.size());
                    }

                    // Setup Dedx
                    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
                    edm::ValueMap<reco::DeDxData> dEdxTrack1;
                    //edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
                    //edm::ValueMap<reco::DeDxData> dEdxTrack2;
                    if(readDedx_) {
                        iEvent.getByToken(Dedx_Token1_, dEdxHandle1);
                        dEdxTrack1 = *dEdxHandle1.product();
                        //iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
                        //dEdxTrack2 = *dEdxHandle2.product();
                    }
					

                    for(std::vector<pat::GenericParticle>::const_iterator tk_it=input_tracks.begin();
                        tk_it != input_tracks.end() ; tk_it++){
                        int tk_hindex = int(tk_it - input_tracks.begin());
                        if(tk_hindex>=int(isNeededTrack.size())) break;
                        if (isNeededTrack[tk_hindex]==false) continue;

                        //Create list of relative xb candidates for later filling
                        std::vector<int> listOfRelativeDCand1;//1~nXb
                        std::vector<int> listOfRelativeDCand2;//1~nXb
                        std::vector<int> listOfRelativeDCand3;//1~nXb
                        std::vector<int> listOfRelativeDCand4;//1~nXb
                        std::vector<int> listOfRelativeDCand5;//1~nXb
                        std::vector<int> listOfRelativeDResCand1;//1~nXb
                        std::vector<int> listOfRelativeDResCand2;//1~nXb
                        std::vector<int> listOfRelativeDResCand3;//1~nXb
                        std::vector<int> listOfRelativeDResCand4;//1~nXb
                        for(int d=0; d < DInfo.size; d++){
                            if(DInfo.rftk1_index[d] == tk_hindex){
                                listOfRelativeDCand1.push_back(d+1);
                            }
                            if(DInfo.rftk2_index[d] == tk_hindex){
                                listOfRelativeDCand2.push_back(d+1);
                            }
                            if(DInfo.rftk3_index[d] == tk_hindex){
                                listOfRelativeDCand3.push_back(d+1);
                            }
                            if(DInfo.rftk4_index[d] == tk_hindex){
                                listOfRelativeDCand4.push_back(d+1);
                            }
                            if(DInfo.rftk5_index[d] == tk_hindex){
                                listOfRelativeDCand5.push_back(d+1);
                            }

                            if(DInfo.tktkRes_rftk1_index[d] == tk_hindex){
                                listOfRelativeDResCand1.push_back(d+1);
                            }
                            if(DInfo.tktkRes_rftk2_index[d] == tk_hindex){
                                listOfRelativeDResCand2.push_back(d+1);
                            }
                            if(DInfo.tktkRes_rftk3_index[d] == tk_hindex){
                                listOfRelativeDResCand3.push_back(d+1);
                            }
                            if(DInfo.tktkRes_rftk4_index[d] == tk_hindex){
                                listOfRelativeDResCand4.push_back(d+1);
                            }
                        }
                        
                        if(dropUnusedTracks_ && listOfRelativeDCand1.size() == 0 && listOfRelativeDCand2.size() == 0 && listOfRelativeDCand3.size() == 0 && listOfRelativeDCand4.size() == 0 && listOfRelativeDCand5.size() == 0 && listOfRelativeDResCand1.size() == 0 && listOfRelativeDResCand2.size() == 0 && listOfRelativeDResCand3.size() == 0 && listOfRelativeDResCand4.size() == 0) continue;//drop unused tracks

                        TrackInfo.index          [TrackInfo.size] = TrackInfo.size;
                        TrackInfo.handle_index   [TrackInfo.size] = tk_hindex;
                        TrackInfo.charge         [TrackInfo.size] = tk_it->charge();
                        TrackInfo.pt             [TrackInfo.size] = tk_it->pt();
                        TrackInfo.eta            [TrackInfo.size] = tk_it->eta();
                        TrackInfo.phi            [TrackInfo.size] = tk_it->phi();
                        TrackInfo.ptErr          [TrackInfo.size] = tk_it->track()->ptError();
                        TrackInfo.etaErr         [TrackInfo.size] = tk_it->track()->etaError();
                        TrackInfo.phiErr         [TrackInfo.size] = tk_it->track()->phiError();
                        //TrackInfo.p              [TrackInfo.size] = tk_it->p();
                        TrackInfo.striphit       [TrackInfo.size] = tk_it->track()->hitPattern().numberOfValidStripHits();
                        TrackInfo.pixelhit       [TrackInfo.size] = tk_it->track()->hitPattern().numberOfValidPixelHits();
                        TrackInfo.nStripLayer    [TrackInfo.size] = tk_it->track()->hitPattern().stripLayersWithMeasurement();
                        TrackInfo.nPixelLayer    [TrackInfo.size] = tk_it->track()->hitPattern().pixelLayersWithMeasurement();
						//TrackInfo.fpbarrelhit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInFirstPixelBarrel();
						//TrackInfo.fpendcaphit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInFirstPixelEndcap();
						TrackInfo.fpbarrelhit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel,1);
						TrackInfo.fpendcaphit    [TrackInfo.size] = tk_it->track()->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelEndcap,1);
						
						TrackInfo.chi2           [TrackInfo.size] = tk_it->track()->chi2();
                        TrackInfo.ndf            [TrackInfo.size] = tk_it->track()->ndof();
                        TrackInfo.d0             [TrackInfo.size] = tk_it->track()->d0();
                        TrackInfo.d0error        [TrackInfo.size] = tk_it->track()->d0Error();
                        TrackInfo.dzPV           [TrackInfo.size] = tk_it->track()->dz(RefVtx);
                        TrackInfo.dxyPV          [TrackInfo.size] = tk_it->track()->dxy(RefVtx);
                        TrackInfo.highPurity     [TrackInfo.size] = tk_it->track()->quality(reco::TrackBase::highPurity);
                        TrackInfo.geninfo_index  [TrackInfo.size] = -1;//initialize for later use
                        if(MVAMapLabelInputTag_.instance() == "MVAVals") 
                        TrackInfo.trkMVAVal      [TrackInfo.size] = (*mvaoutput)[tk_it->track()];
                        if(MVAMapLabelInputTag_.instance() == "MVAValues") 
                        TrackInfo.trkMVAVal      [TrackInfo.size] = mvavector[tk_hindex];
                        TrackInfo.trkAlgo        [TrackInfo.size] = tk_it->track()->algo();
                        TrackInfo.originalTrkAlgo[TrackInfo.size] = tk_it->track()->originalAlgo();
                        if(readDedx_) {
                            TrackInfo.dedx           [TrackInfo.size] = dEdxTrack1[tk_it->track()].dEdx();
                        }else 
                            TrackInfo.dedx           [TrackInfo.size] = -1;

                        if(tk_it->track().isNonnull()){
                            for(int tq = 0; tq < reco::TrackBase::qualitySize; tq++){
                            if (tk_it->track()->quality(static_cast<reco::TrackBase::TrackQuality>(tq))) TrackInfo.trackQuality[TrackInfo.size] += 1 << (tq);
                        }}

                        if (!iEvent.isRealData() && RunOnMC_)
                            genTrackPtr [TrackInfo.size] = tk_it->genParticle();
//cout<<"genTrack"<<genTrackPtr [TrackInfo.size]<<endl;
                        // Fill the same list for DInfo
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand1.size(); iCands++){
                            DInfo.rftk1_index[listOfRelativeDCand1[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand2.size(); iCands++){
                            DInfo.rftk2_index[listOfRelativeDCand2[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand3.size(); iCands++){
                            DInfo.rftk3_index[listOfRelativeDCand3[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand4.size(); iCands++){
                            DInfo.rftk4_index[listOfRelativeDCand4[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDCand5.size(); iCands++){
                            DInfo.rftk5_index[listOfRelativeDCand5[iCands]-1] = TrackInfo.size;
                        }

                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand1.size(); iCands++){
                            DInfo.tktkRes_rftk1_index[listOfRelativeDResCand1[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand2.size(); iCands++){
                            DInfo.tktkRes_rftk2_index[listOfRelativeDResCand2[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand3.size(); iCands++){
                            DInfo.tktkRes_rftk3_index[listOfRelativeDResCand3[iCands]-1] = TrackInfo.size;
                        }
                        for(unsigned int iCands=0; iCands < listOfRelativeDResCand4.size(); iCands++){
                            DInfo.tktkRes_rftk4_index[listOfRelativeDResCand4[iCands]-1] = TrackInfo.size;
                        }
                        TrackInfo.size++;
                    }//end of TrackInfo}}}
                    //printf("-----*****DEBUG:End of TrackInfo.\n");
                }//has nTracks>1
            }//if no Tracks

        }//if no Muons

	


        // GenInfo section{{{
        if (!iEvent.isRealData() && RunOnMC_){
        //if (RunOnMC_){
        //if (1){
            //edm::Handle< std::vector<reco::GenParticle> > gens;
            edm::Handle<reco::GenParticleCollection> gens;
            iEvent.getByToken(genLabel_, gens);

            std::vector<const reco::Candidate *> sel_cands;
            //deprecated
            /*
            std::vector<const reco::Candidate *> cands;
            for(std::vector<reco::GenParticle>::const_iterator it_gen = gens->begin();
                it_gen != gens->end(); it_gen++ ){
                cands.push_back(&*it_gen);
            }
            */

			//fill gen PV, parton, gluon, etc can also be used
			GenInfo.genPVx = -99;
			GenInfo.genPVy = -99;
			GenInfo.genPVz = -99;
            for(std::vector<reco::GenParticle>::const_iterator it_gen=gens->begin();
                it_gen != gens->end(); it_gen++){
				
				//pt 0 particle should be from beam. Apart from proton and neutron, they are also pointing to PV actually
				if( !(it_gen->pt() > 0) ) continue;
				
				//take the vertex of first produced particle, should be enough. Better to check 2-3 particles if want
				GenInfo.genPVx = it_gen->vx();
				GenInfo.genPVy = it_gen->vy();
				GenInfo.genPVz = it_gen->vz();
				//if( abs(it_gen->pdgId()) == 421 ) cout << "DzeroDzeroDzeroDzero!!!!!!" << endl;
				//cout << "Pid: " << it_gen->pdgId() << " status: " << it_gen->status() << endl;
				//cout << " vx: " << it_gen->vx() << " vy: " << it_gen->vy() << " vz: " << it_gen->vz() << endl;
				break;
			}
			//end fill gen PV

            for(std::vector<reco::GenParticle>::const_iterator it_gen=gens->begin();
                it_gen != gens->end(); it_gen++){
		//	cout<<"particles id "<<it_gen->pdgId()<<" , status = "<<it_gen->status()<<endl ;//I tested that the daughter particles are really in the list.
                //if (it_gen->status() > 2 && it_gen->status() != 8) continue;//only status 1, 2, 8(simulated)
                if(GenInfo.size >= MAX_GEN){
                    fprintf(stderr,"ERROR: number of gens exceeds the size of array.\n");
                    break;;
                }
    /////
	



                /*
                if (
                    //(abs(it_gen->pdgId()) == 111 && it_gen->status() == 2) ||//pi 0
                    (abs(it_gen->pdgId()) == 211 && it_gen->status() == 2) ||//pi +-
                    //(abs(it_gen->pdgId()) == 311 && it_gen->status() == 2) ||//K0
                    (abs(it_gen->pdgId()) == 321 && it_gen->status() == 2) //K+-
                ) continue;//only status=1 pi+- and K+-
                */

                bool isGenSignal = false;
                //save target intermediat state particle
                if (
                    abs(int(it_gen->pdgId()/100) % 100) == 4  ||//c menson
                    abs(int(it_gen->pdgId()/100) % 100) == 5  ||//b menson
                    abs(it_gen->pdgId()) == 511 ||//B_0
                    abs(it_gen->pdgId()) == 521 ||//B_+-
                    abs(it_gen->pdgId()) == 531 ||//B_s
                    abs(it_gen->pdgId()) == 130 ||//KL
                    abs(it_gen->pdgId()) == 4122 ||//lamadac
                    //abs(it_gen->pdgId()) == 311 ||//K0
                    //abs(it_gen->pdgId()) == 321 ||//K+
                    //abs(it_gen->pdgId()) == 310 ||//KS
                    //abs(it_gen->pdgId()) == 313 ||//K*0(892)
                    //abs(it_gen->pdgId()) == 323 ||//K*+-(892)
                    //abs(it_gen->pdgId()) == 333 ||//phi(1020)
                    it_gen->pdgId() == 443      ||//Jpsi
                    it_gen->pdgId() == 100443   ||//Psi(2S)
                    it_gen->pdgId() == 553      ||//Upsilon
                    it_gen->pdgId() == 100553     //Upsilon(2S)
                   ) isGenSignal = true;//b, c, s mesons

                if (abs(it_gen->pdgId()) == 13) isGenSignal = true;//all mu

                if (
                    abs(int(it_gen->pdgId()/100) % 100) == 3  ||//s menson
                    abs(it_gen->pdgId()) == 111 || //pi0
                    abs(it_gen->pdgId()) == 211 ||//pi+
                    abs(it_gen->pdgId()) == 2212  ||//proton
					abs(it_gen->pdgId()) == 321 || //K+ actually, I DO NOT need to add this,s meson is in the list
					abs(it_gen->pdgId()) ==3124 || //lambda(1520)0
					abs(it_gen->pdgId()) ==2224 //delta++
                    ){
                    reco::GenParticle _deRef = (*it_gen);
                    reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
                    //std::cout<<Myself->pdgId()<<"-----------"<<std::endl;
                  //  isGenSignal = (Functs.GetAncestor(Myself, 5) | Functs.GetAncestor(Myself, 4));
				    isGenSignal = (Functs.GetAncestor(Myself, 41));
				//	cout<<" isGenSignal "<<isGenSignal<<" ID "<<abs(it_gen->pdgId())<<endl;
                }//all pi and K from b or c meson
				
				if( !isGenSignal){ //should be OK to require the PID > 400
					reco::GenParticle _deRef = (*it_gen);
					reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
					isGenSignal = Functs.GetDescendant(Myself, 41);
				}//other particles (with pid > 400) have D meson in descendant chain
                
///////here I want to test///
//if(abs(it_gen->pdgId())==4122){
//	cout<<" isGenSignal "<< isGenSignal<<endl;
//}
///Till now there is no bugs.

                if (!isGenSignal) continue;

                /*deprecated
                int iMo1 = -1,  iMo2 = -1,  iDa1 = -1,  iDa2 = -1;
                for(std::vector<const reco::Candidate *>::iterator iCands = cands.begin();
                    iCands != cands.end(); iCands++){
                    if (it_gen->numberOfMothers() >= 2){
                        if (it_gen->mother(0) == *iCands)
                            iMo1 = iCands - cands.begin();
                        if (it_gen->mother(1) == *iCands)
                            iMo2 = iCands - cands.begin();
                    }else if(it_gen->numberOfMothers() == 1){
                        if (it_gen->mother(0) == *iCands)
                            iMo1 = iCands - cands.begin();
                    }
                    if (it_gen->numberOfDaughters() >= 2){
                        if (it_gen->daughter(0) == *iCands)
                            iDa1 = iCands - cands.begin();
                        else if (it_gen->daughter(1) == *iCands)
                            iDa2 = iCands - cands.begin();
                    }else if(it_gen->numberOfDaughters() == 1){
                        if (it_gen->daughter(0) == *iCands)
                            iDa1 = iCands - cands.begin();
                    }
                }
                */
                //Find all other particle in TrackInfo
                //printf("-----*****DEBUG:Start of matching.\n");
                //printf("-----*****DEBUG:Entered pion matching block.\n");
                for(int trackIdx = 0; trackIdx < TrackInfo.size; trackIdx++){
                    //match by pat::GenericParticle
                    if (genTrackPtr[trackIdx] == 0 ) continue;
                    if (it_gen->p4() == genTrackPtr[trackIdx]->p4()){
                        TrackInfo.geninfo_index[trackIdx] = GenInfo.size;
                        break;
                    }
                }

                GenInfo.index[GenInfo.size]         = GenInfo.size;
                GenInfo.handle_index[GenInfo.size]  = it_gen-gens->begin();
                GenInfo.pt[GenInfo.size]            = it_gen->pt();
                GenInfo.px[GenInfo.size]            = it_gen->px();
				GenInfo.py[GenInfo.size]            = it_gen->py();
				GenInfo.pz[GenInfo.size]            = it_gen->pz();
				GenInfo.eta[GenInfo.size]           = it_gen->eta();
                GenInfo.phi[GenInfo.size]           = it_gen->phi();
                GenInfo.mass[GenInfo.size]          = it_gen->mass();
                GenInfo.pdgId[GenInfo.size]         = it_gen->pdgId();
                GenInfo.status[GenInfo.size]        = it_gen->status();
                GenInfo.collisionId[GenInfo.size]   = it_gen->collisionId();
                GenInfo.nMo[GenInfo.size]           = it_gen->numberOfMothers();
                GenInfo.nDa[GenInfo.size]           = it_gen->numberOfDaughters();
				GenInfo.vtxX[GenInfo.size]          = it_gen->vx(); //it should be the production vx of the particle, better to double check
				GenInfo.vtxY[GenInfo.size]          = it_gen->vy();
				GenInfo.vtxZ[GenInfo.size]          = it_gen->vz();
                //GenInfo.mo1[GenInfo.size]           = iMo1;//To be matched later.
                //GenInfo.mo2[GenInfo.size]           = iMo2;
                //GenInfo.da1[GenInfo.size]           = iDa1;
                //GenInfo.da2[GenInfo.size]           = iDa2;
                GenInfo.mo1[GenInfo.size]           = -1;//To be matched later.
                GenInfo.mo2[GenInfo.size]           = -1;
                GenInfo.da1[GenInfo.size]           = -1;
                GenInfo.da2[GenInfo.size]           = -1;
                GenInfo.da3[GenInfo.size]           = -1;
                GenInfo.da4[GenInfo.size]           = -1;
                GenInfo.size++;
                sel_cands.push_back(&*it_gen);
            }
            //printf("-----*****DEBUG:End of gens loop.\n");

            int geninfo_idx = 0;
            for(std::vector<const reco::Candidate *>::iterator sCands = sel_cands.begin();
                sCands != sel_cands.end(); sCands++){
                geninfo_idx = int(sCands-sel_cands.begin());
				//cout<<" geninfo_idx "<<geninfo_idx<< " particle "<<(*sCands)->pdgId()<<endl;
				//till now there should be no bugs.
			 for(int nGenMo = 0; nGenMo < std::min(2,int((*sCands)->numberOfMothers())); nGenMo++){

//if((*sCands)->numberOfMothers()==1){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
                        if((*sCands)->mother(nGenMo) == *mCands){
                        //if((*sCands)->mother(0) == *mCands){
                            if(nGenMo == 0) GenInfo.mo1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenMo == 1) GenInfo.mo2[geninfo_idx] = int(mCands-sel_cands.begin());
                        }
                    }
                }
			//	cout<<" geninfo_idx "<<geninfo_idx<< " particle "<<(*sCands)->pdgId()<<" numbers "<<(*sCands)->numberOfDaughters()<<endl;
                for(int nGenDa = 0; nGenDa < std::min(4,int((*sCands)->numberOfDaughters())); nGenDa++){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
						
					//	cout<<" a "<<(*sCands)->daughter(nGenDa)<<" particle "<<(*sCands)->pdgId()<<endl;
						if(abs((*sCands)->pdgId())==4122){
                        //   cout<<" a "<<(*sCands)->daughter(nGenDa)<< " b "<<*mCands<<endl;
						}

                        if((*sCands)->daughter(nGenDa) == *mCands){
						//if((*sCands)->daughter(0)==*mCands){
                            if(nGenDa == 0) GenInfo.da1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 1) GenInfo.da2[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 2) GenInfo.da3[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 3) GenInfo.da4[geninfo_idx] = int(mCands-sel_cands.begin());
		              //  cout<< " geninfoinsidetheloop " << geninfo_idx <<endl;From the test here, it shows that lambdaC not in this loop.
		            //    	cout<<"nGenDa"<<nGenDa<<"   da1 "<<GenInfo.da1[geninfo_idx]<<" geninfo_idx "<< geninfo_idx<<endl;
						//	cout<<"nGenDa"<<nGenDa<<"   da2 "<<GenInfo.da2[geninfo_idx]<<" geninfo_idx "<< geninfo_idx<<endl;
                        //    cout<<"nGenDa"<<nGenDa<<"   da3 "<<GenInfo.da3[geninfo_idx]<<" geninfo_idx "<< geninfo_idx<<endl;
						//	cout<<"nGenDa"<<nGenDa<<"   da4 "<<GenInfo.da4[geninfo_idx]<<" geninfo_idx "<< geninfo_idx<<endl;
							//std::cout << "mCands"<<mCands<<std::endl;
                        }
						//cout<<"daupdgId"<<GenInfo.da1[geninfo_idx]<<endl;
                    }
                }
            }
			///here I add a few lines:
			for(int igen = 0; igen < GenInfo.size; igen++){
                 if (abs(GenInfo.pdgId[igen])==4122){

			//	 cout<<"daupdgId "<<GenInfo.da1[igen]<<" dau2 "<<GenInfo.da2[igen]  <<endl;
		          }
		    }
            /*deprecated
            //Pass handle_index to igen
            for(int igen = 0; igen < GenInfo.size; igen++){
                int iMo1 = GenInfo.mo1[igen];
                int iMo2 = GenInfo.mo2[igen];
                int iDa1 = GenInfo.da1[igen];
                int iDa2 = GenInfo.da2[igen];
                for(int k = 0; k < GenInfo.size; k++){
                    if (iMo1 == GenInfo.handle_index[k])
                        GenInfo.mo1[igen] = k;
                    else if (iMo2 == GenInfo.handle_index[k])
                        GenInfo.mo2[igen] = k;
                    else if (iDa1 == GenInfo.handle_index[k])
                        GenInfo.da1[igen] = k;
                    else if (iDa2 == GenInfo.handle_index[k])
                        GenInfo.da2[igen] = k;
                }
                //In case that GEN particles are omitted from GenInfo
                //handle_index couldn't be the same as igen
                //since the very first proton pair has status 3.
                if (iMo1 == GenInfo.mo1[igen])
                    GenInfo.mo1[igen] = -1;
                if (iMo2 == GenInfo.mo2[igen])
                    GenInfo.mo2[igen] = -1;
                if (iDa1 == GenInfo.da1[igen])
                    GenInfo.da1[igen] = -1;
                if (iDa2 == GenInfo.da2[igen])
                   GenInfo.da2[igen] = -1;
            }
            //printf("-----*****DEBUG:End of IndexToIgen\n");
            */

        }//isRealData}}}
        //printf("-----*****DEBUG:End of GenInfo.\n");
        //std::cout<<"Start to fill!\n";

    }//try
    catch (std::exception & err){
            std::cout  << "Exception during event number: " << iEvent.id()
                << "\n" << err.what() << "\n";
    }//catch 
    root->Fill();
/*
	auto t1 = std::chrono::high_resolution_clock::now(); 
	double dt = 1e-3 * std::chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
	std::cout<<"duration: "<<dt<<std::endl;
*/
/*	
	auto t1 = std::chrono::high_resolution_clock::now(); 
	double dt = 1e-3 * std::chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
	cout<<"duration: "<<dt<<endl;
*/
	//std::cout<<"filled!\n";
 
    //Made a Dntuple on the fly   
    if(makeDntuple_){
        int isDchannel[16];
        isDchannel[0] = 1; //k+pi-
        isDchannel[1] = 1; //k-pi+
        isDchannel[2] = 1; //k-pi+pi+
        isDchannel[3] = 1; //k+pi-pi-
        isDchannel[4] = 1; //k-pi-pi+pi+
        isDchannel[5] = 1; //k+pi+pi-pi-
        isDchannel[6] = 1; 
        isDchannel[7] = 1; 
        isDchannel[8] = 1; 
        isDchannel[9] = 1; 
        isDchannel[10] = 1; 
        isDchannel[11] = 1;
        isDchannel[12] = 1; //B+(D0(k-pi+)pi+)
        isDchannel[13] = 1; //B-(D0(k-pi+)pi-)
        isDchannel[14] = 1; //lambdaC(p+pi+k-)
        isDchannel[15] = 1; //lambdaC(p-k+pi-)
        bool REAL = ((!iEvent.isRealData() && RunOnMC_) ? false:true);
        Dntuple->makeDNtuple(isDchannel, REAL, doDntupleSkim_, &EvtInfo, &VtxInfo, &TrackInfo, &DInfo, &GenInfo, ntD1, ntD2, ntD3, ntD4, ntD5, ntD6, ntD7, ntD8);
        if(!REAL) Dntuple->fillDGenTree(ntGen, &GenInfo);
    }

}


// ------------ method called once each job just after ending the event loop  ------------{{{
void Dfinder::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void Dfinder::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void Dfinder::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void Dfinder::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void Dfinder::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Dfinder::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}//}}}

//Functions{{{
std::vector< std::vector< std::pair<float, int> > > Dfinder::GetPermu(std::vector< std::pair<float, int> > InVec){
    if(InVec.size() == 1){
        std::vector< std::vector< std::pair<float, int> > > OneEntryVec;
        OneEntryVec.push_back(InVec);
        return OneEntryVec;
    }
    std::vector< std::vector< std::pair<float, int> > > NPermu;
    for(unsigned int i = 0; i < InVec.size(); i++){
        std::vector< std::pair<float, int> > copy;
        copy = InVec;
        copy.erase(copy.begin()+i);
        std::vector< std::vector< std::pair<float, int> > > Nminus1Permu;
        Nminus1Permu = GetPermu(copy);
        for(unsigned int j = 0; j < Nminus1Permu.size(); j++){
            Nminus1Permu[j].push_back(InVec[i]);
        }
        NPermu.insert(NPermu.end(), Nminus1Permu.begin(), Nminus1Permu.end());
    }
    return NPermu;
}

std::vector< std::vector< std::pair<float, int> > > Dfinder::DelDuplicate(std::vector< std::vector< std::pair<float, int> > > InVec){
    std::vector< std::vector< std::pair<float, int> > > CleanedVec;
    for(unsigned int i = 0; i < InVec.size(); i++){
        bool IsDuplicate = false;
        for(unsigned int j = 0; j < CleanedVec.size(); j++){
            bool ADuplicate = true;
            for(unsigned int k = 0; k < InVec[i].size(); k++){
                if(InVec[i][k] != CleanedVec[j][k]) {
                    ADuplicate = false;
                    break;
                }
            }
            if(ADuplicate) {
                IsDuplicate = ADuplicate;
                break;
            }
        }
        if(!IsDuplicate) CleanedVec.push_back(InVec[i]);
    }
    return CleanedVec;
}


//long long int cnt0[2], cnt1[2], cnt2[2],cnt3[2];
long long int cnt0=0;
void Dfinder::TkCombinationPermutation_Lc_v3(
		std::vector<pat::GenericParticle> input_tracks, 
		vector<Track> lst,
		vector<TrackXYZP2> lstXYZP2,
		float *mass_window,
		double tktkRes_mass,
		double tktkRes_mass_window,
		vector<TrackSelected> &selectedTkhidxSet,
		int Dchannel_number
		){
	
	int tk1_hindex = -1;
	int tk2_hindex = -1;
	int tk3_hindex = -1;
	
    //auto t0 = std::chrono::high_resolution_clock::now();	
	int number_NeededTrack = (int) lst.size();
	for(int tk1idx = 0; tk1idx < number_NeededTrack; tk1idx++){
		const TrackXYZP2& tr1 = lstXYZP2[tk1idx];
		tk1_hindex = tr1.index;
		int perm1 = tr1.q;  // q1+1
		double p1sq = tr1.p2;
		for(int tk2idx = tk1idx+1; tk2idx < number_NeededTrack; tk2idx++){
			const TrackXYZP2& tr2 = lstXYZP2[tk2idx];
			tk2_hindex = tr2.index;
			int perm2 = (perm1 << 1) + tr2.q; // 2(q1+1) + (q2+1)
			double p2sq = tr2.p2;
			P3 p12(tr2);
			p12 += tr1;
			for(int tk3idx = tk2idx+1; tk3idx < number_NeededTrack; tk3idx++){
				const TrackXYZP2& tr3 = lstXYZP2[tk3idx];
				tk3_hindex = tr3.index;
				int perm3 = (perm2 << 1) + tr3.q; // 4(q1+1) + 2(q2+1) + (q3+1)
				if (perm3 == 0 || perm3 == 14) continue; //this remove the useless permutation
				double p3sq = tr3.p2;
				P3 pD(tr3);
				pD += p12;

			  // pT cut (independent of permutations)
			  double ptD2 = pD.px * pD.px + pD.py * pD.py;
			  if(ptD2 < (dPtCut_[Dchannel_number-1])*(dPtCut_[Dchannel_number-1]))continue;
			  double pzD2 = pD.pz * pD.pz;
			  for (int p = 0; p < 2; p ++) {
			    double p0 = Functs.totE(perm3 + p, p1sq, p2sq, p3sq);
				double mD2 = p0 * p0 - ptD2 - pzD2;
				if(mD2 <(mass_window[0])*(mass_window[0]) || mD2 >(mass_window[1])*(mass_window[1])) continue;
				double mtD2 = mD2 + ptD2;
				//if (pzD2 > shymax2 * mtD2) continue; //this rapdity needs to be changed.
				if (pzD2 > sinh(dRapidityCut_[Dchannel_number-1])*sinh(dRapidityCut_[Dchannel_number-1]) * mtD2) continue;
				cnt0++;
				
				
				
				selectedTkhidxSet.push_back(TrackSelected(tk1_hindex,tk2_hindex,tk3_hindex,perm3+p));//here also need to store the permutation number//here already have all the permutation
				
				
				continue;
			  }//p
			}//tk3id
		}//tk2id
	}//tk1id
	//std::cout<<"TkCombinationPermutation, selectedTkhidxSet.size: "<<selectedTkhidxSet.size()<<std::endl;
/*
	auto t1 = std::chrono::high_resolution_clock::now(); 
	double dt = 1e-3 * std::chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
	std::cout<<"duration: "<<dt<<std::endl;
*/
    if (printFill_Info_)
	{
		//if want to do the counter, this part needs to be modified.
	DMassCutLevel[15]->SetBinContent(1,cnt0);
	//DMassCutLevel[Dchannel_number-1]->SetBinContent(2,cnt1[Dchannel_number-15]);
	//DMassCutLevel[Dchannel_number-1]->SetBinContent(3,cnt2[Dchannel_number-15]);
	//DMassCutLevel[Dchannel_number-1]->SetBinContent(4,cnt3[Dchannel_number-15]);
	}

	return;

}

//BranchOutNTk{{{
void Dfinder::BranchOutNTk(//input 2~4 tracks
    DInfoBranches &DInfo, 
    std::vector<pat::GenericParticle> input_tracks, 
    reco::Vertex thePrimaryV,
	reco::Vertex theBeamSpotV,
	vector<Track> lst,
	vector<TrackXYZP2> lstXYZP2,
    std::vector<int> &D_counter,
    float *mass_window,
    double tktkRes_mass,
    double tktkRes_mass_window,
    bool doConstrainFit,
    bool SequentialFit,
    int Dchannel_number,
	int TkCombinationMethod
){
    if(Dchannel_number > (int)Dchannel_.size()){ printf("Exceeding defined # of channel, exit"); return;}

	vector<TrackSelected > selectedTkhidxSet;
	
	

	if( TkCombinationMethod == 4 ) 
        TkCombinationPermutation_Lc_v3( input_tracks,lst,lstXYZP2, mass_window, tktkRes_mass, tktkRes_mass_window, selectedTkhidxSet, Dchannel_number);
    
	else
		{ printf("unknown method of track combination, exit"); return;}
	
	float chi = 0.;
	float ndf = 0.;
    //auto t0 = std::chrono::high_resolution_clock::now();
    //particle factory: produce transient tracks
    KinematicParticleFactoryFromTransientTrack pFactory;
    VirtualKinematicParticleFactory vFactory;
    //fitter for D
    KinematicParticleVertexFitter   tktk_fitter;
    RefCountedKinematicTree         tktk_VFT;
    RefCountedKinematicParticle     tktk_VFP;
    RefCountedKinematicVertex       tktk_VFPvtx;
    //constrain fit fitter
    KinematicConstrainedVertexFitter kcv_tktk_fitter;
    //fitter for Res
    KinematicParticleVertexFitter   tktkRes_fitter;
    RefCountedKinematicTree         tktkRes_VFT;
    RefCountedKinematicParticle     tktkRes_VFP;
    RefCountedKinematicVertex       tktkRes_VFPvtx;

	//for DCA and it error 
	const MagneticField *field = bField.product();
	AnalyticalImpactPointExtrapolator extrapolator(field);
	TransverseImpactPointExtrapolator transverseExtrapolator(field);
	TrajectoryStateOnSurface tsos;
	TrajectoryStateOnSurface tsos_2D;
    //////
	
    TLorentzVector v4_tk;
    std::vector<TLorentzVector> tktk_4vecs;//fitted tks
    TLorentzVector tktk_4vec;//fitted D
    TLorentzVector unfitted_tktk_4vec;//unfitted D
    std::vector<TLorentzVector> tktkRes_4vecs;//fitted Res tks
    TLorentzVector tktkRes_4vec;//fitted Res
    TLorentzVector unfitted_tktkRes_4vec;//unfitted Res
    std::vector<RefCountedKinematicParticle> tktk_candidate;//input tracks to D fitter
    std::vector<RefCountedKinematicParticle> tktkRes_candidate;//input tracks to Res fitter
    std::vector<RefCountedKinematicParticle> tktkCands;//output tracks from D fitter
    std::vector<RefCountedKinematicParticle> tktkResCands;//output tracks from Res fitter
    TLorentzVector temp_vec;//for temporary usage

    for(int i = 0; i < int(selectedTkhidxSet.size()); i++){
        if (DInfo.size >= MAX_XB) break;


        //clear before using
        v4_tk.Clear();
        tktk_4vecs.clear();
        tktk_4vec.Clear();
        unfitted_tktk_4vec.Clear();
        tktkRes_4vecs.clear();
        tktkRes_4vec.Clear();
        unfitted_tktkRes_4vec.Clear();
        tktk_candidate.clear();
        tktkRes_candidate.clear();
        tktkCands.clear();
        tktkResCands.clear();
        unfitted_tktk_4vec.SetPxPyPzE(0., 0., 0., 0.);
        unfitted_tktkRes_4vec.SetPxPyPzE(0., 0., 0., 0.);

        //push back the Res tracks as first tracks
        ParticleMass tk_mass;
        //keep track of the push_back track index
        std::vector<int> pushbackTrkIdx;
        std::vector<int> pushbackResTrkIdx;
        float tk_sigma;

///////here have to match the permutatiion with the track mass
		double Mass_in_permutation[16][3] = { \
			   {  0.,              0.,                   0.   }, // 0 (- - -)
			   {  0.,              0.,                   0.   },
			   {  PROTON_MASS,     PION_MASS,            KAON_MASS  }, // 2 (- - +)  perm2+p=2
			   {  PION_MASS,       PROTON_MASS,          KAON_MASS  },//perm2+p=3
			   {  PROTON_MASS,     KAON_MASS,            PION_MASS }, // 4 (- + -)
			   {  PION_MASS,       KAON_MASS,            PROTON_MASS  },
			   {  KAON_MASS,       PROTON_MASS,          PION_MASS }, // 6 (- + +)
			   {  KAON_MASS,       PION_MASS,            PROTON_MASS  },
			   {  KAON_MASS,       PROTON_MASS,          PION_MASS }, // 8 (+ - -)
			   {  KAON_MASS,       PION_MASS,            PROTON_MASS  },
			   {  PROTON_MASS,     KAON_MASS,            PION_MASS }, // 10 (+ - +)
			   {  PION_MASS,       KAON_MASS,            PROTON_MASS  },
			   {  PROTON_MASS,     PION_MASS,            KAON_MASS  }, // 12 (+ + -)
			   {  PION_MASS,       PROTON_MASS,          KAON_MASS  },
			   {  0.,              0.,                   0.   }, // 14 (+ + +)
			   {  0.,              0.,                   0.   }
		};

	


		///new version
		//temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk1].pt(), input_tracks[selectedTkhidxSet[i].index_tk1].eta(), input_tracks[selectedTkhidxSet[i].index_tk1].phi(),fabs(TkMassCharge[0].first));
		temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk1].pt(), input_tracks[selectedTkhidxSet[i].index_tk1].eta(), input_tracks[selectedTkhidxSet[i].index_tk1].phi(),Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0]);
		unfitted_tktk_4vec += temp_vec;
		//temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk2].pt(), input_tracks[selectedTkhidxSet[i].index_tk2].eta(), input_tracks[selectedTkhidxSet[i].index_tk2].phi(),fabs(TkMassCharge[1].first));
		temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk2].pt(), input_tracks[selectedTkhidxSet[i].index_tk2].eta(), input_tracks[selectedTkhidxSet[i].index_tk2].phi(),Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1]);
		unfitted_tktk_4vec += temp_vec;
		//temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk3].pt(), input_tracks[selectedTkhidxSet[i].index_tk3].eta(), input_tracks[selectedTkhidxSet[i].index_tk3].phi(),fabs(TkMassCharge[2].first));
		temp_vec.SetPtEtaPhiM(input_tracks[selectedTkhidxSet[i].index_tk3].pt(), input_tracks[selectedTkhidxSet[i].index_tk3].eta(), input_tracks[selectedTkhidxSet[i].index_tk3].phi(),Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2]);
		unfitted_tktk_4vec += temp_vec;




			
			
/////new version			
			reco::TransientTrack tkTT1(input_tracks[selectedTkhidxSet[i].index_tk1].track(), &(*bField) );
			if (tkTT1.isValid())
			{
				//tk_mass = fabs(TkMassCharge[0].first);
				tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][0];
				tk_sigma = Functs.getParticleSigma(tk_mass);
				tktk_candidate.push_back(pFactory.particle(tkTT1,tk_mass,chi,ndf,tk_sigma));
				pushbackTrkIdx.push_back(selectedTkhidxSet[i].index_tk1);
			}
			reco::TransientTrack tkTT2(input_tracks[selectedTkhidxSet[i].index_tk2].track(), &(*bField) );
			if (tkTT2.isValid())
			{
				//tk_mass = fabs(TkMassCharge[1].first);
				tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][1];
				tk_sigma = Functs.getParticleSigma(tk_mass);
				tktk_candidate.push_back(pFactory.particle(tkTT2,tk_mass,chi,ndf,tk_sigma));
				pushbackTrkIdx.push_back(selectedTkhidxSet[i].index_tk2);
			}
			reco::TransientTrack tkTT3(input_tracks[selectedTkhidxSet[i].index_tk3].track(), &(*bField) );
			if (tkTT3.isValid())
			{
				//tk_mass = fabs(TkMassCharge[2].first);
				tk_mass = Mass_in_permutation[selectedTkhidxSet[i].permutation_number][2];
				tk_sigma = Functs.getParticleSigma(tk_mass);
				tktk_candidate.push_back(pFactory.particle(tkTT3,tk_mass,chi,ndf,tk_sigma));
				pushbackTrkIdx.push_back(selectedTkhidxSet[i].index_tk3);
			}


        if (printFill_Info_) {
			DMassCutLevel[Dchannel_number-1]->Fill(5);
		}

        double MaximumDoca = Functs.getMaxDoca(tktk_candidate);
        if (MaximumDoca > MaxDocaCut_[Dchannel_number-1]) continue;
        if (printFill_Info_) {
			DMassCutLevel[Dchannel_number-1]->Fill(6);
		}

        if(tktkRes_mass>0){
            if(doConstrainFit){
                ParticleMass tktkResMass = tktkRes_mass;
                MultiTrackKinematicConstraint *tktkResConstraint = new TwoTrackMassKinematicConstraint(tktkResMass);
                tktk_VFT = kcv_tktk_fitter.fit(tktk_candidate, tktkResConstraint);
            }
            else tktk_VFT = tktk_fitter.fit(tktk_candidate);
        }
        else{
            tktk_VFT = tktk_fitter.fit(tktk_candidate);
        }

        if(!tktk_VFT->isValid()) continue;
        if (printFill_Info_)  {
			DMassCutLevel[Dchannel_number-1]->Fill(7);
		}

        tktk_VFT->movePointerToTheTop();
        tktk_VFP   = tktk_VFT->currentParticle();
        tktk_VFPvtx = tktk_VFT->currentDecayVertex();
        if (!tktk_VFPvtx->vertexIsValid()) continue;
        if (printFill_Info_)  {
			DMassCutLevel[Dchannel_number-1]->Fill(8);
		}

        double chi2_prob_tktk = TMath::Prob(tktk_VFPvtx->chiSquared(),tktk_VFPvtx->degreesOfFreedom());
        if(chi2_prob_tktk < VtxChiProbCut_[Dchannel_number-1]) continue;
        if (printFill_Info_)  {
			DMassCutLevel[Dchannel_number-1]->Fill(9);
		}

        tktkCands  = tktk_VFT->finalStateParticles();

        if (printFill_Info_)  {
			DMassCutLevel[Dchannel_number-1]->Fill(10);
		}

        for(unsigned int k = 0; k < tktkCands.size(); k++){
            v4_tk.SetPxPyPzE(tktkCands[k]->currentState().kinematicParameters().momentum().x(),
                                     tktkCands[k]->currentState().kinematicParameters().momentum().y(),
                                     tktkCands[k]->currentState().kinematicParameters().momentum().z(),
                                     tktkCands[k]->currentState().kinematicParameters().energy());
            tktk_4vecs.push_back(v4_tk);
            v4_tk.Clear();
        }

        tktk_4vec.SetPxPyPzE(tktk_VFP->currentState().kinematicParameters().momentum().x(),
                tktk_VFP->currentState().kinematicParameters().momentum().y(),
                tktk_VFP->currentState().kinematicParameters().momentum().z(),
                tktk_VFP->currentState().kinematicParameters().energy());


        //fit info
        DInfo.index[DInfo.size]           = DInfo.size;
        DInfo.unfitted_mass[DInfo.size]   = unfitted_tktk_4vec.Mag();
        DInfo.unfitted_pt[DInfo.size]     = unfitted_tktk_4vec.Pt();
        DInfo.mass[DInfo.size]            = tktk_4vec.Mag();
        DInfo.pt[DInfo.size]              = tktk_4vec.Pt();
        DInfo.eta[DInfo.size]             = tktk_4vec.Eta();
        DInfo.phi[DInfo.size]             = tktk_4vec.Phi();
        DInfo.px[DInfo.size]              = tktk_4vec.Px();
        DInfo.py[DInfo.size]              = tktk_4vec.Py();
        DInfo.pz[DInfo.size]              = tktk_4vec.Pz();
        DInfo.MaxDoca[DInfo.size]         = MaximumDoca;

        //for 2D distance
		VertexDistanceXY axy;
		VertexDistance3D a3d;
        //https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_0/RecoVertex/VertexTools/src/VertexDistance3D.cc
        DInfo.svpvDistance[DInfo.size] = a3d.distance(thePrimaryV,tktk_VFPvtx->vertexState()).value();
        DInfo.svpvDisErr[DInfo.size] = a3d.distance(thePrimaryV,tktk_VFPvtx->vertexState()).error();
        if( DInfo.pt[DInfo.size] <= dCutSeparating_PtVal_[Dchannel_number-1] && (DInfo.svpvDistance[DInfo.size]/DInfo.svpvDisErr[DInfo.size]) < svpvDistanceCut_lowptD_[Dchannel_number-1]) continue;
        else if( DInfo.pt[DInfo.size] > dCutSeparating_PtVal_[Dchannel_number-1] && (DInfo.svpvDistance[DInfo.size]/DInfo.svpvDisErr[DInfo.size]) < svpvDistanceCut_highptD_[Dchannel_number-1]) continue;
        if (printFill_Info_)  {
			DMassCutLevel[Dchannel_number-1]->Fill(12);
		}

		////for DCA and its error
		tsos = extrapolator.extrapolate(tktk_VFP->currentState().freeTrajectoryState(),RecoVertex::convertPos(thePrimaryV.position()));
		tsos_2D = transverseExtrapolator.extrapolate(tktk_VFP->currentState().freeTrajectoryState(),RecoVertex::convertPos(theBeamSpotV.position()));//for BS DCA
		Measurement1D cur3DIP;
		Measurement1D cur2DIP;
		GlobalPoint refPoint          = tsos.globalPosition();
		GlobalError refPointErr       = tsos.cartesianError().position();
		GlobalPoint vertexPosition    = RecoVertex::convertPos(thePrimaryV.position());
		GlobalError vertexPositionErr = RecoVertex::convertError(thePrimaryV.error());
		cur3DIP =  (a3d.distance(VertexState(vertexPosition,vertexPositionErr), VertexState(refPoint, refPointErr)));
		GlobalPoint refPoint_2D          = tsos_2D.globalPosition();
		GlobalError refPointErr_2D       = tsos_2D.cartesianError().position();
		GlobalPoint vertexPosition_2D    = RecoVertex::convertPos(theBeamSpotV.position());
		GlobalError vertexPositionErr_2D = RecoVertex::convertError(theBeamSpotV.error());
		cur2DIP = axy.distance(VertexState(refPoint_2D, refPointErr_2D),VertexState(vertexPosition_2D,vertexPositionErr_2D));
		
		std::cout<<"DCA value:  "<<cur3DIP.value()<<std::endl;


        //////
		reco::Vertex::Point vp1(thePrimaryV.position().x(), thePrimaryV.position().y(), 0.);
        reco::Vertex::Point vp2(tktk_VFPvtx->vertexState().position().x(), tktk_VFPvtx->vertexState().position().y(), 0.);
        ROOT::Math::SVector<double, 6> sv1(thePrimaryV.covariance(0,0), thePrimaryV.covariance(0,1), thePrimaryV.covariance(1,1), 0., 0., 0.);
        ROOT::Math::SVector<double, 6> sv2(tktk_VFPvtx->vertexState().error().cxx(), tktk_VFPvtx->vertexState().error().cyx(), tktk_VFPvtx->vertexState().error().cyy(), 0., 0., 0.);
        reco::Vertex::Error ve1(sv1);
        reco::Vertex::Error ve2(sv2);
        reco::Vertex v1(vp1, ve1);
        reco::Vertex v2(vp2, ve2);
        DInfo.svpvDistance_2D[DInfo.size] = a3d.distance(v1, v2).value();
        DInfo.svpvDisErr_2D[DInfo.size] = a3d.distance(v1, v2).error();

		//3D dca and it error
		DInfo.ip3d[DInfo.size]            = cur3DIP.value();
		DInfo.ip3derr[DInfo.size]         = cur3DIP.error();
		DInfo.ip2d_BS[DInfo.size]         = cur2DIP.value();
		DInfo.ip2d_BS_err[DInfo.size]     = cur2DIP.error();

        DInfo.vtxX[DInfo.size]            = tktk_VFPvtx->position().x();
        DInfo.vtxY[DInfo.size]            = tktk_VFPvtx->position().y();
        DInfo.vtxZ[DInfo.size]            = tktk_VFPvtx->position().z();
        DInfo.vtxXErr[DInfo.size]         = tktk_VFPvtx->error().cxx();
        DInfo.vtxYErr[DInfo.size]         = tktk_VFPvtx->error().cyy();
        DInfo.vtxZErr[DInfo.size]         = tktk_VFPvtx->error().czz();
        DInfo.vtxYXErr[DInfo.size]        = tktk_VFPvtx->error().cyx();
        DInfo.vtxZXErr[DInfo.size]        = tktk_VFPvtx->error().czx();
        DInfo.vtxZYErr[DInfo.size]        = tktk_VFPvtx->error().czy();
        DInfo.vtxdof[DInfo.size]          = tktk_VFPvtx->degreesOfFreedom();
        DInfo.vtxchi2[DInfo.size]         = tktk_VFPvtx->chiSquared();
        
        TVector3 svpvVec;
        svpvVec.SetXYZ(DInfo.vtxX[DInfo.size]-EvtInfo.PVx, DInfo.vtxY[DInfo.size]-EvtInfo.PVy, DInfo.vtxZ[DInfo.size]-EvtInfo.PVz);
        TVector3 dVec;
        dVec.SetXYZ(DInfo.px[DInfo.size], DInfo.py[DInfo.size], DInfo.pz[DInfo.size]);
        DInfo.alpha[DInfo.size] = svpvVec.Angle(dVec);
        if( DInfo.alpha[DInfo.size] > alphaCut_[Dchannel_number-1]) continue;
        if (printFill_Info_)  {
			DMassCutLevel[Dchannel_number-1]->Fill(13);
		}
        DInfo.rftk1_mass[DInfo.size]      = tktk_4vecs[0].Mag();
        DInfo.rftk1_pt[DInfo.size]        = tktk_4vecs[0].Pt();
        DInfo.rftk1_eta[DInfo.size]       = tktk_4vecs[0].Eta();
        DInfo.rftk1_phi[DInfo.size]       = tktk_4vecs[0].Phi();
        DInfo.rftk2_mass[DInfo.size]      = tktk_4vecs[1].Mag();
        DInfo.rftk2_pt[DInfo.size]        = tktk_4vecs[1].Pt();
        DInfo.rftk2_eta[DInfo.size]       = tktk_4vecs[1].Eta();
        DInfo.rftk2_phi[DInfo.size]       = tktk_4vecs[1].Phi();
        DInfo.rftk3_mass[DInfo.size]      = tktk_4vecs[2].Mag();
		DInfo.rftk3_pt[DInfo.size]        = tktk_4vecs[2].Pt();
		DInfo.rftk3_eta[DInfo.size]       = tktk_4vecs[2].Eta();
		DInfo.rftk3_phi[DInfo.size]       = tktk_4vecs[2].Phi();

        //Index initialize to -2
        DInfo.rftk1_index[DInfo.size] = -2;
        DInfo.rftk2_index[DInfo.size] = -2;
        DInfo.rftk3_index[DInfo.size] = -2;
        DInfo.rftk4_index[DInfo.size] = -2;
        DInfo.rftk5_index[DInfo.size] = -2;

        DInfo.rftk1_index[DInfo.size]     = pushbackTrkIdx[0];
        DInfo.rftk2_index[DInfo.size]     = pushbackTrkIdx[1];
        if( fabs(tktk_4vecs[0].Mag()-PION_MASS) < fabs(tktk_4vecs[0].Mag()-KAON_MASS) ) 
		{DInfo.rftk1_MassHypo[DInfo.size] = 211;
		}
        else DInfo.rftk1_MassHypo[DInfo.size] = 321;
        if( fabs(tktk_4vecs[0].Mag()-PROTON_MASS) < fabs(tktk_4vecs[0].Mag()-KAON_MASS) && fabs(tktk_4vecs[0].Mag()-PROTON_MASS) < fabs(tktk_4vecs[0].Mag()-PION_MASS) )
		{
			DInfo.rftk1_MassHypo[DInfo.size] = 2212;
        }
		//If its a Res particle, save it as D0
        if( DInfo.rftk1_index[DInfo.size] == -1) DInfo.rftk1_MassHypo[DInfo.size] = 421;
        if( fabs(tktk_4vecs[1].Mag()-PION_MASS) < fabs(tktk_4vecs[1].Mag()-KAON_MASS) ) DInfo.rftk2_MassHypo[DInfo.size] = 211;
        else DInfo.rftk2_MassHypo[DInfo.size] = 321;
        if( fabs(tktk_4vecs[1].Mag()-PROTON_MASS) < fabs(tktk_4vecs[1].Mag()-KAON_MASS) && fabs(tktk_4vecs[1].Mag()-PROTON_MASS) < fabs(tktk_4vecs[1].Mag()-PION_MASS) ) DInfo.rftk2_MassHypo[DInfo.size] = 2212;
		
        if(tktkCands.size()>2){
            DInfo.rftk3_mass[DInfo.size]  = tktk_4vecs[2].Mag();
            DInfo.rftk3_pt[DInfo.size]    = tktk_4vecs[2].Pt();
            DInfo.rftk3_eta[DInfo.size]   = tktk_4vecs[2].Eta();
            DInfo.rftk3_phi[DInfo.size]   = tktk_4vecs[2].Phi();
            DInfo.rftk3_index[DInfo.size] = pushbackTrkIdx[2];
            if( fabs(tktk_4vecs[2].Mag()-PION_MASS) < fabs(tktk_4vecs[2].Mag()-KAON_MASS) ) DInfo.rftk3_MassHypo[DInfo.size] = 211;
            else DInfo.rftk3_MassHypo[DInfo.size] = 321;
            if( fabs(tktk_4vecs[2].Mag()-PROTON_MASS) < fabs(tktk_4vecs[2].Mag()-KAON_MASS) && fabs(tktk_4vecs[2].Mag()-PROTON_MASS) < fabs(tktk_4vecs[2].Mag()-PION_MASS) ) DInfo.rftk3_MassHypo[DInfo.size] = 2212;
        }

        DInfo.type[DInfo.size] = Dchannel_number;
        D_counter[Dchannel_number-1]++;

        tktk_candidate.clear();
        tktkCands.clear();
        DInfo.size++;

    }
	/*
	auto t1 = std::chrono::high_resolution_clock::now(); 
	double dt = 1e-3 * std::chrono::duration_cast<chrono::microseconds>(t1 - t0).count();
	std::cout<<"duration: "<<dt<<std::endl;
*/
}
//}}}

//define this as a plug-in
DEFINE_FWK_MODULE(Dfinder);
